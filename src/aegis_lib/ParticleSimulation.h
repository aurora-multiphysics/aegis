#ifndef ParticleSimulation__
#define ParticleSimulation__

#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>  // for min/max values
#include <set>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <time.h>
#include <any>

#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyLine.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkCompositeDataSet.h>
#include <vtkInformation.h>
#include <vtkSTLReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkAppendFilter.h>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "Inputs.h"
#include "SimpleLogger.h"
#include "EquilData.h"
#include "Source.h"
#include "Integrator.h"
#include "CoordTransform.h"
#include "alglib/interpolation.h"
#include "Particle.h"
#include "VtkInterface.h"
#include "AegisBase.h"


using namespace moab;


enum class meshWriteOptions
{
  FULL, // write out all mesh sets
  TARGET, // write out only the target mesh set
  BOTH, // write out target and full mesh to seperate files
  PARTIAL // write out multiple specified mesh sets to a single file
};

enum class ExecuteOptions
{
  SERIAL,
  MPI,
  MPI_PADDED,
  MPI_DYNAMIC
};

class ParticleSimulation : public AegisBase
{
  public:
  ParticleSimulation(std::shared_ptr<JsonHandler> configFile, std::shared_ptr<EquilData> equil);
  void Execute(); // switch between runs
  void Execute_serial(); // serial
  void Execute_mpi(); // MPI_Gatherv
  void Execute_dynamic_mpi(); // dynamic load balancing
  void Execute_padded_mpi(); // padded MPI_Gather
  void init_geometry();
  std::vector<std::pair<double,double>> psiQ_values; // for l2 norm test
  int target_num_facets();

  protected:

  private:
  void worker();
  void handler(std::vector<double> &handlerQVals);
  void dynamic_task_init();
  void implicit_complement_testing(); 
  moab::Range select_target_surface(); // get target surfaces of interest from aegis_settings.json
  std::vector<double> loop_over_facets(int startFacet, int endFacet); // loop over facets in target surfaces
  
  void setup_sources();

  void cartesian_track();
  void polar_track();
  void flux_track();
  
  terminationState loop_over_particle_track(TriangleSource&tri, std::unique_ptr<ParticleBase> &particle); // loop over individual particle tracks
  void terminate_particle_depositing(TriangleSource&tri); // end particle track
  void terminate_particle_shadow(TriangleSource&tri); // end particle track
  void terminate_particle_lost(TriangleSource&tri); // end particle track
  void terminate_particle_maxlength(TriangleSource&tri); // end particle track
  void test_cyl_ray_fire(std::unique_ptr<ParticleBase> &particle);  
  
  
  void ray_hit_on_launch(std::unique_ptr<ParticleBase> &particle); // particle hit on initial launch from surface
  void print_particle_stats(std::array<int, 5> particleStats); // print number of particles that reached each termination state
  void mpi_particle_stats(); // get inidividual particle stats for each process
  void read_params(const std::shared_ptr<JsonHandler> &inputs); // read parameters from aegis_settings.json
  void attach_mesh_attribute(const std::string &tagName, moab::Range &entities, std::vector<double> &dataToAttach);
  void attach_mesh_attribute(const std::string &tagName, moab::Range &entities, std::vector<std::vector<double>> &dataToAttach);
  
  void write_out_mesh(meshWriteOptions option, moab::Range rangeOfEntities = {});
  void mesh_coord_transform(coordinateSystem coordSys);
  void select_coordinate_system();
  void update_progress_indicator(int batchesComplete, int totalBatches);

  // Simulation config parameters
  std::string executionType = "dynamic";
  std::string dagmcInputFile; // parameter required to be specified by user
  double powerSOL = 0.0;
  double lambdaQ = 0.0; 
  double trackStepSize = 0.001;
  int maxTrackSteps = 0;
  std::string particleLaunchPos = "fixed";
  double userROutrBdry = 0.0;
  int dynamicBatchSize = 16;
  std::string coordinateConfig = "cart";
  bool workerProfiling = false;
  bool workerDebug = false;
  coordinateSystem coordSys = coordinateSystem::CARTESIAN; // default cartesian 
  bool noMidplaneTermination = false;
  bool printMpiParticleStats = false;

  std::vector<TriangleSource> listOfTriangles;
  int totalNumberOfFacets = 0;

  std::vector<int> vectorOfTargetSurfs;
  unsigned int numberOfRayFireCalls = 0;
  unsigned int numberOfClosestLocCalls = 0;
  unsigned int iterationCounter=0;


  // DAGMC variables
  std::unique_ptr<moab::DagMC> DAG;
  moab::Range surfsList; // list of moab::EntityHandle of surfaces in mesh
  moab::Range volsList; // list of moab::EntityHandle of volumes in mesh
  moab::Range facetsList; // list of moab::EntityHandle of facets in mesh
  moab::Range targetFacets; // list of moab::EntityHandle of facets in target mesh
  std::vector<EntityHandle> nodesList; // list of moab::EntityHandle of nodes in mesh
  std::vector<double> nodeCoords; // 1D flattended list of node XYZ coordinates 
  std::vector<double> nodeCoordsPol; // 1D flattended list of node polar coordinates
  int targetNumFacets = 0;
  int numNodes = 0;
  moab::EntityHandle prevSurf;
  moab::EntityHandle nextSurf;
  moab::EntityHandle implComplementVol; 
  double nextSurfDist = 0.0; // distance to next surface
  moab::EntityHandle intersectedFacet;
  int rayOrientation = 1; // rays are fired along surface normals
  double trackLength = 0.0;

  double startTime;
  double facetsLoopTime;

  std::shared_ptr<EquilData> equilibrium;
  
  std::unique_ptr<SurfaceIntegrator> integrator;

  std::unique_ptr<VtkInterface> vtkInterface;
  const std::string branchShadowedPart = "Shadowed Particles";
  const std::string branchLostPart = "Lost Particles";
  const std::string branchDepositingPart = "Depositing Particles";
  const std::string branchMaxLengthPart = "Max Length Particles";
  
};

#endif

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


enum class meshWriteOptions{
FULL, // write out all mesh sets
TARGET, // write out only the target mesh set
BOTH, // write out target and full mesh to seperate files
PARTIAL // write out multiple specified mesh sets to a single file
};

class ParticleSimulation : public AegisBase
{
  public:
  ParticleSimulation(std::string filename);
  void Execute(); // serial
  void Execute_mpi(); // MPI_Gatherv
  void Execute_dynamic_mpi(); // dynamic load balancing
  void Execute_padded_mpi(); // padded MPI_Gather
  void init_geometry();
  int num_facets();
  std::vector<std::pair<double,double>> psiQ_values; // for l2 norm test


  protected:

  private:
  void dynamic_task_init();
  void implicit_complement_testing(); 
  moab::Range select_target_surface(); // get target surfaces of interest from aegis_settings.json
  std::vector<double> loop_over_facets(int startFacet, int endFacet); // loop over facets in target surfaces
  terminationState loop_over_particle_track(const moab::EntityHandle &facet, ParticleBase &particle, DagMC::RayHistory &history); // loop over individual particle tracks
  void terminate_particle(const moab::EntityHandle &facet, DagMC::RayHistory &history, terminationState termination); // end particle track
  void ray_hit_on_launch(ParticleBase &particle, DagMC::RayHistory &history); // particle hit on initial launch from surface
  void print_particle_stats(std::array<int, 5> particleStats); // print number of particles that reached each termination state
  void mpi_particle_stats(); // get inidividual particle stats for each process
  void read_params(const std::shared_ptr<InputJSON> &inputs); // read parameters from aegis_settings.json
  void attach_mesh_attribute(const std::string &tagName, moab::Range &entities, std::vector<double> &dataToAttach);
  void write_out_mesh(meshWriteOptions option, moab::Range rangeOfEntities = {});

  
  int nFacets;
  
  std::string settingsFileName;
  std::shared_ptr<InputJSON> JSONsettings; 
  std::string dagmcInputFile;
  std::string vtkInputFile;
  std::string eqdskInputFile;
  double powerSOL = 0.0;
  double lambdaQ = 0.0; 
  double trackStepSize;
  int maxTrackSteps = 0;
  std::string particleLaunchPos;
  double userROutrBdry;
  std::string drawParticleTracks;
  int dynamicTaskSize;
  double rmove = 0.0;
  double zmove = 0.0;
  double fscale = 1.0;
  double psiref = 0.0;
  bool noMidplaneTermination = false;
  std::vector<int> vectorOfTargetSurfs;
  unsigned int numberOfRayFireCalls = 0;
  unsigned int numberOfClosestLocCalls = 0;
  unsigned int iterationCounter=0;

  std::unique_ptr<moab::DagMC> DAG;
  std::unique_ptr<moab::DagMC> polarDAG;
  moab::Range surfsList;
  moab::Range volsList;
  moab::Range facetsList;
  moab::Range targetFacets;
  int numFacets = 0;
  int numNodes = 0;
  moab::EntityHandle prevSurf;
  moab::EntityHandle nextSurf;
  moab::EntityHandle implComplementVol; 
  double nextSurfDist = 0.0; // distance to next surface
  moab::EntityHandle intersectedFacet;
  int rayOrientation = 1; // rays are fired along surface normals
  double trackLength = 0.0;
  int facetCounter = 0;

  std::vector<double> qValues;
  std::vector<double> psiValues;
  double startTime;
  double endTime;
  double facetsLoopTime;

  EquilData equilibrium;
  bool plotBFieldRZ = false;
  bool plotBFieldXYZ = false;
  
  std::unique_ptr<SurfaceIntegrator> integrator;
  double BdotN = 0.0; // dot product of magnetic field and triangle surface normal
  double Q = 0.0; // heatflux incident on surface triangle
  double psi = 0.0; // value of psi at current position
  double psid = 0.0; // psi - psi_m
  double psiOnSurface = 0.0;

  std::unique_ptr<VtkInterface> vtkInterface;
  const std::string branchShadowedPart = "Shadowed Particles";
  const std::string branchLostPart = "Lost Particles";
  const std::string branchDepositingPart = "Depositing Particles";
  const std::string branchMaxLengthPart = "Max Length Particles";
  
};

#endif

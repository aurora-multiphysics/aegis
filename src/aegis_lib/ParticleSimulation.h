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



using namespace moab;

class ParticleSimulation  
{
  public:
  ParticleSimulation(std::string filename);
  void Execute(); 
  void init_geometry();
  int num_facets();
  std::vector<std::pair<double,double>> psiQ_values; // for l2 norm test


  protected:

  private:
  void dynamic_task_init();
  moab::Range select_target_surface();
  void loop_over_facets(int startFacet, int endFacet,const  moab::Range targetSurfaceList);
  bool loop_over_particle_track(const moab::EntityHandle &facet, ParticleBase &particle, DagMC::RayHistory &history);
  void terminate_particle(const moab::EntityHandle &facet, DagMC::RayHistory &history, terminationState termination);
  void ray_hit_on_launch(ParticleBase &particle, DagMC::RayHistory &history);
  void print_particle_stats(std::array<int, 4> particleStats);
  void mpi_particle_stats();
  int rank, nprocs;
  int nFacets;
  void read_params(const std::shared_ptr<InputJSON> &inputs);
  
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
  double rmove = 0.0;
  double zmove = 0.0;
  double fscale = 1.0;
  double psiref = 0.0;
  bool noMidplaneTermination = false;
  std::vector<int> vectorOfTargetSurfs;

  std::stringstream stringToPrint;

  std::unique_ptr<moab::DagMC> DAG;
  std::unique_ptr<moab::DagMC> polarDAG;
  moab::Range surfsList;
  moab::Range volsList;
  moab::Range facetsList;
  int numFacets = 0;
  int numNodes = 0;
  moab::EntityHandle prevSurf;
  moab::EntityHandle nextSurf;
  moab::EntityHandle volID; 
  double nextSurfDist = 0.0; // distance to next surface
  moab::EntityHandle intersectedFacet;
  int rayOrientation = 1; // rays are fired along surface normals
  double trackLength = 0.0;
  int facetCounter = 0;
  std::vector<double> qValues;
  std::vector<double> psiValues;

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

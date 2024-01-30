#ifndef Aegis__
#define Aegis__

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
#include "settings.hpp"
#include "simpleLogger.h"
#include "equData.h"
#include "source.h"
#include "integrator.h"
#include "coordtfm.h"
#include "alglib/interpolation.h"
#include "particle.h"
#include "vtkAegis.h"



using namespace moab;

class AegisClass  
{
  public:

  void Execute(); 
  void init_solve();
  void init_geometry();
  int num_facets();
  std::vector<std::pair<double,double>> psiQ_values; // for l2 norm test
  moab::Range select_target_surface();
  void max_length_termination(const moab::EntityHandle &facet, DagMC::RayHistory &history);
  void midplane_termination(const moab::EntityHandle &facet, DagMC::RayHistory &history);
  void lost_termination(const moab::EntityHandle &facet, DagMC::RayHistory &history);
  void shadowed_termination(const moab::EntityHandle &facet, DagMC::RayHistory &history);
  void ray_hit_on_launch(particleBase &particle, DagMC::RayHistory &history);
  protected:


  private:
  settings runSettings;
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

  std::unique_ptr<moab::DagMC> DAG;
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

  equData bFieldData;
  bool plotBFieldRZ = false;
  bool plotBFieldXYZ = false;
  
  std::unique_ptr<surfaceIntegrator> integrator;
  double BdotN = 0.0; // dot product of magnetic field and triangle surface normal
  double Q = 0.0; // heatflux incident on surface triangle
  double psi = 0.0; // value of psi at current position
  double psid = 0.0; // psi - psi_m
  bool traceEnded = false;
  double psiOnSurface = 0.0;

  std::unique_ptr<vtkAegis> vtkInterface;
  const std::string branchShadowedPart = "Shadowed Particles";
  const std::string branchLostPart = "Lost Particles";
  const std::string branchDepositingPart = "Depositing Particles";
  const std::string branchMaxLengthPart = "Max Length Particles";

  

};














#endif

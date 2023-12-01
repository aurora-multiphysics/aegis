#ifndef source__
#define source__

#include <stdio.h>
#include <iostream>
#include <cctype>
#include <random>
#include "equData.h"
#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>

//create_source()

// OpenMC has function IndependantSource() that generates an independant source
// It appears to do a number of different things:
// 1) Check particle type 
// 2) Check source strength - relevant if multiple sources present (i.e to pick one source over the other)
// 3) Check for external source file
// 4) Determine the spatial distribution of the source (if None then create point source)
// 5) Determine the angular distribution (if None create isotropic source)
// 6) Determine the source energy disitribution
// 7) Determine the source time distribution (defualt to constant time T=0)

// 7 Things that make up the definition of a radiation source. 

// First attempt could be something like:
// 1) Neutron
// 2) 1.0
// 3) No external Source file
// 4) Spherical dsitribution 
// 5) Isotropic
// 6) OpenMC defualt energy distribution
// 7) Constant time T=0

// Next steps: 
// Create a particleType class 
// Create a source class
// 


using moab::DagMC;
using moab::OrientedBoxTreeTool;
moab::DagMC* DAGInstance;


class pointSource{
  public:
  double r[3];
  double dir[3];

  pointSource(std::vector<double> xyz);
  void set_dir(double newDir[3]);
  void get_isotropic_dir();
  void get_hemisphere_surface_dir(double surfaceNormal[3]);
};


class sphereSource{
  public:
  double sample[3]; //sampled position
  double origin[3]; // origin of sphere
  double r; // radius of sphere
};

class boxSource{
  public:
  double pA[3];
  double pB[3];
  double dir[3];
  double pR[3];
  boxSource(double xyz3[3], double xyz4[3]);
  void get_pt();
  void get_dir();

};

class triSource
{

  private:
  std::vector<double> xyzA;
  std::vector<double> xyzB;
  std::vector<double> xyzC;
  double D;

  public:
  std::vector<double> normal;
  std::vector<double> unitNormal;
  double Bn; // B.n of current triangle 
  moab::EntityHandle entityHandle; // moab EntityHandle of triangle
  triSource(std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyz3, moab::EntityHandle handle);
  void dagmcInstance(moab::DagMC* DAG); 
  std::vector<double> random_pt();
  std::vector<double> centroid();



};

#endif
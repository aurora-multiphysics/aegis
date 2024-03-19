#ifndef Source__
#define Source__

#include <stdio.h>
#include <iostream>
#include <cctype>
#include <random>
#include "EquilData.h"
#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>



using moab::DagMC;
using moab::OrientedBoxTreeTool;


class pointSource
{
  public:
  double r[3];
  double dir[3];

  pointSource(std::vector<double> xyz);
  void set_dir(double newDir[3]);
  void get_isotropic_dir();
  void get_hemisphere_surface_dir(double surfaceNormal[3]);
};


class sphereSource
{ 
  public:
  double sample[3]; //sampled position
  double origin[3]; // origin of sphere
  double r; // radius of sphere
};

class boxSource
{
  public:
  double pA[3];
  double pB[3];
  double dir[3];
  double pR[3];
  boxSource(double xyz3[3], double xyz4[3]);
  void get_pt();
  void get_dir();

};

class TriSource : public AegisBase
{
  private:
  std::vector<double> xyzA;
  std::vector<double> xyzB;
  std::vector<double> xyzC;
  std::vector<double> launchPos;
  double Q = 0.0; // Q of triangle
  double psi = 0.0; // psi at particle start
  double Bn = 0.0; // B.n of Triangle
  moab::EntityHandle entityHandle;
  double D;


  public:
  std::vector<double> normal;
  std::vector<double> unitNormal;

  TriSource(std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyz3,
            moab::EntityHandle handle, std::string launchType);
  
  void dagmcInstance(moab::DagMC* DAG); 
  std::vector<double> random_pt();
  std::vector<double> centroid();
  // set BdotN and Q for the triangle
  void set_heatflux_params(EquilData &equilibrium, const std::string formula);
  void update_heatflux(double newHeatflux);
  std::vector<double> launch_pos(); // return launch position on triangle
  double BdotN();
  double heatflux();
  moab::EntityHandle entity_handle();

};

#endif
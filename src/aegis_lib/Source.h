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


class Sources : public AegisBase
{
  public:
  void set_heatflux_params(const std::shared_ptr<EquilData> &equilibrium, const std::string formula);
  void update_heatflux(double newHeatflux);
  std::vector<double> launch_pos(); // return launch position on triangle
  double BdotN();
  double heatflux();
  void add_heatflux(double heatflux);
  moab::EntityHandle entity_handle();
  double get_psi();
  std::vector<double> get_normal();

  protected:
  std::vector<double> launchPos;
  double Q = 0.0; // Q of triangle
  double _totalHeatflux = 0.0;
  double psi = 0.0; // psi at particle start
  double Bn = 0.0; // B.n of Triangle
  moab::EntityHandle entityHandle;
  std::vector<double> normal;
  std::vector<double> unitNormal;

  private:

};



class TriangleSource : public Sources
{
  private:
  std::vector<double> xyzA;
  std::vector<double> xyzB;
  std::vector<double> xyzC;
  double D;


  public:


  TriangleSource(std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyz3,
            moab::EntityHandle handle, std::string launchType);
  
  void dagmcInstance(moab::DagMC* DAG); 
  std::vector<double> random_pt();
  std::vector<double> centroid();
  // set BdotN and Q for the triangle


};

#endif
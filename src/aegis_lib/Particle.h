#ifndef Particle__
#define Particle__

#include <iostream>
#include <stdio.h>
#include <cctype>
#include <vector>
#include <map>
#include <set>
#include <variant>
#include <unordered_map>
#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>

#include "CoordTransform.h"
#include "EquilData.h"



class ParticleBase 
{
  private:
  double power = 0.0; // power associated with particle
  double lengthTravelled = 0.0; // total length travelled by particle
  bool directionUp = false; 
  std::vector<double> previousPos;
  //DagMC::RayHistory history; // stores entity handles of surfaces crossed and facets hit

  public:
  std::vector<double> Bfield; // magnetic field at current pos
  std::vector<double> dir; // unit direction vector of particle at current position
  std::vector<double> BfieldXYZ; // cartesian mangetic field
  std::vector<double> pos; // current position of particle 
  std::vector<double> launchPos; // initial starting position of particle on triangle
  int atMidplane; // 0 == Not at midplane; 1 == At Inner-Midplane; 2 == At Outer-Midplane  
  bool outOfBounds = false;

  void set_pos(std::vector<double> newPosition); // set new particle position
  void set_pos(double newPosition[]); // overload for C-style array
  std::vector<double> get_pos(std::string coordType); // return STL vector of the current position
  double get_psi(EquilData &EquData); // get psi value at current pos
  void set_dir(EquilData &EquData); // set unit drection vector and Bfield at current position
  std::vector<double> get_dir(std::string coordType); // return STL vector of unit direction
  void check_if_in_bfield(EquilData &Equdata); // check if in magnetic field
  void align_dir_to_surf(double Bn); // align particle dir to surface normal
  void update_vectors(double distanceTravelled); // update position  
  void update_vectors(double distanceTravelled, EquilData &EquData); // overload to update dir as well
  void check_if_midplane_reached(const std::array<double, 3> &midplaneParameters); // check if particle has reached inner or outer midplane and set value of atMidplane

};


#endif
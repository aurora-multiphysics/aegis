#ifndef particle__
#define particle__

#include <iostream>
#include <stdio.h>
#include <cctype>
#include <vector>
#include <map>
#include <set>
#include <variant>
#include <unordered_map>
#include <moab/Core.hpp>

#include "coordtfm.h"
#include "equData.h"



class particleBase 
{
  private:
  double power; // power associated with particle

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
  void set_dir(equData &EquData); // set unit drection vector and Bfield at current position
  std::vector<double> get_dir(std::string coordType); // return STL vector of unit direction
  void check_if_in_bfield(equData &Equdata); // check if in magnetic field
  void align_dir_to_surf(double Bn); // align particle dir to surface normal
  void update_vectors(double distanceTravelled); // update position  
  void update_vectors(double distanceTravelled, equData &EquData); // overload to update dir as well
  void check_if_midplane_reached(double zcen, double rInrBdry, double rOutrBdry); // check if particle has reached inner or outer midplane and set value of atMidplane
};


#endif
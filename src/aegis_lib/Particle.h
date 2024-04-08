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



class ParticleBase : public AegisBase
{
  private:
  double power = 0.0; // power associated with particle
  double lengthTravelled = 0.0; // total length travelled by particle
  bool directionUp = false; 
  std::vector<double> previousPos;
  double thresholdDistanceThreshold = 0.0;
  double euclidDistTravelled = 0.0;
  bool thresholdDistanceCrossed = false;
  coordinateSystem coordSystem = coordinateSystem::CARTESIAN;
  //DagMC::RayHistory history; // stores entity handles of surfaces crossed and facets hit

  public:
  std::vector<double> Bfield; // magnetic field at current pos
  std::vector<double> dir; // unit direction vector of particle at current position
  std::vector<double> BfieldXYZ; // cartesian mangetic field
  std::vector<double> posXYZ; // current position of particle in xyz
  std::vector<double> posPolar;
  std::vector<double> posFlux; 
  std::vector<double> launchPos; // initial starting position of particle on triangle
  int atMidplane; // 0 == Not at midplane; 1 == At Inner-Midplane; 2 == At Outer-Midplane  
  bool outOfBounds = false;
  bool thresholdDistanceSet = false;



  ParticleBase(coordinateSystem coordSys);

  // set new particle position
  void set_pos(std::vector<double> newPosition); 
  
  // overload for C-style array
  void set_pos(double newPosition[]); 

  // return STL vector of the current position
  std::vector<double> get_pos();

  // get psi value at current pos
  double get_psi(const std::shared_ptr<EquilData> &EquData); 

  // set unit drection vector and Bfield at current position
  void set_dir(const std::shared_ptr<EquilData> &EquData); 

  // return STL vector of unit direction
  std::vector<double> get_dir(); 
  void check_if_in_bfield(const std::shared_ptr<EquilData> &Equdata); // check if in magnetic field
  void align_dir_to_surf(double Bn); // align particle dir to surface normal
  void update_vectors(double distanceTravelled); // update position  
  void update_vectors(double distanceTravelled, const std::shared_ptr<EquilData> &EquData); // overload to update dir as well
  void check_if_midplane_crossed(const std::array<double, 3> &midplaneParameters); // check if particle has reached inner or outer midplane and set value of atMidplane
  void set_intersection_threshold(double distanceThreshold);
  bool check_if_threshold_crossed();

};


#endif
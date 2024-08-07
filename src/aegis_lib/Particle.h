#ifndef Particle__
#define Particle__

#include <iostream>
#include <stdio.h>
#include <cctype>
#include <vector>
#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>

#include "CoordTransform.h"
#include "EquilData.h"

enum class magneticFieldDirection
{
  FORWARDS,
  BACKWARDS
};

class ParticleBase : public AegisBase
{
  private:
  Position3D _position;
  std::vector<double> _pos;
  double _heatflux = 0.0; // heatflux associated with particle
  double lengthTravelled = 0.0; // total length travelled by particle
  std::vector<double> previousPos; // previous position in particle history
  moab::DagMC::RayHistory facetHistory; // stores entity handles of surfaces crossed and facets hit
  moab::EntityHandle _parentEntity;
  magneticFieldDirection fieldDir = magneticFieldDirection::BACKWARDS; // set particle along or against magnetic field 
  std::vector<double> _direction; // unit direction vector of particle at current position
  std::vector<double> _launchPos; // initial starting position of particle on triangle
  coordinateSystem _coordSystem = coordinateSystem::CARTESIAN;

  double _rInnerMidplane = 0.0;
  double _rOuterMidplane = 0.0;
  double _zMidplane = 0.0;


  public:
  int atMidplane; // 0 == Not at midplane; 1 == At Inner-Midplane; 2 == At Outer-Midplane  

  // constructors
  ParticleBase(coordinateSystem coordSystem, std::vector<double> startingPosition, double & heatflux, moab::EntityHandle entityHandle, Position3D position);
  ParticleBase(coordinateSystem coordSystem, std::vector<double> startingPosition, double heatflux);  
  ParticleBase(coordinateSystem coordSystem, std::vector<double> startingPosition); 

  void set_facet_history(moab::DagMC::RayHistory history);
  Position3D& posi();
  // set new particle position
  void set_pos(std::vector<double> newPosition); 

  // get psi value at current pos
  double get_psi(const std::shared_ptr<EquilData> &EquData); 

  // set unit drection vector and Bfield at current position
  void set_dir(std::vector<double> &dir); 

  // set coordinate system flag of particle
  void set_coordinate_system(const coordinateSystem &coordSys)
  {
    _coordSystem = coordSys;
  };

  // return coordinate system flag of particle
  coordinateSystem coord_sys()
  {
    return _coordSystem;
  }

  // return STL vector of unit direction
  std::vector<double> dir()
  {
    return _direction;
  };

  // return position vector of particle 
  std::vector<double> pos()
  {
    return _pos;
  };

  // return previous position of particle
  std::vector<double> previous_pos()
  {
    return previousPos;
  };
  
  // return the position the particle was launched from
  std::vector<double> launchPos()
  {
    return _launchPos;
  };
  
  // return heatflux assigned to particle
  double heatflux()
  {
    return _heatflux;
  };

  // return parent entity particle launched from
  moab::EntityHandle
  parent_entity_handle()
  {
    return _parentEntity;
  }


  void align_dir_to_surf(double Bn); // align particle dir to surface normal
  void step(const double &distance); // step particle forwards 
  void update_vectors(const double &distanceTravelled, const std::shared_ptr<EquilData> &EquData); // overload to update dir as well
  void check_if_midplane_crossed(const std::array<double, 3> &midplaneParameters); // check if particle has reached inner or outer midplane and set value of atMidplane
  void set_intersection_threshold(double distanceThreshold);
  bool check_if_threshold_crossed();
  void set_midplane_parameters(const std::array<double, 3> & midplaneParameters);
};


#endif
#include "Particle.h"
#include "SimpleLogger.h"
#include "moab/Interface.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <moab/Core.hpp>
#include <moab/OrientedBoxTreeTool.hpp>

ParticleBase::ParticleBase(coordinateSystem coordSys, std::vector<double> startingPosition,
                           double & heatflux, moab::EntityHandle entityHandle, Position3D position)
{
  _coordSystem = coordSys;
  _pos = startingPosition;
  _launchPos = startingPosition;
  _heatflux = heatflux;
  _parentEntity = entityHandle;
  _position = position;
}

ParticleBase::ParticleBase(coordinateSystem coordSys, std::vector<double> startingPosition,
                           double heatflux)
{
  _coordSystem = coordSys;
  _pos = startingPosition;
  _heatflux = heatflux;
}

void
ParticleBase::set_facet_history(moab::DagMC::RayHistory history)
{
  facetHistory = history;
}
// set current particle position
void
ParticleBase::set_pos(std::vector<double> newPosition)
{
  if (_pos.empty())
  {
    _launchPos = newPosition;
    _pos = _launchPos;
    return;
  }
  previousPos = _pos;
  _pos = newPosition; // set the new position
}

Position3D &
ParticleBase::posi()
{
  return _position;
}

// return double of psi at current pos
double
ParticleBase::get_psi(const std::shared_ptr<EquilData> & equilibrium)
{
  std::vector<double> polarPos;
  double psi;
  switch (_coordSystem)
  {
    case coordinateSystem::CARTESIAN:
      polarPos = CoordTransform::cart_to_polar_xy(_pos);
      psi = equilibrium->get_psi(polarPos[0], polarPos[1]);
      break;
    case coordinateSystem::POLAR:
      psi = equilibrium->get_psi(_pos[0], _pos[1]);
  }
  return psi;
}

// set unit drection vector and Bfield at current position
void
ParticleBase::set_dir(std::vector<double> & dir)
{
  // flip sign depending on surface normal particle launched from
  if (fieldDir == magneticFieldDirection::BACKWARDS)
  {
    dir[0] = -dir[0];
    dir[1] = -dir[1];
    dir[2] = -dir[2];
  }

  _direction = dir;
}

// align particle direction along surface normal of facet particle launched from
void
ParticleBase::align_dir_to_surf(double Bn)
{
  fieldDir = (Bn < 0) ? magneticFieldDirection::BACKWARDS : magneticFieldDirection::FORWARDS;
}

// update position vector of particle with set distance travelled
void
ParticleBase::step(const double & distance)
{
  // This will update the position
  std::vector<double> newPt(3);

  double x2 = pow((newPt[0] - _pos[0]), 2);
  double y2 = pow((newPt[1] - _pos[1]), 2);
  double z2 = pow((newPt[2] - _pos[2]), 2);

  newPt[0] = _pos[0] + _direction[0] * distance;
  newPt[1] = _pos[1] + _direction[1] * distance;
  newPt[2] = _pos[2] + _direction[2] * distance;

  set_pos(newPt);
}

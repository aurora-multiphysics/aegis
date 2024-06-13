
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
                           double heatflux, moab::EntityHandle entityHandle)
{
  coordSystem = coordSys;
  pos = startingPosition;
  launchPos = startingPosition;
  _heatflux = heatflux;
  parentEntity = entityHandle;
}

ParticleBase::ParticleBase(coordinateSystem coordSys, std::vector<double> startingPosition,
                           double heatflux)
{
  coordSystem = coordSys;
  pos = startingPosition;
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
  if (pos.empty())
  {
    launchPos = newPosition;
    pos = launchPos;
    return;
  }
  previousPos = pos;
  pos = newPosition; // set the new position
}

// return STL vector of current particle position
std::vector<double>
ParticleBase::get_pos()
{
  std::vector<double> currentPosition = pos;
  return currentPosition;
}

std::vector<double>
ParticleBase::get_xyz_pos()
{
  switch (coordSystem)
  {
    case coordinateSystem::CARTESIAN:
      return pos;

    case coordinateSystem::POLAR:
      std::vector<double> cartPos = CoordTransform::polar_to_cart(pos);
      return cartPos;
  }
  return pos;
}

// return double of psi at current pos
double
ParticleBase::get_psi(const std::shared_ptr<EquilData> & equilibrium)
{
  std::vector<double> polarPos;
  double psi;
  switch (coordSystem)
  {
    case coordinateSystem::CARTESIAN:
      polarPos = CoordTransform::cart_to_polar_xy(pos);
      psi = equilibrium->get_psi(polarPos[0], polarPos[1]);
      break;
    case coordinateSystem::POLAR:
      psi = equilibrium->get_psi(pos[0], pos[1]);
  }
  return psi;
}

// set unit drection vector and Bfield at current position
void
ParticleBase::set_dir(const std::shared_ptr<EquilData> & equilibrium)
{
  std::vector<double> polarPos; // polar position
  std::vector<double> magn(3);

  switch (coordSystem)
  {
    case coordinateSystem::CARTESIAN:
      polarPos = CoordTransform::cart_to_polar(pos);
      magn = equilibrium->b_field(pos, "cart"); // magnetic field
      magn = equilibrium->b_field_cart(magn, polarPos[2]);
      break;

    case coordinateSystem::POLAR:
      magn = equilibrium->b_field(pos, "polar"); // magnetic field
      break;

    case coordinateSystem::FLUX:
      break;
  }

  // normalise particle direction vector
  double norm = sqrt(pow(magn[0], 2) + pow(magn[1], 2) + pow(magn[2], 2));
  Bfield = magn;
  std::vector<double> normB(3);
  normB[0] = magn[0] / norm;
  normB[1] = magn[1] / norm;
  normB[2] = magn[2] / norm;
  dir = normB;

  // flip sign depending on surface normal particle launched from
  if (fieldDir == magneticFieldDirection::BACKWARDS)
  {
    dir[0] = -dir[0];
    dir[1] = -dir[1];
    dir[2] = -dir[2];
  }
}

// get current particle direction (along magnetic field)
std::vector<double>
ParticleBase::get_dir()
{
  std::vector<double> currentDirection;
  currentDirection = dir;
  return currentDirection;
}

// align particle direction along surface normal of facet particle launched from
void
ParticleBase::align_dir_to_surf(double Bn)
{
  fieldDir = (Bn < 0) ? magneticFieldDirection::BACKWARDS : magneticFieldDirection::FORWARDS;
}

// update position vector of particle with set distance travelled
void
ParticleBase::update_vectors(const double & distanceTravelled)
{
  // This will update the position
  std::vector<double> newPt(3);

  double x2 = pow((newPt[0] - pos[0]), 2);
  double y2 = pow((newPt[1] - pos[1]), 2);
  double z2 = pow((newPt[2] - pos[2]), 2);
  double thresholdDistUpdate = sqrt(x2 + y2 + z2);

  euclidDistTravelled += thresholdDistUpdate;

  newPt[0] = pos[0] + dir[0] * distanceTravelled;
  newPt[1] = pos[1] + dir[1] * distanceTravelled;
  newPt[2] = pos[2] + dir[2] * distanceTravelled;

  set_pos(newPt);
}

// update position vector of particle with set distance travelled and update
// the particle direction at new position
void
ParticleBase::update_vectors(const double & distanceTravelled,
                             const std::shared_ptr<EquilData> & equilibrium)
{
  // This will update the position and direction vector of the particle
  std::vector<double> newPt(3);

  newPt[0] = pos[0] + dir[0] * distanceTravelled;
  newPt[1] = pos[1] + dir[1] * distanceTravelled;
  newPt[2] = pos[2] + dir[2] * distanceTravelled;

  set_pos(newPt);
  set_dir(equilibrium);
}

// check if the particle has reached the outer-midplane of plasma
// (TerminationState = DEPOSITING)
void
ParticleBase::check_if_midplane_crossed(const std::array<double, 3> & midplaneParameters)
{

  double x, y, r;
  double currentZ, prevZ;

  switch (coordSystem)
  {
    case coordinateSystem::CARTESIAN:
      x = pos[0];
      y = pos[1];
      r = sqrt(pow(x, 2) + pow(y, 2));
      currentZ = pos[2];
      prevZ = previousPos[2];
      break;
    case coordinateSystem::POLAR:
      r = pos[0];
      currentZ = pos[1];
      prevZ = previousPos[1];
  }

  atMidplane = 0;

  double rInnerMidplane = midplaneParameters[0]; // R at inner midplane
  double rOuterMidplane = midplaneParameters[1]; // R at outer midplane
  double zMidplane = midplaneParameters[2];      // Z at midplane

  if (prevZ > currentZ)
  {
    if (currentZ <= zMidplane)
    {
      if (r <= rInnerMidplane)
      {
        atMidplane = 1;
      }
      else if (r >= rOuterMidplane)
      {
        atMidplane = 2;
      }
    }
  }
  else if (prevZ < currentZ)
  {
    if (currentZ >= zMidplane)
    {
      if (r <= rInnerMidplane)
      {
        atMidplane = 1;
      }
      else if (r >= rOuterMidplane)
      {
        atMidplane = 2;
      }
    }
  }
}

// set distance threshold where if particle has travlled less, ray_fire() calls
// do not happen
void
ParticleBase::set_intersection_threshold(double distanceThreshold)
{
  thresholdDistanceThreshold = distanceThreshold;
  thresholdDistanceCrossed = false;
  thresholdDistanceSet = true;
}

// check if the distance threshold for ray_fire() calls has been called
// if true returned ray_fire() should be called
// if false returned then particle position should update without ray_fire()
bool
ParticleBase::check_if_threshold_crossed()
{
  if (euclidDistTravelled > thresholdDistanceThreshold)
  {
    thresholdDistanceThreshold = 0.0;
    euclidDistTravelled = 0.0;
    thresholdDistanceSet = false;
    return true;
  }
  else
  {
    return false;
  }
}

double
ParticleBase::heatflux()
{
  return _heatflux;
}

moab::EntityHandle
ParticleBase::parent_entity_handle()
{
  return parentEntity;
}
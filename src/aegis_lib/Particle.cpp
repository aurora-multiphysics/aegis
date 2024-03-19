
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <cctype>
#include "Particle.h"
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "SimpleLogger.h"


ParticleBase::ParticleBase(coordinateSystem coordSys)
{
  coordSystem = coordSys;
}

// set current particle posXYZition
void ParticleBase::set_pos(std::vector<double> newPosition) 
{
  if (posXYZ.empty()) 
  { 
    launchPos = newPosition;
    posXYZ = launchPos;
    return;
  }
  previousPos = posXYZ;
  posXYZ = newPosition; // set the new posXYZition
}

// overload set_pos() for C-style array
void ParticleBase::set_pos(double newPosition[]) 
{
  std::vector<double> tempVector(3);
  for (int i=0; i<3; i++)
  {
    tempVector[i] = newPosition[i];
  }

  if (posXYZ.empty())
  {
    launchPos = tempVector;
  } 
  posXYZ = tempVector; // set the new posXYZition  
}

// Check if particle is currently in magnetic field
void ParticleBase::check_if_in_bfield(EquilData &equilibrium)
{
  std::vector<double> currentBfield = equilibrium.b_field(posXYZ, "cart");
  if (currentBfield.empty()) {outOfBounds = true;}
  else {outOfBounds = false;}
}

// return STL vector of current particle posXYZition
std::vector<double> ParticleBase::get_pos() 
{
  std::vector<double> currentPosition = posXYZ;
  return currentPosition;
}

// return double of psi at current posXYZ
double ParticleBase::get_psi(EquilData &equilibrium) 
{
  std::vector<double> currentPosition;
  double psi;
  currentPosition = CoordTransform::cart_to_polar(posXYZ, "forwards");
  currentPosition = CoordTransform::polar_to_flux(currentPosition, "forwards", equilibrium);
  psi = currentPosition[0];

  return psi;
}

// set unit drection vector and Bfield at current posXYZition
void ParticleBase::set_dir(EquilData &equilibrium) 
{
  check_if_in_bfield(equilibrium);
  if (outOfBounds) {return;} // exit function early if out of bounds

  double norm; // normalisation constant for magnetic field
  std::vector<double> normB(3); // normalised magnetic field
  std::vector<double> polarPos = CoordTransform::cart_to_polar(posXYZ, "forwards"); // polar posXYZition
  std::vector<double> magn(3); 
  magn = equilibrium.b_field(posXYZ, "cart"); // magnetic field

  switch (coordSystem)
  {
  case coordinateSystem::CARTESIAN:
    magn = equilibrium.b_field_cart(magn, polarPos[2], 0);
    Bfield = magn;
    norm = sqrt(pow(magn[0],2) + pow(magn[1],2) + pow(magn[2],2));
    break;
  
  case coordinateSystem::POLAR:
    Bfield = magn;
    norm = sqrt(pow(magn[0],2) + pow(magn[1],2) + pow(magn[2],2));
    break;
  
  case coordinateSystem::FLUX:
    break;
  }

  normB[0] = magn[0]/norm;
  normB[1] = magn[1]/norm;
  normB[2] = magn[2]/norm;
  dir = normB;

}

// get current particle direction (along magnetic field)
std::vector<double> ParticleBase::get_dir()
{
  std::vector<double> currentDirection;
  currentDirection = dir;
  return currentDirection;
}

// align particle direction along surface normal of facet particle launched from
void ParticleBase::align_dir_to_surf(double Bn)
{
  if (Bn < 0)
  {
    dir[0] = -dir[0];
    dir[1] = -dir[1];
    dir[2] = -dir[2];
  }
}

// update posXYZition vector of particle with set distance travelled
void ParticleBase::update_vectors(double distanceTravelled)
{
  // This will update the posXYZition 
  std::vector<double> newPt(3);

  double x2= pow((newPt[0]-posXYZ[0]), 2);
  double y2= pow((newPt[1]-posXYZ[1]), 2);
  double z2= pow((newPt[2]-posXYZ[2]), 2);
  double thresholdDistUpdate = sqrt(x2 + y2 + z2);

  euclidDistTravelled += thresholdDistUpdate;

  newPt[0] = posXYZ[0] + dir[0] * distanceTravelled;
  newPt[1] = posXYZ[1] + dir[1] * distanceTravelled;
  newPt[2] = posXYZ[2] + dir[2] * distanceTravelled;

  set_pos(newPt);
}

// update posXYZition vector of particle with set distance travelled and update the particle direction at new posXYZition
void ParticleBase::update_vectors(double distanceTravelled, EquilData &equilibrium)
{
  // This will update the posXYZition and direction vector of the particle 
  std::vector<double> newPt(3);

  newPt[0] = posXYZ[0] + dir[0] * distanceTravelled;
  newPt[1] = posXYZ[1] + dir[1] * distanceTravelled;
  newPt[2] = posXYZ[2] + dir[2] * distanceTravelled;

  set_pos(newPt);
  set_dir(equilibrium);
}

// check if the particle has reached the outer-midplane of plasma (TerminationState = DEPOSITING)
void ParticleBase::check_if_midplane_crossed(const std::array<double, 3> &midplaneParameters)
{

  double x = posXYZ[0];
  double y = posXYZ[1];
  double currentZ = posXYZ[2];
  double r = sqrt(pow(x,2) + pow(y,2));
  double startZ = launchPos[2];
  double prevZ = previousPos[2]; 
  atMidplane = 0;

  double rInnerMidplane = midplaneParameters[0]; // R at inner midplane
  double rOuterMidplane = midplaneParameters[1]; // R at outer midplane
  double zMidplane = midplaneParameters[2]; // Z at midplane

  if (prevZ > currentZ)
  {
    if (currentZ <= zMidplane)
    {
      if (r <= rInnerMidplane) {atMidplane = 1;}
      else if (r >= rOuterMidplane) {atMidplane = 2;}
    }
  }
  else if (prevZ < currentZ)
  {
    if (currentZ >= zMidplane)
    {
      if (r <= rInnerMidplane) {atMidplane = 1;}
      else if (r >= rOuterMidplane) {atMidplane = 2;}
    }
  }
}

// set distance threshold where if particle has travlled less, ray_fire() calls do not happen
void ParticleBase::set_intersection_threshold(double distanceThreshold)
{
  thresholdDistanceThreshold = distanceThreshold;
  thresholdDistanceCrossed = false;
  thresholdDistanceSet = true;
}

// check if the distance threshold for ray_fire() calls has been called
// if true returned ray_fire() should be called
// if false returned then particle posXYZition should update without ray_fire()
bool ParticleBase::check_if_threshold_crossed()
{
  if (euclidDistTravelled > thresholdDistanceThreshold)
    {
      thresholdDistanceThreshold = 0.0;
      euclidDistTravelled = 0.0;
      thresholdDistanceSet = false;
      return true;
    }
  return false;
}
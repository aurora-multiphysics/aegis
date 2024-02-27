
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <cctype>
#include "Particle.h"
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "SimpleLogger.h"

// set current particle position
void ParticleBase::set_pos(std::vector<double> newPosition) 
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

// overload set_pos() for C-style array
void ParticleBase::set_pos(double newPosition[]) 
{
  std::vector<double> tempVector(3);
  for (int i=0; i<3; i++)
  {
    tempVector[i] = newPosition[i];
  }

  if (pos.empty())
  {
    launchPos = tempVector;
  } 
  pos = tempVector; // set the new position  
}

// Check if particle is currently in magnetic field
void ParticleBase::check_if_in_bfield(EquilData &equilibrium)
{
  std::vector<double> currentBfield = equilibrium.b_field(pos, "cart");
  if (currentBfield.empty()) {outOfBounds = true;}
  else {outOfBounds = false;}
}

// return STL vector of current particle position
std::vector<double> ParticleBase::get_pos(std::string coordType) 
{
  // make coordType lowercase
  std::transform(coordType.begin(), coordType.end(), coordType.begin(), 
                  [](unsigned char c){return std::tolower(c);});

  if (coordType == "pol" || "polar")
  {
    std::vector<double> currentPosition = CoordTransform::cart_to_polar(pos, "forwards");
    return currentPosition;
  }
  else 
  {
    std::vector<double> currentPosition = pos;
    return currentPosition;
  }
}

// return double of psi at current pos
double ParticleBase::get_psi(EquilData &equilibrium) 
{
  std::vector<double> currentPosition;
  double psi;
  currentPosition = CoordTransform::cart_to_polar(pos, "forwards");
  currentPosition = CoordTransform::polar_to_flux(currentPosition, "forwards", equilibrium);
  psi = currentPosition[0];

  return psi;
}

// set unit drection vector and Bfield at current position
void ParticleBase::set_dir(EquilData &equilibrium) 
{
  check_if_in_bfield(equilibrium);
  if (outOfBounds) {return;} // exit function early if out of bounds

  double norm; // normalisation constant for magnetic field
  std::vector<double> normB(3); // normalised magnetic field
  std::vector<double> polarPos = CoordTransform::cart_to_polar(pos, "forwards"); // polar position
  std::vector<double> magn(3); 
  magn = equilibrium.b_field(pos, "cart"); // magnetic field
  Bfield = magn;
  magn = equilibrium.b_field_cart(magn, polarPos[2], 0);
  BfieldXYZ = magn;
  norm = sqrt(pow(magn[0],2) + pow(magn[1],2) + pow(magn[2],2));
  normB[0] = magn[0]/norm;
  normB[1] = magn[1]/norm;
  normB[2] = magn[2]/norm;

  dir = normB;
}

// get current particle direction (along magnetic field)
std::vector<double> ParticleBase::get_dir(std::string coordType)
{
  std::vector<double> currentDirection = dir;
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

// update position vector of particle with set distance travelled
void ParticleBase::update_vectors(double distanceTravelled)
{
  // This will update the position 
  std::vector<double> newPt(3);

  newPt[0] = pos[0] + dir[0] * distanceTravelled;
  newPt[1] = pos[1] + dir[1] * distanceTravelled;
  newPt[2] = pos[2] + dir[2] * distanceTravelled;

  set_pos(newPt);
}

// update position vector of particle with set distance travelled and update the particle direction at new position
void ParticleBase::update_vectors(double distanceTravelled, EquilData &equilibrium)
{
  // This will update the position and direction vector of the particle 
  std::vector<double> newPt(3);

  newPt[0] = pos[0] + dir[0] * distanceTravelled;
  newPt[1] = pos[1] + dir[1] * distanceTravelled;
  newPt[2] = pos[2] + dir[2] * distanceTravelled;

  set_pos(newPt);
  set_dir(equilibrium);
}

// check if the particle has reached the outer-midplane of plasma (TerminationState = DEPOSITING)
void ParticleBase::check_if_midplane_reached(const std::array<double, 3> &midplaneParameters)
{
  double x = pos[0];
  double y = pos[1];
  double currentZ = pos[2];
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
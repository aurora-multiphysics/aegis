
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <cctype>
#include "particle.h"
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "simpleLogger.h"


void particleBase::set_pos(std::vector<double> newPosition) // set current particle position
{
  if (this->pos.empty()) 
  { 
    this->launchPos = newPosition;
  }
  this->pos = newPosition; // set the new position
}

void particleBase::set_pos(double newPosition[]) // overload set_pos() for C-style array
{
  std::vector<double> tempVector(3);
  for (int i=0; i<3; i++)
  {
    tempVector[i] = newPosition[i];
  }

  if (this->pos.empty())
  {
    this->launchPos = tempVector;
  } 
  this->pos = tempVector; // set the new position  
}

void particleBase::check_if_in_bfield(equData &EquData)
{
  std::vector<double> currentBfield = EquData.b_field(this->pos, "cart");
  if (currentBfield.empty()) {this->outOfBounds = true;}
  else {this->outOfBounds = false;}
}

std::vector<double> particleBase::get_pos(std::string coordType) // return STL vector of current particle position
{
  // make coordType lowercase
  std::transform(coordType.begin(), coordType.end(), coordType.begin(), 
                  [](unsigned char c){return std::tolower(c);});

  if (coordType == "pol" || "polar")
  {
    std::vector<double> currentPosition = coordTfm::cart_to_polar(this->pos, "forwards");
    return currentPosition;
  }
  else 
  {
    std::vector<double> currentPosition = this->pos;
    return currentPosition;
  }
}

double particleBase::get_psi(equData &EquData) // return double of psi at current pos
{
  std::vector<double> currentPosition;
  double psi;
  currentPosition = coordTfm::cart_to_polar(this->pos, "forwards");
  currentPosition = coordTfm::polar_to_flux(currentPosition, "forwards", EquData);
  psi = currentPosition[0];

  return psi;
}

void particleBase::set_dir(equData &EquData) // Set unit direction vector along cartesian magnetic field vector
{
  this->check_if_in_bfield(EquData);
  if (this->outOfBounds) {return;} // exit function early if out of bounds

  double norm; // normalisation constant for magnetic field
  std::vector<double> normB(3); // normalised magnetic field
  std::vector<double> polarPos = coordTfm::cart_to_polar(this->pos, "forwards"); // polar position
  std::vector<double> magn(3); 
  magn = EquData.b_field(this->pos, "cart"); // magnetic field
  this->Bfield = magn;
  magn = EquData.b_field_cart(magn, polarPos[2], 0);
  this->BfieldXYZ = magn;
  norm = sqrt(pow(magn[0],2) + pow(magn[1],2) + pow(magn[2],2));
  normB[0] = magn[0]/norm;
  normB[1] = magn[1]/norm;
  normB[2] = magn[2]/norm;

  this->dir = normB;
}

std::vector<double> particleBase::get_dir(std::string coordType)
{
  std::vector<double> currentDirection = this->dir;
  return currentDirection;
}

void particleBase::align_dir_to_surf(double Bn)
{
  if (Bn < 0)
  {
    this->dir[0] = -this->dir[0];
    this->dir[1] = -this->dir[1];
    this->dir[2] = -this->dir[2];
  }
}

void particleBase::update_vectors(double distanceTravelled)
{
  // This will update the position 
  std::vector<double> newPt(3);

  newPt[0] = this->pos[0] + this->dir[0] * distanceTravelled;
  newPt[1] = this->pos[1] + this->dir[1] * distanceTravelled;
  newPt[2] = this->pos[2] + this->dir[2] * distanceTravelled;

  this->set_pos(newPt);
}

void particleBase::update_vectors(double distanceTravelled, equData &EquData)
{
  // This will update the position and direction vector of the particle 
  std::vector<double> newPt(3);

  newPt[0] = this->pos[0] + this->dir[0] * distanceTravelled;
  newPt[1] = this->pos[1] + this->dir[1] * distanceTravelled;
  newPt[2] = this->pos[2] + this->dir[2] * distanceTravelled;

  this->set_pos(newPt);
  this->set_dir(EquData);
}

void particleBase::check_if_midplane_reached(double zcen, double rInrBdry, double rOutrBdry)
{
  double x = this->pos[0];
  double y = this->pos[1];
  double z = this->pos[2];
  double r = sqrt(pow(x,2) + pow(y,2));
  double startZ = this->launchPos[2]; 
  this->atMidplane = 0;

  if (startZ > zcen)
  {
    if (r <= rInrBdry) {this->atMidplane = 1;}
    else if (r >= rOutrBdry) {this->atMidplane = 2;}
  }
  else if (startZ < zcen)
  {
    if (r <= rInrBdry) {this->atMidplane = 1;}
    else if (r >= rOutrBdry) {this->atMidplane = 2;}
  }
}


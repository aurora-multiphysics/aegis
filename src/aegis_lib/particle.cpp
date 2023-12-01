
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
  this->pos = newPosition; // set the new position
}

void particleBase::set_pos(double newPosition[]) // overload set_pos() for C-style array
{
  std::vector<double> tempVector(3); 
  for (int i=0; i<3; i++)
  {
    tempVector[i] = newPosition[i];
  }
  this->pos = tempVector; // set the new position  
}

bool particleBase::check_if_in_bfield(equData &EquData)
{
  std::vector<double> currentBfield = EquData.b_field(this->pos, "cart");
  if (currentBfield.empty()) {return true;}
  else {return false;}
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

void particleBase::set_dir(equData &EquData) // Set unit direction vector along cartesian magnetic field vector
{
  double norm; // normalisation constant for magnetic field
  std::vector<double> normB(3); // normalised magnetic field
  std::vector<double> polarPos = coordTfm::cart_to_polar(this->pos, "forwards"); // polar position
  std::vector<double> magn = EquData.b_field(this->pos, "cart"); // magnetic field
  this->Bfield = magn;
  magn = EquData.b_field_cart(magn, polarPos[2], 0);
  this->BfieldXYZ = magn;
  norm = sqrt(pow(magn[0],2) + pow(magn[1],2) + pow(magn[2],2));
  normB[0] = magn[0]/norm;
  normB[1] = magn[1]/norm;
  normB[2] = magn[2]/norm;

  this->dir = normB;
  this->dirCS[0] = this->dir[0];
  this->dirCS[1] = this->dir[1];
  this->dirCS[2] = this->dir[2];

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

  // update C style array for direction vector

  this->dirCS[0] = this->dir[0];
  this->dirCS[1] = this->dir[1];
  this->dirCS[2] = this->dir[2];
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

void particleBase::update_cs_arrays() 
{
  this->dirCS[0] = this->dir[0];
  this->dirCS[1] = this->dir[1];
  this->dirCS[2] = this->dir[2];

  this->posCS[0] = this->pos[0];
  this->posCS[1] = this->pos[1];
  this->posCS[2] = this->pos[2];
}
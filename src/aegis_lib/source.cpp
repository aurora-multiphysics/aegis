#include <stdio.h>
#include <iostream>
#include <cctype>
#include <random>
#include <math.h>
#include "source.h"

//create_source()

// OpenMC has function IndependantSource() that generates an independant source
// It appears to do a number of different things:
// 1) Check particle type 
// 2) Check source strength - relevant if multiple sources present (i.e to pick one source over the other)
// 3) Check for external source file
// 4) Determine the spatial distribution of the source (if None then create point source)
// 5) Determine the angular distribution (if None create isotropic source)
// 6) Determine the source energy disitribution
// 7) Determine the source time distribution (defualt to constant time T=0)

// 7 Things that make up the definition of a radiation source. 

// First attempt could be something like:
// 1) Neutron
// 2) 1.0
// 3) No external Source file
// 4) Spherical dsitribution 
// 5) Isotropic
// 6) OpenMC defualt energy distribution
// 7) Constant time T=0

// Next steps: 
// Create a particleType class 
// Create a source class
// 


pointSource::pointSource(std::vector<double> xyz)
{
  for (int i=0; i<xyz.size(); i++)
  {
    r[i] = xyz[i];
  }
}

void pointSource::set_dir(double newDir[3])
{
  dir[0] = newDir[0];
  dir[1] = newDir[1];
  dir[2] = newDir[2];
}


// sample a random direction isotropically from point 
void pointSource::get_isotropic_dir()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0,1.0);
  double theta = 2 * M_PI * dist(rng);
  double phi = acos(1 - 2 * dist(rng));
  dir[0] = sin(phi) * cos(theta);
  dir[1] = sin(phi) * sin(theta);
  dir[2] = cos(phi);
}

void pointSource::get_hemisphere_surface_dir(double surfaceNormal[3])
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0,1.0);
  double theta = 2 * M_PI * dist(rng);
  double phi = acos(1 - 2 * dist(rng));
  dir[0] = sin(phi) * cos(theta);
  dir[1] = sin(phi) * sin(theta);
  dir[2] = cos(phi);
}


// boxSource constructor, a square plane can be described by providing position vectors
// with one dimension that is constant 
boxSource::boxSource(double xyz1[3], double xyz2[3])
{
  for (int i=0; i<3; i++)
  {  
    pA[i] = xyz1[i];
    pB[i] = xyz2[i];
  }  
}

void boxSource::get_pt()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0,1);
  double random[3]={dist(rng), dist(rng), dist(rng)};

  for (int i=0; i<3; i++){
    pR[i] = pA[i] + random[i] * (pB[i]-pA[i]);
  }
}

void boxSource::get_dir()
{
  std::vector<double> test(3);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0,1.0);
  double theta = 2 * M_PI * dist(rng);
  double phi = acos(1 - 2 * dist(rng));
  dir[0] = sin(phi) * cos(theta);
  dir[1] = sin(phi) * sin(theta);
  dir[2] = cos(phi);
}



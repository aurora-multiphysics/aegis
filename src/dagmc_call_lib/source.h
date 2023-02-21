#ifndef source__
#define source__

#include <stdio.h>
#include <iostream>
#include <cctype>
#include <random>

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



class Particle{
  public:
  std::string nameType; // i.e neutron
  double weight; 
  double dir[3]; // direction
  double pos[3]; // position
  double energy;
};

class pointSource{
  public:
  double r[3];
  double dir[3];
  Particle neutron;

  pointSource(double xyz[3]);
  void get_isotropic_dir();
};


class sphereSource{
  public:
  double sample[3]; //sampled position
  double origin[3]; // origin of sphere
  Particle neutron;
  double r; // radius of sphere
};

class boxSource{
  public:
  double pA[3];
  double pB[3];
  double dir[3];
  double pR[3];
  Particle neutron;
  boxSource(double xyz3[3], double xyz4[3]);
  void get_pt();
  void get_dir();

};

#endif
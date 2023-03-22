#ifndef integrator__
#define integrator__

#include <stdio.h>
#include <iostream>
#include <cctype>
#include <vector>
#include <map>

#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>

using namespace moab;


using moab::OrientedBoxTreeTool;


// integrator .h file

class surfaceIntegrator 
{
  public:
  int raysTotal=0; // Total number of rays fired into geometry
  int raysHit=0; // Number of rays hit
  int nSurfaces=0; // Number of surfaces in geometry
  std::map<EntityHandle, int> nRays; // Number of rays intersected with a given surface EntityHandle
  
  
  surfaceIntegrator(std::vector<EntityHandle> surface_list); // constructor
  void count_hit(EntityHandle surface_hit); // count hits belonging to each surface

  void ray_reflect_dir(double prev_dir[3], double surface_normal[3], // get reflected dir
                         double reflected_dir[3]); 


};

#endif
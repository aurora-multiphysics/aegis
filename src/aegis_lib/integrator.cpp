#include <stdio.h>
#include <iostream>
#include <cctype>
#include <vector>
#include <map>
#include <math.h>
#include "integrator.h"


#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>

// cpp file for integrator. Count nuber of ray intersections belonging to a given surface

surfaceIntegrator::surfaceIntegrator(std::vector<EntityHandle> surface_list)
{
  nSurfaces = surface_list.size();
  for (auto i:surface_list)
  {
    nRays[i] = 0;
  }
}

void surfaceIntegrator::count_hit(EntityHandle surface_hit)
{
  nRays[surface_hit] +=1;
  raysHit +=1;
}


void surfaceIntegrator::ray_reflect_dir(double prev_dir[3], double surface_normal[3], 
                                        double reflected_dir[3])
{
  double reflect_dot;
  reflect_dot = prev_dir[0]*surface_normal[0] + prev_dir[1]*surface_normal[1] 
                + prev_dir[2]*surface_normal[2]; 
  
  for (int i=0; i<3; i++)
  {
    reflected_dir[i] = prev_dir[i] - 2*reflect_dot*surface_normal[i];
  }
}
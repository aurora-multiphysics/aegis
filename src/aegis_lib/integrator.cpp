
#include "integrator.h"
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>




// surface integrator class constructor 
// initialise list of EntityHandles and maps associated with 
surfaceIntegrator::surfaceIntegrator(moab::Range const &Facets)
{
  LOG_TRACE << "-----surfaceIntegrator CONSTRUCTOR()-----";

  nFacets = Facets.size();
  for (auto i:Facets)
  {
    facetEntities.push_back(i);
    nRays[i] = 0;
    powFac[i] = 0;
  }

}

void surfaceIntegrator::count_hit(EntityHandle facet_hit)
{
  nRays[facet_hit] +=1;
  raysHit +=1;
}


void surfaceIntegrator::ray_reflect_dir(double const prev_dir[3], double const surface_normal[3], 
                                        double reflected_dir[3])
{
  double reflect_dot = 0 ;
  reflect_dot = prev_dir[0]*surface_normal[0] + prev_dir[1]*surface_normal[1] 
                + prev_dir[2]*surface_normal[2]; 
  
  for (int i=0; i<3; i++)
  {
    reflected_dir[i] = prev_dir[i] - 2*reflect_dot*surface_normal[i];
  }
}

int_sorted_map surfaceIntegrator::sort_map(std::unordered_map<moab::EntityHandle, int> const &map)
{
  int_sorted_map sorted_map(map.begin(), map.end());
  return sorted_map;
}

dbl_sorted_map surfaceIntegrator::sort_map(std::unordered_map<moab::EntityHandle, double> const &map)
{
  dbl_sorted_map sorted_map(map.begin(), map.end());
  return sorted_map;
}
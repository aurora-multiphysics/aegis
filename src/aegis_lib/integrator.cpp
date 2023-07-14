


#include "integrator.h"
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "simpleLogger.h"



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

void surfaceIntegrator::count_lost_ray()
{
  raysLost += 1;
}

void surfaceIntegrator::store_heat_flux(EntityHandle facet, double heatflux)
{
  powFac[facet] += heatflux; // tally accumulated heatflux on facet
  if (heatflux > 0 ) {raysHeatDep += 1;} 
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

// return sorted map of EntityHandles against ints 
void surfaceIntegrator::facet_values(std::unordered_map<moab::EntityHandle, int> const &map)
{
    int_sorted_map sortedMap(map.begin(), map.end());
    for (auto const &pair: sortedMap)
    {
      if (pair.second > 0)
      {
        std::cout << "EntityHandle: " << pair.first << "[" << pair.second << "] rays hit" << std::endl;
      }
    }
}

// return sorted map of EntityHandles against doubles 
void surfaceIntegrator::facet_values(std::unordered_map<moab::EntityHandle, double> const &map)
{
    std::cout << std::endl;
    std::cout << "Heat flux deposited per facet" << std::endl;
    std::cout << std::endl;

    dbl_sorted_map sortedMap(map.begin(), map.end());
    for (auto const &pair: sortedMap)
    {
      if (pair.second > 0)
      {
        std::cout << "EntityHandle:" << pair.first << " [" << pair.second << "] MW/m^2" << std::endl;
      }
    }
}

// output values of unordered_map to csv format
void surfaceIntegrator::csv_out(std::unordered_map<moab::EntityHandle, double> const &map)
{
  moab::EntityHandle tri;
  double value;
  double x,y,z;
  std::ofstream heat_out("heat_out.csv");

  for (const auto &i:map) // i -> EntityHandle 
  {
    tri = i.first;
    value = i.second;
    x = launchPositions[tri][0];
    y = launchPositions[tri][1];
    z = launchPositions[tri][2];

    heat_out << x << ", " << y << ", " << z << ", " << value << std::endl; 
  }

}

void surfaceIntegrator::piecewise_multilinear_out(std::unordered_map<moab::EntityHandle, double> const &map)
{
  moab::EntityHandle tri;
  std::vector<double> value; // value associated with EntityHandle
  std::vector<double> x,y,z; // vectors of x,y,z coordinates
  std::ofstream piecewise_multilinear_out("piecewise_multilinear.txt");
  for (const auto &i:map) // i -> EntityHandle 
  {
    tri = i.first;
    value.push_back(i.second);
    x.push_back(launchPositions[tri][0]);
    y.push_back(launchPositions[tri][1]);
    z.push_back(launchPositions[tri][2]);


  }

  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  std::sort(z.begin(), z.end());
  
  piecewise_multilinear_out << "AXIS X" << std::endl;
  for (auto i:x)
  {
    piecewise_multilinear_out << i << " ";
  }
  piecewise_multilinear_out << std::endl << std::endl;

  piecewise_multilinear_out << "AXIS Y" << std::endl;
  for (auto i:y)
  {
    piecewise_multilinear_out << i << " ";
  }
  piecewise_multilinear_out << std::endl << std::endl;

  piecewise_multilinear_out << "AXIS Z" << std::endl;
  for (auto i:z)
  {
    piecewise_multilinear_out << i << " ";
  }
  piecewise_multilinear_out << std::endl << std::endl;

  piecewise_multilinear_out << "DATA" << std::endl;
  for (auto i:value)
  {
    piecewise_multilinear_out << i << std::endl;
  }


}
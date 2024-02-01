#ifndef integrator__
#define integrator__

#include <iostream>
#include <stdio.h>
#include <cctype>
#include <vector>
#include <map>
#include <set>
#include <variant>
#include <unordered_map>
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "equData.h"



using namespace moab;


using moab::OrientedBoxTreeTool;

struct comp
{
    template<typename T>
    bool operator()(const T &l, const T &r) const
    {
        if (l.second != r.second) {
            return l.second > r.second;
        }
        return l.first > r.first;
    }
};

typedef std::set<std::pair<EntityHandle, int>, comp> int_sorted_map;
typedef std::set<std::pair<EntityHandle, double>, comp> dbl_sorted_map;


// integrator .h file

class surfaceIntegrator 
{
  public:

  int nParticlesTotal=0; // Total number of rays fired into geometry
  int nParticlesShadowed=0; // Number of particles shadowed
  int nParticlesLost=0; // number of particles lost from magnetic domain
  int nParticlesHeatDep=0; // number of particles despositing heat
  int nParticlesMaxLength=0; // number of particles reached max termination length

  int nFacets=0; // Number of facets in geometry
  std::vector<EntityHandle> facetEntities; // list of all entity handles in geometry
  std::unordered_map<EntityHandle, int> particlesShadowedBy; // Number of rays intersected with a given surface EntityHandle
  std::unordered_map<EntityHandle, double> powFac; // power assigned to each facet
  std::unordered_map<EntityHandle, std::vector<double>> launchPositions; // launch positions on each facet

  // Methods
  surfaceIntegrator(moab::Range const &Facets); // constructor (with moab::Range)
  surfaceIntegrator(std::vector<EntityHandle> const &Facets); // constructor (with std::vector<EntityHandle>)

  void count_hit(EntityHandle facet_hit); // count hits belonging to each facet
  void count_lost_ray();
  void store_heat_flux(EntityHandle facet, double heatflux); // store the power associated with a particular facet
  void count_particle(EntityHandle facet, int terminationState, double heatflux);

  void ray_reflect_dir(double const prev_dir[3], double const surface_normal[3], // get reflected dir
                         double reflected_dir[3]); 
  
  // return sorted map of entityhandles against chosen value
  void facet_values(std::unordered_map<EntityHandle, int> const &map);
  void facet_values(std::unordered_map<EntityHandle, double> const &map);

  void csv_out(std::unordered_map<moab::EntityHandle, double> const &map);
  void piecewise_multilinear_out(std::unordered_map<moab::EntityHandle, double> const &map);

  void print_particle_stats(); // return number of particles depositing, shadowed, lost, etc...


/// TODO - Create a class that holds the unordered map and a string for the name of the map
///      - Potentially also any other variables that may be useful (i.e units for power values) 

};

#endif
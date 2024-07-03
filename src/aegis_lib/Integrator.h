#ifndef Integrator__
#define Integrator__

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
#include "EquilData.h"
#include "AegisBase.h"
#include "Source.h"



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

// terminate particle depending on 1 of 4 termination states:
//
// DEPOSITING - Particle reaches outer midplane and deposits power on facet.
// Heatflux = Q SHADOWED - Particle hits another piece of geometry. Heatflux = 0
// LOST - Particle leaves magnetic field so trace stops. Heatflux = 0
// MAX LENGTH - Particle reaches maximum user set length before anything else.
// Heatflux = 0
enum class terminationState{
DEPOSITING, // midplane reached 
SHADOWED, // shadow geometry hit
LOST, // left magnetic domain
MAXLENGTH, // max length reached
PADDED // padded null particle
};


// integrator .h file

class SurfaceIntegrator : public AegisBase
{
  public:
  SurfaceIntegrator();
  SurfaceIntegrator(moab::Range const &Facets); // constructor (with moab::Range)
  SurfaceIntegrator(std::vector<EntityHandle> const &Facets); // constructor (with std::vector<EntityHandle>)
  
//   void q_values(); // return list of qvalues for each triangle
//   void psi_values(); // return list of psi values for each triangle
  void set_facets(moab::Range &Facets);
  void count_hit(EntityHandle facet_hit); // count hits belonging to each facet
  void count_lost_ray();
  void store_heat_flux(EntityHandle facet, double heatflux); // store the power associated with a particular facet
  void count_particle(const EntityHandle &facet, terminationState termination,const double &heatflux);

  void ray_reflect_dir(double const prev_dir[3], double const surface_normal[3], // get reflected dir
                         double reflected_dir[3]); 
  
  // return sorted map of entityhandles against chosen value
  void facet_values(std::unordered_map<EntityHandle, int> const &map);
  void facet_values(std::unordered_map<EntityHandle, double> const &map);

  void csv_out(std::unordered_map<moab::EntityHandle, double> const &map);
  void piecewise_multilinear_out(std::unordered_map<moab::EntityHandle, double> const &map);

  void print_particle_stats(); // return number of particles depositing, shadowed, lost, etc...

  void set_launch_position(const moab::EntityHandle &facet, const std::vector<double> &position);
  std::array<int, 4> particle_stats(); // return array of integers corresponding to particle stats
  void clear_stats();


  private:
  std::vector<TriangleSource> listOfElements; // list of all elements in aegis run

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

};

#endif
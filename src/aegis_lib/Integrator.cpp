

#include "Integrator.h"
#include <moab/Core.hpp>
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "SimpleLogger.h"

SurfaceIntegrator::SurfaceIntegrator() {}

// surface integrator class constructor
// initialise list of EntityHandles and maps associated with
SurfaceIntegrator::SurfaceIntegrator(moab::Range const & Facets)
{

  nFacets = Facets.size();

  for (auto i : Facets)
  {
    facetEntities.push_back(i);
    particlesShadowedBy[i] = 0;
    powFac[i] = 0;
  }
}

void
SurfaceIntegrator::set_facets(moab::Range & Facets)
{
  nFacets = Facets.size();
  for (const auto i : Facets)
  {
    facetEntities.push_back(i);
    particlesShadowedBy[i] = 0;
    powFac[i] = 0;
  }
}

// Overload for constructor use with STL vector
SurfaceIntegrator::SurfaceIntegrator(std::vector<EntityHandle> const & Facets)
{

  nFacets = Facets.size();

  for (auto i : Facets)
  {
    facetEntities.push_back(i);
    particlesShadowedBy[i] = 0;
    powFac[i] = 0;
  }
}

void
SurfaceIntegrator::count_hit(EntityHandle facet_hit)
{
  particlesShadowedBy[facet_hit] += 1;
  nParticlesShadowed += 1;
}

void
SurfaceIntegrator::count_lost_ray()
{
  nParticlesLost += 1;
}

void
SurfaceIntegrator::store_heat_flux(EntityHandle facet, double heatflux)
{
  powFac[facet] += heatflux; // tally accumulated heatflux on facet
  if (heatflux > 0)
  {
    nParticlesHeatDep += 1;
  }
}

void
SurfaceIntegrator::count_particle(const EntityHandle & facet, terminationState termination,
                                  const double & heatflux)
{

  if (heatflux < 0.0)
  {
    throw std::invalid_argument("Heatflux cannot be a negative value");
  }

  // count number of particles that reach each state and store the associated heatflux
  switch (termination)
  {
    case terminationState::DEPOSITING:
      nParticlesHeatDep++;
      break;

    case terminationState::SHADOWED:
      nParticlesShadowed++;
      if (heatflux > 0.0)
      {
        throw std::invalid_argument("Heatflux should be 0 for SHADOWED particles");
      }
      break;

    case terminationState::LOST:
      nParticlesLost++;
      if (heatflux > 0.0)
      {
        throw std::invalid_argument("Heatflux should be 0 for LOST particles");
      }
      break;

    case terminationState::MAXLENGTH:
      nParticlesMaxLength++;
      if (heatflux > 0.0)
      {
        throw std::invalid_argument(
            "Heatflux should be 0 for particles that have reached MAX LENGTH");
      }
      break;

    default:
      throw std::invalid_argument("Invalid particle termination state provided.");
  }

  if (facet != 0)
  {
    powFac[facet] += heatflux;
  }
}

void
SurfaceIntegrator::ray_reflect_dir(double const prev_dir[3], double const surface_normal[3],
                                   double reflected_dir[3])
{
  double reflectDot = 0;
  reflectDot = prev_dir[0] * surface_normal[0] + prev_dir[1] * surface_normal[1] +
               prev_dir[2] * surface_normal[2];

  for (int i = 0; i < 3; i++)
  {
    reflected_dir[i] = prev_dir[i] - 2 * reflectDot * surface_normal[i];
  }
}

// return sorted map of EntityHandles against ints
void
SurfaceIntegrator::facet_values(std::unordered_map<moab::EntityHandle, int> const & map)
{
  if (rank == 0)
  {
    int_sorted_map sortedMap(map.begin(), map.end());
    for (auto const & pair : sortedMap)
    {
      if (pair.second > 0)
      {
        std::cout << "EntityHandle: " << pair.first << "[" << pair.second << "] rays hit"
                  << std::endl;
      }
    }
  }
}

// return sorted map of EntityHandles against doubles
void
SurfaceIntegrator::facet_values(std::unordered_map<moab::EntityHandle, double> const & map)
{
  if (rank == 0)
  {
    std::cout << std::endl;
    std::cout << "Heat flux deposited per facet" << std::endl;
    std::cout << std::endl;

    dbl_sorted_map sortedMap(map.begin(), map.end());
    for (auto const & pair : sortedMap)
    {
      if (pair.second > 0)
      {
        std::cout << "EntityHandle:" << pair.first << " [" << pair.second << "] MW/m^2"
                  << std::endl;
      }
    }
  }
}

// output values of unordered_map to csv format
void
SurfaceIntegrator::csv_out(std::unordered_map<moab::EntityHandle, double> const & map)
{
  moab::EntityHandle tri;
  double value;
  double x, y, z;
  std::ofstream heatOut("heat_out.txt");
  heatOut << "x,"
          << "y,"
          << "z,"
          << "heat" << std::endl;

  for (const auto & i : map) // i -> EntityHandle
  {
    tri = i.first;
    value = i.second;
    x = launchPositions[tri][0];
    y = launchPositions[tri][1];
    z = launchPositions[tri][2];

    heatOut << x << "," << y << "," << z << "," << value << std::endl;
  }
}

void
SurfaceIntegrator::piecewise_multilinear_out(
    std::unordered_map<moab::EntityHandle, double> const & map)
{
  moab::EntityHandle tri;
  std::vector<double> value;   // value associated with EntityHandle
  std::vector<double> x, y, z; // vectors of x,y,z coordinates
  std::ofstream piecewiseMultilinearOut("piecewise_multilinear.txt");
  for (const auto & i : map) // i -> EntityHandle
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

  piecewiseMultilinearOut << "AXIS X" << std::endl;
  for (auto i : x)
  {
    piecewiseMultilinearOut << i << " ";
  }
  piecewiseMultilinearOut << std::endl << std::endl;

  piecewiseMultilinearOut << "AXIS Y" << std::endl;
  for (auto i : y)
  {
    piecewiseMultilinearOut << i << " ";
  }
  piecewiseMultilinearOut << std::endl << std::endl;

  piecewiseMultilinearOut << "AXIS Z" << std::endl;
  for (auto i : z)
  {
    piecewiseMultilinearOut << i << " ";
  }
  piecewiseMultilinearOut << std::endl << std::endl;

  piecewiseMultilinearOut << "DATA" << std::endl;
  for (auto i : value)
  {
    piecewiseMultilinearOut << i << std::endl;
  }
}

// return number of particles depositing, shadowed, lost etc.
void
SurfaceIntegrator::print_particle_stats()
{
  std::cout << "DEPOSITING - " << nParticlesHeatDep << "\n";
  std::cout << "SHADOWED - " << nParticlesShadowed << "\n";
  std::cout << "LOST - " << nParticlesLost << "\n";
  std::cout << "MAX LENGTH - " << nParticlesTotal << "\n";
  int totalParticlesHandled =
      nParticlesHeatDep + nParticlesShadowed + nParticlesLost + nParticlesMaxLength;
  std::cout << "TOTAL - " << totalParticlesHandled << "\n";
};

void
SurfaceIntegrator::set_launch_position(const moab::EntityHandle & facet,
                                       const std::vector<double> & position)
{
  launchPositions[facet] = position;
}

std::array<int, 4>
SurfaceIntegrator::particle_stats()
{

  std::array<int, 4> particleStats;
  particleStats[0] = nParticlesHeatDep;
  particleStats[1] = nParticlesShadowed;
  particleStats[2] = nParticlesLost;
  particleStats[3] = nParticlesMaxLength;

  return particleStats;
}

void
SurfaceIntegrator::clear_stats()
{
  nParticlesHeatDep = 0;
  nParticlesShadowed = 0;
  nParticlesLost = 0;
  nParticlesShadowed = 0;
  nParticlesTotal = 0;
}
#include "ParticleSimulation.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <array>
#include <mpi.h>
#include <unistd.h>
#include "Integrator.h"

int
main(int argc, char ** argv)
{
  // Setup MPI
  int rank, nprocs;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  double startTime = MPI_Wtime(); // get start time across all processes
  // ---------------------------------------------------------------------------------

  // Get Config file
  std::string configFileName;
  if (argc > 1)
  {
    configFileName = argv[1];
  }
  else
  {
    if (rank == 0)
    {
      std::cout << "Error - No config file provided, defaulting to "
                   "'aegis_settings.json'"
                << std::endl;
    }
    configFileName = "aegis_settings.json";
  }
  // ---------------------------------------------------------------------------------

  // equilibrium setup
  auto configFile = std::make_shared<JsonHandler>(configFileName);
  double equilibriumInstantiationStart = MPI_Wtime();
  auto equilibrium = std::make_shared<EquilData>(configFile);
  equilibrium->move();
  equilibrium->psiref_override();
  equilibrium->init_interp_splines();
  equilibrium->centre(1);
  double equilibriumInstantionTime = MPI_Wtime();
  // ---------------------------------------------------------------------------------

  // Integrator setup
  auto integrator = std::make_shared<SurfaceIntegrator>();
  // ---------------------------------------------------------------------------------

  // Main Particle simulation
  ParticleSimulation simulation(configFile, equilibrium, integrator);
  simulation.Execute();
  // ---------------------------------------------------------------------------------

  // Print individual process stats
  for (int i = 1; i < nprocs; ++i)
  {
    if (rank == i)
    {
      std::cout << "\nProcess " << i << " handled particles: \n";
      integrator->print_particle_stats();
      double endTime = MPI_Wtime();
      double totalTime = endTime - startTime;
      std::cout << "Elapsed wall Time on process " << i << " = " << totalTime << "\n";
    }
  }
  std::cout << "\n";
  // ---------------------------------------------------------------------------------

  // Print wall time and other profiling times
  if (rank == 0)
  {
    std::map<std::string, double> profilingTimes = simulation.get_profiling_times();
    profilingTimes.insert(
        std::make_pair("Equilibrium Instantiation Time = ", equilibriumInstantionTime));
    std::cout << "------------------------------------ \n";
    for (auto const & [key, value] : profilingTimes)
    {
      std::cout << key << value << "\n";
    }
    std::cout << "Total wall time = " << MPI_Wtime() - startTime << "\n";
    std::cout << "------------------------------------" << std::endl;
  }
  // ---------------------------------------------------------------------------------

  MPI_Finalize();
  return 0;
}

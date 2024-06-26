#include "ParticleSimulation.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <mpi.h>
#include <unistd.h>

int
main(int argc, char ** argv)
{

  int rank, nprocs;
  MPI_Init(&argc, &argv);

  double startTime = MPI_Wtime(); // get start time on across all processes

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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

  auto configFile = std::make_shared<JsonHandler>(configFileName);

  double equilibriumInstantiationStart = MPI_Wtime();
  auto equilibrium = std::make_shared<EquilData>(configFile);
  equilibrium->move();
  equilibrium->psiref_override();
  equilibrium->init_interp_splines();
  equilibrium->centre(1);
  double equilibriumInstantionTime = MPI_Wtime();

  ParticleSimulation simulation(configFile, equilibrium);
  simulation.Execute();

  // print wall times for each process
  for (int i = 1; i < nprocs; ++i)
  {
    if (rank == i)
    {
      double endTime = MPI_Wtime();
      double totalTime = endTime - startTime;
      std::cout << "Elapsed wall Time on process " << i << " = " << totalTime << std::endl;
      std::cout << "----------------------------" << std::endl << std::endl;
    }
  }

  MPI_Finalize();
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
  return 0;
}

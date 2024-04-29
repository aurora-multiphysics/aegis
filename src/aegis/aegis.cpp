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

  double startTime = MPI_Wtime();
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

  auto equilibrium = std::make_shared<EquilData>(configFile);
  equilibrium->move();
  equilibrium->psiref_override();
  equilibrium->init_interp_splines();
  equilibrium->centre(1);
  equilibrium->write_bfield();

  ParticleSimulation simulation(configFile, equilibrium);
  simulation.Execute();

  auto mpiTimings = configFile->data()["global_params"]["mpi_timings"];

  if (mpiTimings) // print individual timings for each MPI rank
  {
    for (int i = 1; i < nprocs; ++i)
    {
      if (rank == i)
      {
        double endTime = MPI_Wtime();
        double totalTime = endTime - startTime;
        std::cout << "Elapsed wall Time on process " << i << " = " << totalTime << "\n";
        std::cout << "---------------------------- \n \n";
      }
    }
  }

  if (rank == 0) // print final wall time once rank 0 is complete
  {
    double endTime = MPI_Wtime();
    double totalTime = endTime - startTime;
    std::cout << "\n---------------------------- \n";
    std::cout << "Total elapsed wall Time = " << totalTime << "\n";
    std::cout << "---------------------------- \n" << std::endl;
  }

  MPI_Finalize();

  return 0;
}

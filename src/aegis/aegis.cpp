#include "ParticleSimulation.h"
#include <cstdio>
#include <cstring>
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

  auto configFile = std::make_shared<InputJSON>(configFileName);

  auto equilibrium = std::make_shared<EquilData>();
  equilibrium->setup(configFile);
  equilibrium->move();
  equilibrium->psiref_override();
  equilibrium->init_interp_splines();
  equilibrium->centre(1);
  equilibrium->write_bfield();

  ParticleSimulation simulation(configFile, equilibrium);
  simulation.Execute();

  for (int i = 0; i < nprocs; ++i)
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

  return 0;
}

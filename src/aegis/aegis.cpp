#include "ParticleSimulation.h"
#include <mpi.h>

int main(int argc, char **argv) {

  int rank, nprocs;
  MPI_Init(&argc, &argv);

  double startTime = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::string settingsFile;
  if (argc > 1) 
  {
    settingsFile = argv[1];
  }
  else 
  {
    if (rank == 0) {std::cout << "Error - No config file provided, defaulting to 'aegis_settings.json'" << std::endl;}
    settingsFile = "aegis_settings.json";
  }

  ParticleSimulation simulation(settingsFile);
  
  if (nprocs < 2)
  {
    std::cout << "Running on a single core in serial mode... \n";
    simulation.Execute();
  }
  else
  {
    simulation.Execute_mpi();
  }


  for (int i=0; i<nprocs; ++i){
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



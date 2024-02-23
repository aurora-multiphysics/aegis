#include "ParticleSimulation.h"
#include <mpi.h>

int main(int argc, char **argv) {

  int rank, nprocs;
  MPI_Init(&argc, &argv);

  double startTime = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (nprocs < 2)
  {
      std::cout << "Error! This program requires at least two processes to run. Exiting... \n";
      std::exit(1);
  }

  std::string settingsFile;
  if (argc > 1) 
  {
    settingsFile = argv[1];
  }
  else 
  {
    std::cout << "Error - No config file provided, defaulting to 'aegis_settings.json'" << std::endl;
    settingsFile = "aegis_settings.json";
  }
 
  ParticleSimulation simulation(settingsFile);
  simulation.Execute();

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



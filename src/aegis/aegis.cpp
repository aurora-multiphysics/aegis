#include "aegisClass.h"
#include <mpi.h>

int main(int argc, char **argv) {

  int rank, nprocs;

  MPI_Init(&argc, &argv);
  double startTime = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  AegisClass aegis;
  aegis.Execute("runSettings.txt");




  double endTime = MPI_Wtime();
  double totalTime = endTime - startTime;

  if (rank == 0)
  {
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Elapsed W Time = " << totalTime << " on number of processors: " << nprocs << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
  }
 

  MPI_Finalize();

  return 0;
}



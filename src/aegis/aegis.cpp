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


  for (int i=0; i<nprocs; ++i){
    if (rank == i)
    {
      double endTime = MPI_Wtime();
      double totalTime = endTime - startTime;
      std::cout << "Elapsed wall Time on process " << i << " = " << totalTime << std::endl; 
    }
  }

  MPI_Finalize();

  return 0;
}



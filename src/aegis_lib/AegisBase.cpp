#include "AegisBase.h"
#include <mpi.h>




// print error message and MPI_Finalize() exit out of code
void AegisBase::error_exit_mpi(std::string inString, std::string file, std::string func, int line, int rank)
{
  std::cout << file << ":" << func << ":" << line << ":RANK[" << rank << "] --- " 
            << inString << std::endl;
  if (MPI_Finalize() != MPI_SUCCESS) {
    std::cout << "ERROR: MPI_Finalize != MPI_SUCCESS" << std::endl;
  }
  std::exit(1);
}

// cout for strings only on rank 0
void AegisBase::print_mpi(std::string inString)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    std::cout << inString << std::endl;
  }
}

void AegisBase::setup_profiling()
{
  profiling.open("aegis_profiling.txt");
  profiling << "Object to time | " << " elapsed wall time | " << "Rank of process" << std::endl;
  log_string(LogLevel::WARNING, "saving aegis profiling data to aegis_profiling.txt");
}

void AegisBase::out_mpi_wtime(std::string inString, double totalTime)
{
  profiling << inString << " | " << totalTime << "s | " << rank << std::endl;
}

// wrapper for boost LOG macros to test for MPI rank
void AegisBase::log_string(LogLevel level, std::string inString)
{

  if (rank == 0)
  {
    switch (level)
    {
      case LogLevel::INFO:
        LOG_INFO << inString; 
        break;

      case LogLevel::TRACE:
        LOG_TRACE << inString; 
        break;

      case LogLevel::WARNING:
        LOG_WARNING << inString; 
        break;

      case LogLevel::ERROR:
        LOG_ERROR << inString; 
        break;

      case LogLevel::FATAL:
        LOG_FATAL << inString; 
        break;

      default:
        LOG_WARNING << inString;
    }
  }

}

// return MPI rank
int AegisBase::get_mpi_rank()
{
  return rank; 
}

// return MPI size in MPI_COMM_WORLD
int AegisBase::get_mpi_size()
{
  return nprocs;
}

// call MPI_Comm_rank and MPI_Comm_size 
void AegisBase::set_mpi_params()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
}


void AegisBase::setup_logger()
{

}

void AegisBase::string_to_lowercase(std::string &inputStr)
{
  std::transform(inputStr.begin(), inputStr.end(), inputStr.begin(), 
                [](unsigned char c){return std::tolower(c);});
}
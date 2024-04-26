#include "AegisBase.h"
#include <mpi.h>

// print error message and MPI_Abort() exit out of code
void
AegisBase::error_abort_mpi()
{
  std::cerr << "Terminating program... \n";
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void
AegisBase::debug_error_exit(std::string errorString, char * file, char * func, int line, int rank)
{
  std::cout << file << ":" << func << ":" << line << ":RANK[" << rank << "] --- " << errorString
            << std::endl;
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

// cout for strings only on rank 0
void
AegisBase::print_mpi(std::string inString)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    std::cout << inString << std::endl;
  }
}

// wrapper for boost LOG macros that only prints on rank 0
void
AegisBase::log_string(LogLevel level, std::string inString)
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
  else
  {
    return;
  }
}

// return MPI rank
int
AegisBase::get_mpi_rank()
{
  return rank;
}

// return MPI size in MPI_COMM_WORLD
int
AegisBase::get_mpi_size()
{
  return nprocs;
}

// call MPI_Comm_rank and MPI_Comm_size
void
AegisBase::set_mpi_params()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
}

void
AegisBase::string_to_lowercase(std::string & inputStr)
{
  std::transform(inputStr.begin(), inputStr.end(), inputStr.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}
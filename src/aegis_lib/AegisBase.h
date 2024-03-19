#ifndef AegisBase__
#define AegisBase__

#include <iostream>
#include <string>
#include <fstream>
#include <SimpleLogger.h>
#include <boost/log/trivial.hpp>
#include <boost/log/sources/global_logger_storage.hpp>

enum class LogLevel
{
  INFO, 
  TRACE,
  WARNING,
  ERROR,
  FATAL
};

enum class coordinateSystem
{
  CARTESIAN, // cartesian (x,y,z)
  POLAR, // cylindrical polar (r,z,phi)
  FLUX // flux (psi,theta,phi)
};

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;

// Base class for all AEGIS objects with some utility functions (mostly around MPI)
class AegisBase 
{
  public:

  protected:
  void error_exit_mpi(std::string inString, std::string file, std::string func, int line, int rank);
  void print_mpi(std::string inString);
  void log_string(LogLevel level, std::string inString);
  void set_mpi_params();
  int get_mpi_rank();
  int get_mpi_size();
  void string_to_lowercase(std::string &inputStr);
  
  src::severity_logger_mt<boost::log::trivial::severity_level> logger;

  int rank; 
  int nprocs;

  private:


};

#endif
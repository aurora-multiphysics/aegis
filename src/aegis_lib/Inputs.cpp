#include "Inputs.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <sstream>

JsonHandler::JsonHandler() = default;

JsonHandler::JsonHandler(std::string filename)
{
  filepath = filename;

  read_json();
}

JsonHandler::JsonHandler(nlohmann::json existingJSON) { jsonData = existingJSON; }

nlohmann::json
JsonHandler::data()
{
  return jsonData;
}

void
JsonHandler::read_json()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ifstream file(filepath);
  if (file.is_open())
  {
    if (rank == 0)
    {
      std::cout << "Settings found in JSON file '" << filepath << "'" << std::endl;
    }
    jsonData = nlohmann::json::parse(file);
    return;
  }
  else
  {
    if (rank == 0)
    {
      std::cout << "Error! The requested JSON file '" << filepath << "' does not exist."
                << std::endl;
    }
    std::exit(EXIT_FAILURE);
  }
}

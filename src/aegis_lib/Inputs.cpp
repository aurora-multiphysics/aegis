#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include "Inputs.h"
#include <mpi.h>



InputJSON::InputJSON()
{

}

InputJSON::InputJSON(std::string filename){
  filepath = filename;

  data = read_json();
}

json InputJSON::read_json(){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ifstream file(filepath);
  if (file.is_open())
  {
    if (rank == 0) {std::cout << "Settings found in JSON file '" << filepath << "'" << std::endl;}
    json data = json::parse(file); 
    return data;
  }
  else
  {
    if (rank == 0) {std::cout << "Error! The requested JSON file '" << filepath << "' does not exist." << std::endl;}  
    std::exit(EXIT_FAILURE);
  }
} 
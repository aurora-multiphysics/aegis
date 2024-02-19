#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include "inputs.h"



InputJSON::InputJSON()
{

}

InputJSON::InputJSON(std::string filename){
  filepath = filename;

  data = read_json();
}

json InputJSON::read_json(){
  std::ifstream file(filepath);
  if (file.is_open())
  {
    std::cout << "Settings found in JSON file '" << filepath << "'" << std::endl;
    json data = json::parse(file); 
    return data;
  }
  else
  {
    std::cout << "Error! The requested JSON file '" << filepath << "' does not exist." << std::endl;  
    exit(EXIT_FAILURE);
  }
} 
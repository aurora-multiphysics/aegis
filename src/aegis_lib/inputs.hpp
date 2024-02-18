#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <stdio.h>
#include <vector>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <any>
#include <nlohmann/json.hpp>

/**
 * Class for storing input settings from JSON file
 * More structured format to replace settings class
*/

using json = nlohmann::json;

class InputJSON{
  public:
  InputJSON();
  InputJSON(std::string filename);
  json read_json();
  json data;

  private:
  std::string filepath; 

};

#endif

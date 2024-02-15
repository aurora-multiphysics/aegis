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


// Object for storing user input settings 

class settings{
  public: 
    std::string ray_qry;
    std::string dagmc_input;
    std::string vtk_input;
    std::string runcase;
    std::string eqdsk_file;
    std::string trace;
    std::unordered_map<std::string, std::string> sValues;
    std::unordered_map<std::string, int> iValues;
    std::unordered_map<std::string, double> dValues;
    std::unordered_map<std::string, std::vector<int>> vValues;
    std::unordered_map<std::string, std::any> settingsVals;
    std::vector<int> surfs;
    void load_params(std::string fileName);
    void print_params();
    void load_any_params();


};

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

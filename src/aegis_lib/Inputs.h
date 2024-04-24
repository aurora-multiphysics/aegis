#ifndef inputs__
#define inputs__

#include <stdio.h>
#include <vector>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <any>
#include <nlohmann/json.hpp>
#include <variant>

/**
 * Class for storing input settings from JSON file
 * More structured format to replace settings class
*/

using json = nlohmann::json;

class JsonHandler
{
  public:
  JsonHandler();
  JsonHandler(std::string filename);
  JsonHandler(json existingJSON);
  json data();
  void read_json();

  template <class T> std::optional<T> get(std::string parameter)
  {
    if (jsonData.contains(parameter))
    {
      return (jsonData[parameter]);
    }
    else 
    {
      return std::nullopt;
    }
  }
  

  private:
  std::string filepath; 
  json jsonData;

};

#endif

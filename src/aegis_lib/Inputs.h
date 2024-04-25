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
#include <AegisBase.h>

/**
 * Class for storing input settings from JSON file
 * More structured format to replace settings class
*/

using json = nlohmann::json;

class JsonHandler : public AegisBase
{
public:

  JsonHandler();
  JsonHandler(std::string filename);
  JsonHandler(json existingJSON);
  json data();
  void read_json();

  template <class T> std::optional<T> get_optional(std::string paramName)  
  {
    if (jsonData.contains(paramName))
    {
      return (jsonData[paramName]);
    }
    else 
    {
      std::stringstream ss;
      ss << "Default value used for parameter " << paramName;
      log_string(LogLevel::INFO, ss.str());
      return std::nullopt;
    }  
  }
  
private:
  std::string filepath; 
  json jsonData;

};

#endif

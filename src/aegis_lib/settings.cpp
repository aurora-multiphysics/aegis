#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include "settings.hpp"

// Object for storing 

void settings::load_params()
{
  std::ifstream settingsFile("runSettings.txt");
  std::string lineStr;
  while (getline(settingsFile, lineStr, '\n'))
  {
    std::istringstream lineSS(lineStr);
    std::string word, type, paramName;
    lineSS >> word; // pull parameter name from line
    paramName = word;
    lineSS >> type;

    if (type == "[s]")
    {
      lineSS >> word;
      sValues[paramName] = word;
    }
    else if (type == "[i]")
    {
      lineSS >> word;
      iValues[paramName] = std::stoi(word);
    }
    else if (type == "[d]")
    {
      lineSS >> word;
      dValues[paramName] = std::stod(word);
    }
    else if (type == "[v]")
    {
      while (lineSS >> word)
      {
        vValues[paramName].push_back(std::stoi(word));
      }
    }
  }
}

void settings::print_params()
{
   
   std::cout << std::endl;
   std::cout << "----- String Valued Parameters -----" << std::endl;
  for (auto const &pair: sValues) {
    std::cout << "{ " << pair.first << ": " << pair.second << "}\n";
  }
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "----- Int Valued Parameters -----" << std::endl;
  for (auto const &pair: iValues) {
    std::cout << "{ " << pair.first << ": " << pair.second << "}\n";
  }
  std::cout << std::endl;
  
  std::cout << std::endl;
  std::cout << "----- Double Valued Parameters -----" << std::endl;
  for (auto const &pair: dValues) {
    std::cout << "{ " << pair.first << ": " << pair.second << "}\n";
  }
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "----- Vector<int> Valued Parameters -----" << std::endl;
  std::vector<int> temp;
  std::string ttemp;
  for(auto const &pair: vValues)
  { 
    ttemp = pair.first;
    temp = pair.second;
  }
  std::cout << "{ " << ttemp << ": ";
  for (auto i:temp)
  {
    std::cout << i << " ";
  }
  std::cout << "}" << std::endl;
  std::cout << std::endl;

}
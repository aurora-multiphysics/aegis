#include <stdio.h>
#ifndef SETTINGS_HPP
#define SETTINGS_HPP
// Object for storing user input settings 

class settings{
  public: 
    std::string ray_qry;
    std::string geo_input;
    std::string runcase;
    std::string eqdsk_file;
    std::string source_power;
    std::string reflections;
    std::string nSample;
    int cenopt;
    void load_settings();
};

#endif

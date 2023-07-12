#include <stdio.h>
#ifndef SETTINGS_HPP
#define SETTINGS_HPP
// Object for storing user input settings 

class settings{
  public: 
    std::string ray_qry;
    std::string dagmc_input;
    std::string vtk_input;
    std::string runcase;
    std::string eqdsk_file;
    std::string trace;
    double Psol;
    double lambda_q;
    std::string reflections;
    int nSample;
    int cenopt;
    int surf;
    void load_settings();
};

#endif

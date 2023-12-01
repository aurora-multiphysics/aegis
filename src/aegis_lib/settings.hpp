#include <stdio.h>
#include <vector>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
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
    std::unordered_map<std::string, std::string> sValues;
    std::unordered_map<std::string, int> iValues;
    std::unordered_map<std::string, double> dValues;
    std::unordered_map<std::string, std::vector<int>> vValues;
    double Psol;
    double lambda_q;
    std::string reflections;
    int nSample;
    int cenopt;
    int surf;
    std::vector<int> surfs;
    int nTrack;
    double dsTrack;
    double psiref; // reference psi value for psibdry
    void load_settings(); // DEPRECIATED
    void load_params();
    void print_params();
};

#endif

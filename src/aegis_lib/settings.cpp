#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include "settings.hpp"

// Object for storing 


    void settings::load_settings()
    {
      std::ifstream in("settings.txt");
      if (!in.is_open())
      {
        std::cout << "Cannot open settings file" << std::endl;
        return;
      }

      int valueInt;

      std::string param;
      std::string valueStr;

      while (!in.eof())
      {
        in >> param;
        in >> valueStr;

        if (param == "DAGMC_input")
        {
          dagmc_input = valueStr; // string
        }
        if (param == "VTK_input")
        {
          vtk_input = valueStr;
        }
        else if (param == "ray_qry")
        {
          ray_qry = valueStr; // string
        }
        else if (param == "eqdsk_file")
        {
          eqdsk_file = valueStr; // string
        }
        else if (param == "runcase")
        {
          runcase = valueStr; // string
        }
        else if (param == "trace")
        {
          trace = valueStr;
        }
	      else if (param == "Psol")
	      {
	        Psol = std::stod(valueStr); // double
	      }
        else if (param == "lambda_q")
        {
          lambda_q = std::stod(valueStr);
        }
	      else if (param == "reflections")
	      {
	      reflections = valueStr; // int
	      }
        else if (param == "number_of_rays_launched_per_tri")
	      {
	        nSample = std::stoi(valueStr); // int
	      }
        else if (param == "cenopt")
        {
          cenopt = std::stoi(valueStr);
        }
        else if (param == "surf")
        {
          surf = std::stoi(valueStr);
        }
      }
      in.close();
    }


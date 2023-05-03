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

        if (param == "geo_input")
        {
          geo_input = valueStr; // string
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
	      else if (param == "source_power")
	      {
	      source_power = valueStr; // double
	      }
	      else if (param == "reflections")
	      {
	      reflections = valueStr; // int
	      }
        else if (param == "number_of_rays_fired")
	      {
	      nSample = valueStr; // int
	      }
        else if (param == "cenopt")
        cenopt = std::stoi(valueStr);
      }
      in.close();
    }


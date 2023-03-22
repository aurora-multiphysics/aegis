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

      std::string param;
      std::string value;

      while (!in.eof())
      {
        in >> param;
        in >> value;

        if (param == "geo_input")
        {
          geo_input = value; // string
        }
        else if (param == "ray_qry")
        {
          ray_qry = value; // string
        }
        else if (param == "eqdsk_file")
        {
          eqdsk_file = value; // string
        }
        else if (param == "runcase")
        {
          runcase = value; // string
        }
	      else if (param == "source_power")
	      {
	      source_power = value; // double
	      }
	      else if (param == "reflections")
	      {
	      reflections = value; // int
	      }
        else if (param == "number_of_rays_fired")
	      {
	      nSample = value; // int
	      }
      }
      in.close();
    }


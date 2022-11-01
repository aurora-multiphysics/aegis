#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <limits>  // for min/max values
#include <set>
#include <vector>
//#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>



using namespace moab;

using moab::DagMC;
using moab::OrientedBoxTreeTool;
moab::DagMC* DAG;

static const char input_file[] = "triple_block.h5m";

int main() {
  DAG = new DagMC();
  DAG->load_file(input_file); // open test dag file 
  DAG->init_OBBTree();
  DagMC::RayHistory history;
  int vol_idx = 1;
  int n_surfs = DAG->num_entities(2);
  int n_vols = DAG->num_entities(3);
  std::cout << "No. of surfaces - " << n_surfs << std::endl;
  std::cout << "No. of Volumes - " << n_vols << std::endl;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 1.0, 0.0};       // ray along x direction
  double origin[3] = {-5.0, 0.0, 0.0};  // origin at 0 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  int invol;
  int nrayfire = 1;
  double ray_start[3];
  EntityHandle vol_n;
  double angle[3];
  EntityHandle* nearest_surf;
  double result;

  std::ofstream ray_coords;
  ray_coords.open("ray_coords.txt");

  ray_coords << origin[0] << ' ' << origin[1] << ' ' << origin[2] << ' ' << std::endl; 

  // check if in volume before ray_fire
  DAG->point_in_volume(vol_h, origin, invol); 
  std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
  std::cout << "Volume ID - " << vol_h << std::endl;
  history.reset();


  while (nrayfire < 10){

    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      origin[i] = origin[i] + (next_surf_dist * dir[i]); 
      ray_coords << origin[i] << ' ';
    }
    ray_coords << std::endl;
    std::cout << std::endl;

    std::cout << "Firing ray" << std::endl;
    DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
    DAG->point_in_volume(vol_h, origin, invol); 
    std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
    std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
    std::cout << "Next Surface id - " << next_surf << std::endl;
    DAG->get_angle(next_surf, NULL, angle, &history);
    std::cout << "Angle calculated from get_angle - ";
    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      std::cout << angle[i] << " ";
    }
    if (next_surf == 0) {
      DAG->next_vol(next_surf, vol_h, vol_h);
      //std::cout << "Ran next_vol" <<std::endl;
      dir[1]=-dir[1]; 
      }
    std::cout << "Next Volume ID - " << vol_h << std::endl; 
    nrayfire += 1;
    
  }

  std::cout << "No. of ray_fire calls - " << nrayfire << std::endl;
    return 0;
}

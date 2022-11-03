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

static const char input_file[] = "big_pipe.h5m";

int main() {
  DAG = new DagMC();
  DAG->load_file(input_file); // open test dag file
  DAG->init_OBBTree(); // initialise OBBTree 
  DagMC::RayHistory history; // initialise RayHistory object
  int vol_idx = 1; 
  int n_surfs = DAG->num_entities(2); // No. of surfaces in loaded geometry
  int n_vols = DAG->num_entities(3); // No. of volumes in loaded geometry
  std::cout << "No. of surfaces - " << n_surfs << std::endl;
  std::cout << "No. of Volumes - " << n_vols << std::endl;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 1.0};       // ray launch direction
  double origin[3] = {-6.0, 0.0, -6.0};  // ray launch origin
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; // previous next_surf_dist before the next ray_fire call
  EntityHandle next_surf; // surface id 
  int invol; // in volume result
  int nrayfire=0; // No. of ray_fire calls 
  EntityHandle prev_surf; // previous surface id
  int n_sample=1000;  
  double sample_step_length;
  double sample_point[3];


  std::ofstream ray_coords; // define stream ray_coords
  ray_coords.open("ray_coords.txt"); // write out stream ray_coords to "ray_coords.txt" file 




  // check if in volume before ray_fire
  DAG->point_in_volume(vol_h, origin, invol);
  std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
  std::cout << "Volume ID - " << vol_h << std::endl;
  history.reset(); // reset history before launching any rays

  // launch
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


  sample_step_length = next_surf_dist/n_sample;
  std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));

  for (int i=0; i<=n_sample; i++){
    ray_coords << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
    for (int j=0; j<3; ++j) { // calculate next sample point along ray
      sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    }
  }

  nrayfire +=1;
  for (int i=0; i<3; ++i) { // calculate next ray launch point
    origin[i] = origin[i] + (next_surf_dist * dir[i]);
  }
    std::cout << std::endl;
    std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
    std::cout << "Next Surface id - " << next_surf << std::endl;
    DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
    DAG->point_in_volume(vol_h, origin, invol);
    std::cout << "Volume ID - " << vol_h << std::endl;
    std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;

    // launch
    DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


  sample_step_length = next_surf_dist/n_sample;
  std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;

  for (int i=0; i<=n_sample; i++){
    ray_coords << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
    for (int j=0; j<3; ++j) { // calculate next sample point along ray
      sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    }
  }

    nrayfire +=1;
    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      origin[i] = origin[i] + (next_surf_dist * dir[i]);
    }
    std::cout << std::endl;
    DAG->point_in_volume(vol_h, origin, invol);
    std::cout << "Volume ID - " << vol_h << std::endl;
    std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
    std::cout << "Next Surface id - " << next_surf << std::endl;
    DAG->point_in_volume(vol_h, origin, invol);
    std::cout << "Volume ID - " << vol_h << std::endl;
    std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;

  while (next_surf !=0){
    dir[0]=-dir[0]; // flip Z-component to keep ray in pipe
    prev_surf = next_surf;
    prev_surf_dist = next_surf_dist;
    DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);

  std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;
    nrayfire +=1;
    if (next_surf_dist < 1e-3){ // fix distance if too small 
      next_surf_dist = prev_surf_dist;
    }  
  sample_step_length = next_surf_dist/n_sample;
  for (int i=0; i<=n_sample; i++){
    ray_coords << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
    for (int j=0; j<3; ++j) { // calculate next sample point along ray
      sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    }
  }

    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      origin[i] = origin[i] + (next_surf_dist * dir[i]);
    }
    std::cout << std::endl;
    std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
    std::cout << "Next Surface id - " << next_surf << std::endl;
    if (prev_surf == next_surf){
      DAG->next_vol(next_surf, vol_h, vol_h);
      DAG->point_in_volume(vol_h, origin, invol);
      std::cout << "Volume ID - " << vol_h << std::endl;
      std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
    }

  }
  std::cout << "No. of ray_fire calls - " << nrayfire << std::endl; // print the No. of ray_fire calls
    return 0;
}

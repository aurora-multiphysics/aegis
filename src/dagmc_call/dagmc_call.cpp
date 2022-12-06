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
#include <sstream>
#include <limits>  // for min/max values
#include <set>
#include <vector>
#include <array>
//#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>


using namespace moab;

using moab::DagMC;
using moab::OrientedBoxTreeTool;
moab::DagMC* DAG;

void next_pt(double prev_pt[3], double origin[3], double next_surf_dist, 
                          double dir[3], std::ofstream &ray_intersect);
double dot_product(double vector_a[], double vector_b[]);
void reflect(double dir[3], double prev_dir[3], EntityHandle next_surf);
double * vecNorm(double vector[3]);


static const char input_file[] = "sduct.h5m";
static const char ray_qry_exps[] = "exps00000200.qry";
static const char ray_qry_unif[] = "unif0000100.qry";


int main() {
  DAG = new DagMC();
  DAG->load_file(input_file); // open test dag file
  DAG->init_OBBTree(); // initialise OBBTree 
  DagMC::RayHistory history; // initialise RayHistory object
  int vol_idx = 1; 
  int n_surfs = DAG->num_entities(2); // No. of surfaces in loaded geometry
  int n_vols = DAG->num_entities(3); // No. of volumes in loaded geometry
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {-1.0, 1.0, 3.0};       // ray launch direction
  double origin[3] = {0.0, -15.0, -25.0};  // ray launch origin
  double prev_dir[3];
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; // previous next_surf_dist before the next ray_fire call
  EntityHandle next_surf; // surface id 
  int invol; // point in volume result
  int nrayfire=0; // No. of ray_fire calls 
  EntityHandle prev_surf; // previous surface id
  int n_sample=1000;  // number of samples along ray 
  double sample_step_length; // step length between sample positions along ray
  double sample_point[3]; // position of sample point along ray
  double dir_mag; // magnitude of ray direction vector
  double normal[3]; // vector of surface facet normal intersected by last ray_fire
  double reflect_dot; // dot product of ray dir vector and normal vector  


  std::cout << "No. of surfaces - " << n_surfs << std::endl;
  std::cout << "No. of Volumes - " << n_vols << std::endl;

  std::vector<std::vector<double>> rayqry; // xyz data from rayqry file

  // Read in qry data 
  std::ifstream ray_input(ray_qry_exps);
  std::string line;
  double word;

  if (ray_input) {
        while(getline(ray_input, line, '\n'))        
        {
            //create a temporary vector that will contain all the columns
            std::vector<double> tempVec;
            std::istringstream ss(line);
            //read word by word(or int by int) 
            while(ss >> word)
            {
                //std::cout<<"word:"<<word<<std::endl;
                //add the word to the temporary vector 
                tempVec.push_back(word);
            }             
            //now all the words from the current line has been added to the temporary vector 
            rayqry.emplace_back(tempVec);
        }    
    }
    else 
    {
        std::cout<<"file cannot be opened"<<std::endl;
    }
    ray_input.close();

    //now you can do the whatever processing you want on the vector
    int j=0;
    int k;
    int qrymax = rayqry.size();
    double dir_array[rayqry.size()][3]; 
    std::cout << rayqry.size() << " - rayqry.size()" << std::endl;
    std::cout << qrymax << " - qrymax" << std::endl;

    double raydirs[qrymax][3];
    //next_dir(raydir, qrymax, dir_array);
     for (int j=0; j<qrymax; j+=2)
     {
      for (int k=0; k<3; k++)
      {
        raydirs[j][k] = dir_array[j+1][k] - dir_array[j][k]; 
        std::cout << raydirs[j][k]<<std::endl;
      }
      std::cout << std::endl;
    }
    // define qrydata origin
    double qryorigin[3];
    double intersect[3];
    qryorigin[0] = dir_array[0][0];
    qryorigin[1] = dir_array[0][1];
    qryorigin[2] = dir_array[0][2];
    std::ofstream ray_intersect; // stream for ray-surface intersection points
    ray_intersect.open("ray_intersect.txt"); // write out stream to "ray_intersect.txt" file
    ray_intersect << std::setprecision(5) << std::fixed;
    ray_intersect << qryorigin[0] << ' ' << qryorigin[1] << ' ' << qryorigin[2] << ' ' << std::endl; // store first ray origin

    // check if in volume before ray_fire
    DAG->point_in_volume(vol_h, qryorigin, invol);
    std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
    std::cout << "Volume ID - " << vol_h << std::endl;
    double *normdir;
    //loop through qrydata
    for (int qrycount=0; qrycount < qrymax; qrycount+=2){
      normdir = vecNorm(raydirs[qrycount]);
      std::cout << normdir[0] << std::endl;
      std::cout << normdir[1] << std::endl;
      std::cout << normdir[2] << std::endl;
      std::cout << std::endl;
      history.reset(); // reset history before launching ray
      // launch
      DAG->ray_fire(vol_h, qryorigin, normdir, next_surf, next_surf_dist, &history, 0, 1);

      // calculate next intersection point and write out to textfile
      for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
        intersect[i] = qryorigin[i] + (next_surf_dist * normdir[i]);
        ray_intersect << intersect[i] << ' ';
        }
        ray_intersect << std::endl;


    }
    







  // dir_mag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  // std::cout << "---- DIR_mag = " << dir_mag << " ----" << std::endl;

  // for (int i=0; i<3; ++i) { // calculate next ray launch point
  //   dir[i] = dir[i]/dir_mag;
  //   std::cout << "---- DIR[" << i << "] = " << dir[i] << " ----" << std::endl;
  // }

  // std::ofstream ray_coords1; // define stream ray_coords1
  // std::ofstream ray_coords2; // define stream ray_coords2
  // std::ofstream ray_intersect; // stream for ray-surface intersection points
  // ray_coords1.open("ray_coords1.txt"); // write out stream ray_coords1 to "ray_coords1.txt" file 
  // ray_coords2.open("ray_coords2.txt"); // write out stream ray_coords2 to "ray_coords1.txt" file 
  // ray_intersect.open("ray_intersect.txt"); // write out stream to "ray_intersect.txt" file

  // ray_intersect << origin[0] << ' ' << origin[1] << ' ' << origin[2] << ' ' << std::endl; // store first ray origin



  // // check if in volume before ray_fire
  // DAG->point_in_volume(vol_h, origin, invol);
  // std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
  // std::cout << "Volume ID - " << vol_h << std::endl;
  // history.reset(); // reset history before launching any rays

  // // launch
  // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


  // sample_step_length = next_surf_dist/n_sample;
  // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));

  // for (int i=0; i<=n_sample; i++){
  //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
  //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //   }
  // }

  // nrayfire +=1;

  // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
  // std::cout << std::endl;
  // std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
  // std::cout << "Next Surface id - " << next_surf << std::endl;
  // DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
  // DAG->point_in_volume(vol_h, origin, invol);
  // std::cout << "Volume ID - " << vol_h << std::endl;
  // std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;

  // //launch
  // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


  // sample_step_length = next_surf_dist/n_sample;
  // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;

  // for (int i=0; i<=n_sample; i++){
  //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
  //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //   }
  // }

  //   nrayfire +=1;

  // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);    
  // std::cout << std::endl;
  // DAG->point_in_volume(vol_h, origin, invol);
  // std::cout << "Volume ID - " << vol_h << std::endl;
  // std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
  // std::cout << "Next Surface id - " << next_surf << std::endl;
  // DAG->point_in_volume(vol_h, origin, invol);
  // std::cout << "Volume ID - " << vol_h << std::endl;
  // std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
  
  
  // while (next_surf !=0){
  //   prev_dir[0] = dir[0];
  //   prev_dir[1] = dir[1];
  //   prev_dir[2] = dir[2];

  //   DAG->get_angle(next_surf, NULL, normal, &history);
  //   reflect_dot = dir[0]*normal[0] + dir[1]*normal[1] + dir[2]*normal[2];
  //   for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
  //     dir[i] = dir[i] - 2*reflect_dot*normal[i];
  //   }

  //   prev_surf = next_surf;
  //   prev_surf_dist = next_surf_dist;

  //   if (next_surf_dist < 1e-3 || next_surf_dist == 0){ // fix distance if too small 
  //     dir[0] = prev_dir[0];
  //     dir[1] = prev_dir[1];
  //     dir[2] = prev_dir[2];
  //   } 

  //   if (dir != prev_dir){
  //     history.reset(); // reset if direction changed
  //   }

  //   DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  //   sample_step_length = next_surf_dist/n_sample;

  // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;
  //   nrayfire +=1;
  // for (int i=0; i<=n_sample; i++){
  //   if (dir[2] < 0){
  //     ray_coords2 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
  //     for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //       sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //     }
  //     }
  //   else {
  //     ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] << std::endl;
  //     for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //       sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //     }
  //   }
  // }
  
  //   next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
  //   std::cout << std::endl;
  //   std::cout << "Distance to next surface = " << next_surf_dist << std::endl;
  //   std::cout << "Next Surface id - " << next_surf << std::endl;
  //   std::cout << "Volume ID - " << vol_h << std::endl;
  //   if (prev_surf == next_surf){
  //     DAG->next_vol(next_surf, vol_h, vol_h);
  //     DAG->point_in_volume(vol_h, origin, invol);
  //     std::cout << "Volume ID - " << vol_h << std::endl;
  //     std::cout << invol << " = in volume (1 for yes, 0 for no)" << std::endl;
  //   }

  // }
  // std::cout << "No. of ray_fire calls - " << nrayfire << std::endl; // print the No. of ray_fire calls
    
    
    return 0;
}

// Updates origin array for use in ray_fire 
void next_pt(double prev_pt[3], double origin[3], double next_surf_dist,
            double dir[3], std::ofstream &ray_intersect){
  // prev_pt -> Array of previous ray launch point
  // next_surf_dist -> Input distance to next surface (from ray_fire)
  // dir -> Array of direction vector of ray
  // write_stream -> filestream to write out ray_intersection points
  for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
    origin[i] = prev_pt[i] + (next_surf_dist * dir[i]);
    ray_intersect << origin[i] << ' ';
  }
  ray_intersect << std::endl;
  return;
}

double * vecNorm(double vector[3]){
  static double normalised_vector[3];
  double vector_mag;

  vector_mag = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
  for (int i=0; i<3; i++){
    normalised_vector[i] = vector[i]/vector_mag;
  }
  return normalised_vector;
}


// void next_dir(double raydir[3], const int qrymax, double dir_array[][3]){
//   for (int j=0; j<qrymax; j++)
//   {
//     for (int k=0; k<3; k++)
//   {
//     raydir[k] = dir_array[j+1][k] - dir_array[j][k]; 
//     //std::cout << raydir[k]<<std::endl;
//   }
//   //std::cout << std::endl;
//   }
//   return;
// }



// double dot_product(double vector_a[], double vector_b[]){
//    double product = 0;
//    for (int i = 0; i < 3; i++)
//    product = product + vector_a[i] * vector_b[i];
//    return product;
// }

// void reflect(double dir[3], double prev_dir[3], EntityHandle next_surf){
//   double normal[3];
//   prev_dir[0] = dir[0];
//   prev_dir[1] = dir[1];
//   prev_dir[2] = dir[2];
//   DAG->get_angle(next_surf, NULL, normal);
//   double reflect_dot = dot_product(dir, normal);
//   for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
//     dir[i] = dir[i] - 2*reflect_dot*normal[i];
//   }
//   return;
// }
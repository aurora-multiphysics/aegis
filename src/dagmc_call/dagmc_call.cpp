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
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/keywords/format.hpp>

//#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "settings.hpp"
#include "simpleLogger.h"

using namespace moab;

using moab::DagMC;
using moab::OrientedBoxTreeTool;
moab::DagMC* DAG;

void next_pt(double prev_pt[3], double origin[3], double next_surf_dist, 
                          double dir[3], std::ofstream &ray_intersect);
double dot_product(double vector_a[], double vector_b[]);
void reflect(double dir[3], double prev_dir[3], EntityHandle next_surf);
double * vecNorm(double vector[3]);


// logging to a file

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;

// LOG macros

//LOG_TRACE << "this is a trace message";
//LOG_DEBUG << "this is a debug message";
//LOG_WARNING << "this is a warning message";
//LOG_ERROR << "this is an error message";
//LOG_FATAL << "this is a fatal error message";

int main() {
  settings settings;

  settings.load_settings();
  LOG_WARNING << "h5m Faceted Geometry file to be used = " << settings.geo_input;
  //static const char input_file[] = settings.geo_input;
  //static const char ray_qry_exps[] = settings.ray_qry;
  static const char* input_file = settings.geo_input.c_str();
  static const char* ray_qry_exps = settings.ray_qry.c_str();
  

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


  LOG_WARNING << "No. of surfaces - " << n_surfs;
  LOG_WARNING << "No. of Volumes - " << n_vols;

  std::vector<std::vector<double>> rayqry; // xyz data from rayqry file




  // Read in qry data 
  std::ifstream ray_input(ray_qry_exps);
  double word;
  std::string line;

  if (ray_input) {
        while(getline(ray_input, line, '\n'))        
        {
            //create a temporary vector that will contain all the columns
            std::vector<double> tempVec;
            std::istringstream ss(line);
            //read word by word(or int by int) 
            while(ss >> word)
            {
                //LOG_WARNING<<"word:"<<word;
                //add the word to the temporary vector 
                tempVec.push_back(word);
            }             
            //now all the words from the current line has been added to the temporary vector 
            rayqry.emplace_back(tempVec);
        }    
    }
    else 
    {
        LOG_ERROR<<"file cannot be opened";
    }
    ray_input.close();

    //now you can do the whatever processing you want on the vector
    int j=0;
    int k;
    int qrymax = rayqry.size();
    double dir_array[rayqry.size()][3]; 
    LOG_WARNING << rayqry.size() << " - rayqry.size()" ;
    LOG_WARNING << qrymax << " - qrymax" ;

    //lets check out the elements of the 2D vector so the we can confirm if it contains all the right elements(rows and columns)
    for(std::vector<double> &newvec: rayqry)
    {
      k=0;
        for(const double &elem: newvec)
        {
          dir_array[j][k]=elem;
          //LOG_WARNING << elem;
          k+=1;
        }
        j+=1;
    }

    double raydirs[qrymax][3];
     for (int j=0; j<qrymax; j+=2)
     {
      for (int k=0; k<3; k++)
      {
        raydirs[j][k] = dir_array[j+1][k] - dir_array[j][k]; 
      }
      LOG_TRACE << raydirs[j][0] << ", " << raydirs[j][1] << ", " << raydirs[j][2]; // print out all ray directions from ray_qry

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
    ray_intersect << qryorigin[0] << ' ' << qryorigin[1] << ' ' << qryorigin[2] << ' ' ; // store first ray origin

    // check if in volume before ray_fire
    DAG->point_in_volume(vol_h, qryorigin, invol);
    LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
    LOG_WARNING << "Volume ID - " << vol_h ;
    double *normdir;
    //loop through qrydata
    for (int qrycount=0; qrycount < qrymax; qrycount+=2){
      normdir = vecNorm(raydirs[qrycount]);
      // print normalised ray directions

      //LOG_WARNING << normdir[0] ;
      //LOG_WARNING << normdir[1] ;
      //LOG_WARNING << normdir[2] ;
      //LOG_WARNING ;

      history.reset(); // reset history before launching ray
      // launch
      DAG->ray_fire(vol_h, qryorigin, normdir, next_surf, next_surf_dist, &history, 0, 1);

      // calculate next intersection point and write out to textfile
      for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
        intersect[i] = qryorigin[i] + (next_surf_dist * normdir[i]);
        ray_intersect << intersect[i] << ' ';
        }
        ray_intersect ;


    }
    







  // dir_mag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  // LOG_WARNING << "---- DIR_mag = " << dir_mag << " ----" ;

  // for (int i=0; i<3; ++i) { // calculate next ray launch point
  //   dir[i] = dir[i]/dir_mag;
  //   LOG_WARNING << "---- DIR[" << i << "] = " << dir[i] << " ----" ;
  // }

  // std::ofstream ray_coords1; // define stream ray_coords1
  // std::ofstream ray_coords2; // define stream ray_coords2
  // std::ofstream ray_intersect; // stream for ray-surface intersection points
  // ray_coords1.open("ray_coords1.txt"); // write out stream ray_coords1 to "ray_coords1.txt" file 
  // ray_coords2.open("ray_coords2.txt"); // write out stream ray_coords2 to "ray_coords1.txt" file 
  // ray_intersect.open("ray_intersect.txt"); // write out stream to "ray_intersect.txt" file

  // ray_intersect << origin[0] << ' ' << origin[1] << ' ' << origin[2] << ' ' ; // store first ray origin



  // // check if in volume before ray_fire
  // DAG->point_in_volume(vol_h, origin, invol);
  // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
  // LOG_WARNING << "Volume ID - " << vol_h ;
  // history.reset(); // reset history before launching any rays

  // // launch
  // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


  // sample_step_length = next_surf_dist/n_sample;
  // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));

  // for (int i=0; i<=n_sample; i++){
  //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
  //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //   }
  // }

  // nrayfire +=1;

  // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
  // LOG_WARNING ;
  // LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
  // LOG_WARNING << "Next Surface id - " << next_surf ;
  // DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
  // DAG->point_in_volume(vol_h, origin, invol);
  // LOG_WARNING << "Volume ID - " << vol_h ;
  // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;

  // //launch
  // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


  // sample_step_length = next_surf_dist/n_sample;
  // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;

  // for (int i=0; i<=n_sample; i++){
  //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
  //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //   }
  // }

  //   nrayfire +=1;

  // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);    
  // LOG_WARNING ;
  // DAG->point_in_volume(vol_h, origin, invol);
  // LOG_WARNING << "Volume ID - " << vol_h ;
  // LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
  // LOG_WARNING << "Next Surface id - " << next_surf ;
  // DAG->point_in_volume(vol_h, origin, invol);
  // LOG_WARNING << "Volume ID - " << vol_h ;
  // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
  
  
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
  //     ray_coords2 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
  //     for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //       sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //     }
  //     }
  //   else {
  //     ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
  //     for (int j=0; j<3; ++j) { // calculate next sample point along ray
  //       sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
  //     }
  //   }
  // }
  
  //   next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
  //   LOG_WARNING ;
  //   LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
  //   LOG_WARNING << "Next Surface id - " << next_surf ;
  //   LOG_WARNING << "Volume ID - " << vol_h ;
  //   if (prev_surf == next_surf){
  //     DAG->next_vol(next_surf, vol_h, vol_h);
  //     DAG->point_in_volume(vol_h, origin, invol);
  //     LOG_WARNING << "Volume ID - " << vol_h ;
  //     LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
  //   }

  // }
  // LOG_WARNING << "No. of ray_fire calls - " << nrayfire ; // print the No. of ray_fire calls
    
    
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
  ray_intersect ;
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
//     //LOG_WARNING << raydir[k];
//   }
//   //LOG_WARNING ;
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
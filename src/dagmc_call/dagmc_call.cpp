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
#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "settings.hpp"
#include "simpleLogger.h"
#include "equData.h"
#include "source.h"

using namespace moab;

using moab::DagMC;
using moab::OrientedBoxTreeTool;
moab::DagMC* DAG;

void next_pt(double prev_pt[3], double origin[3], double next_surf_dist, 
                          double dir[3], std::ofstream &ray_intersect);
double dot_product(double vector_a[], double vector_b[]);
void reflect(double dir[3], double prev_dir[3], EntityHandle next_surf);
double * vecNorm(double vector[3]);


// LOG macros

//LOG_TRACE << "this is a trace message";             ***WRITTEN OUT TO LOGFILE***
//LOG_DEBUG << "this is a debug message";             ***WRITTEN OUT TO LOGFILE***
//LOG_WARNING << "this is a warning message";         ***WRITTEN OUT TO CONSOLE AND LOGFILE***
//LOG_ERROR << "this is an error message";            ***WRITTEN OUT TO CONSOLE AND LOGFILE***
//LOG_FATAL << "this is a fatal error message";       ***WRITTEN OUT TO CONSOLE AND LOGFILE***

int main() {
  settings settings;
  settings.load_settings();
  LOG_WARNING << "h5m Faceted Geometry file to be used = " << settings.geo_input;
  LOG_WARNING << "ray directions file to be used = " << settings.ray_qry;
  static const char* input_file = settings.geo_input.c_str();
  static const char* ray_qry_exps = settings.ray_qry.c_str();
  std::string eqdsk_file = settings.eqdsk_file;

  DAG = new DagMC(); // New DAGMC instance
  DAG->load_file(input_file); // open test dag file
  DAG->init_OBBTree(); // initialise OBBTree 
  DagMC::RayHistory history; // initialise RayHistory object
  int vol_idx = 1; 
  int n_surfs = DAG->num_entities(2); // No. of surfaces in loaded geometry
  int n_vols = DAG->num_entities(3); // No. of volumes in loaded geometry
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double prev_dir[3];
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; // previous next_surf_dist before the next ray_fire call
  EntityHandle next_surf; // surface id 
  int invol; // point in volume result
  int nrayfire=0; // No. of ray_fire calls 
  EntityHandle prev_surf; // previous surface id
  double dir_mag; // magnitude of ray direction vector
  double normal[3]; // vector of surface facet normal intersected by last ray_fire
  double reflect_dot; // dot product of ray dir vector and normal vector 
  double *normdir;
  int lost_rays=0;


  LOG_WARNING << "No. of surfaces = " << n_surfs;
  LOG_WARNING << "No. of Volumes = " << n_vols;


  // --------------------------------------------------------------------------------------------------------
  if (settings.runcase=="rayqry") // rayqry test case
  {
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "--------------------RAYQRY CASE-----------------------" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  
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
          LOG_FATAL << "rayqry file cannot be opened";
      }
      ray_input.close();

      //now you can do the whatever processing you want on the vector
      int j=0;
      int k;
      int qrymax = rayqry.size();
      double dir_array[rayqry.size()][3]; 

      //lets check out the elements of the 2D vector so the we can confirm if it contains all the right elements(rows and columns)
      for(std::vector<double> &newvec: rayqry)
      {
        k=0;
          for(const double &elem: newvec)
          {
            dir_array[j][k]=elem;
            //LOG_TRACE << elem;
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
        //LOG_TRACE << raydirs[j][0] << ", " << raydirs[j][1] << ", " << raydirs[j][2]; // print out all ray directions from ray_qry

      }
      // define qrydata origin
      double qryorigin[3];
      double intersect[3];
      qryorigin[0] = dir_array[0][0];
      qryorigin[1] = dir_array[0][1];
      qryorigin[2] = dir_array[0][2];
      std::ofstream ray_intersect("ray_intersect.txt"); // stream for ray-surface intersection points
      ray_intersect << std::setprecision(5) << std::fixed;
      //ray_intersect << qryorigin[0] << ' ' << qryorigin[1] << ' ' << qryorigin[2] << ' ' ; // store first ray origin

      // check if in volume before ray_fire
      DAG->point_in_volume(vol_h, qryorigin, invol);
      LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
      LOG_WARNING << "Volume ID - " << vol_h ;
      //loop through qrydata
      for (int qrycount=0; qrycount < qrymax; qrycount+=2){
        normdir = vecNorm(raydirs[qrycount]);
        // print normalised ray directions

        LOG_TRACE << normdir[0] << ", " << normdir[1] << ", " << normdir[2]; // print out all ray directions from ray_qry

        history.reset(); // reset history before launching ray
        // launch
        DAG->ray_fire(vol_h, qryorigin, normdir, next_surf, next_surf_dist, &history, 0, 1);
        // Count number of lost rays
        if (next_surf == 0)
        {
          lost_rays += 1;
        }

        // calculate next intersection point and write out to textfile
        for (int i=0; i<3; ++i) 
        { // Calculate ray intersect point
          intersect[i] = qryorigin[i] + (next_surf_dist * normdir[i]);
          ray_intersect << intersect[i] << ' ';
        }
        ray_intersect << std::endl;
      }
      LOG_ERROR << "Lost ray count = " << lost_rays;
  }

  // --------------------------------------------------------------------------------------------------------
  else if (settings.runcase == "specific") // specific ray case
  {
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "--------------------SPECIFIC CASE---------------------" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    // Set the specific launch direction and origin
    double dir[3] = {0, 0, 1}; // ray launch direction
    double origin[3] = {0.0, 0.0, 0.0};  // ray launch origin

    int n_sample=10;  // number of samples along ray 
    double sample_step_length; // step length between sample positions along ray
    double sample_point[3]; // position of sample point along ray

    normdir = vecNorm(dir);
    dir[0] = normdir[0];
    dir[1] = normdir[1];
    dir[2] = normdir[2];        

    std::ofstream ray_coords1; // define stream ray_coords1
    std::ofstream ray_coords2; // define stream ray_coords2
    std::ofstream ray_intersect; // stream for ray-surface intersection points
    ray_coords1.open("ray_coords1.txt"); // write out stream ray_coords1 to "ray_coords1.txt" file 
    ray_coords2.open("ray_coords2.txt"); // write out stream ray_coords2 to "ray_coords1.txt" file 
    ray_intersect.open("ray_intersect.txt"); // write out stream to "ray_intersect.txt" file

    //ray_intersect << origin[0] << ' ' << origin[1] << ' ' << origin[2] << ' ' ; // store first ray origin

    EntityHandle surface;
    for (int i=1; i <=n_surfs; i++)
    { 
      surface = DAG->entity_by_index(2, i);
      std::cout << surface << std::endl;
    }

    DAG->write_mesh("dag.out", 1);
    
    //double pA[3] = {-50, 50, -100};
    //double pB[3] = {-50, 50, -121};
    double intersect_pt[3];
    
    double pSource[3] = {0, 0, 100};

    //boxSource spatialSource(pA, pB);
    pointSource spatialSource(pSource);
    int nSample = 10000;
    int ray_hit = 0;
    for (int i=0; i<nSample; i++)
    {
      spatialSource.get_isotropic_dir();

      DAG->ray_fire(vol_h, spatialSource.r, spatialSource.dir, next_surf, next_surf_dist, &history, 0, -1);
      history.reset();
      if (next_surf == 0)
      {
        next_surf_dist = 0;
      }
      else
      {
        ray_hit += 1; 
        for (int j=0; j<3; j++)
        {
          intersect_pt[j] = spatialSource.r[j] + next_surf_dist*spatialSource.dir[j];
        }
      ray_intersect << intersect_pt[0] << ' ' << intersect_pt[1] << ' ' << intersect_pt[2] << std::endl;
      }
      ray_coords1 << spatialSource.dir[0] << ' ' << spatialSource.dir[1] << ' ' << spatialSource.dir[2] << std::endl;
    
    }

    LOG_WARNING << "Number of rays hit = " << ray_hit;



    // // check if in volume before ray_fire
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // history.reset(); // reset history before launching any rays

    // // launch
    // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
    // LOG_WARNING << "Next Surface Distance = " << next_surf_dist;
    // sample_step_length = next_surf_dist/n_sample;
    // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));

    // for (int i=0; i<=n_sample; i++){
    //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
    //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
    //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    //   }
    //   ray_coords1 << std::endl;
    // }

    // nrayfire +=1;


    // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
    // LOG_WARNING << "Next Surface id - " << next_surf ;
    // DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)";

    // //launch
    // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
    // LOG_WARNING << "Next Surface Distance = " << next_surf_dist;
    // DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)


    // sample_step_length = next_surf_dist/n_sample;
    // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;

    // for (int i=0; i<=n_sample; i++){
    //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
    //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
    //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    //   }
    //   ray_coords1 << std::endl;
    // }

    // nrayfire +=1;

    // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);    
    // LOG_WARNING ;
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
    // LOG_WARNING << "Next Surface id - " << next_surf ;
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
    // LOG_WARNING << "No. of ray_fire calls - " << nrayfire ; // print the No. of ray_fire calls

    // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
    // LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
    // LOG_WARNING << "Next Surface id - " << next_surf ;
    // DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)";

    // //launch
    // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


    // sample_step_length = next_surf_dist/n_sample;
    // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;

    // for (int i=0; i<=n_sample; i++){
    //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
    //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
    //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    //   }
    //   ray_coords1 << std::endl;
    // }

    //   nrayfire +=1;

    //       next_pt(origin, origin, next_surf_dist, dir, ray_intersect);    
    // LOG_WARNING ;
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
    // LOG_WARNING << "Next Surface id - " << next_surf ;
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
    // LOG_WARNING << "No. of ray_fire calls - " << nrayfire ; // print the No. of ray_fire calls

    // next_pt(origin, origin, next_surf_dist, dir, ray_intersect);
    // LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
    // LOG_WARNING << "Next Surface id - " << next_surf ;
    // DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
    // DAG->point_in_volume(vol_h, origin, invol);
    // LOG_WARNING << "Volume ID - " << vol_h ;
    // LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)";

    // //launch
    // DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);


    // sample_step_length = next_surf_dist/n_sample;
    // std::copy(std::begin(origin), std::end(origin), std::begin(sample_point));origin;

    // for (int i=0; i<=n_sample; i++){
    //   ray_coords1 << sample_point[0] << ' ' << sample_point[1] << ' ' << sample_point[2] ;
    //   for (int j=0; j<3; ++j) { // calculate next sample point along ray
    //     sample_point[j] = sample_point[j] + (sample_step_length * dir[j]);
    //   }
    //   ray_coords1 << std::endl;
    // }

    //   nrayfire +=1;


  }

  // --------------------------------------------------------------------------------------------------------
  
  else if (settings.runcase=="eqdsk"){
    equData EquData;
    EquData.read_eqdsk(eqdsk_file);
    EquData.write_eqdsk_out();


    
    


  }
  
  else // No runcase specified
  {
    LOG_FATAL << "No runcase specified - please set runcase parameter as either 'specific' or 'rayqry'";
  }

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


// REFLECTION CODE 

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
    //   LOG_WARNING << "Distance to next surface = " << next_surf_dist ;
    //   LOG_WARNING << "Next Surface id - " << next_surf ;
    //   LOG_WARNING << "Volume ID - " << vol_h ;
    //   if (prev_surf == next_surf){
    //     DAG->next_vol(next_surf, vol_h, vol_h);
    //     DAG->point_in_volume(vol_h, origin, invol);
    //     LOG_WARNING << "Volume ID - " << vol_h ;
    //     LOG_WARNING << invol << " = in volume (1 for yes, 0 for no)" ;
    //   }

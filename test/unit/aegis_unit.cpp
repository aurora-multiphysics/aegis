#include <gtest/gtest.h>
#include <stdio.h>
#include "DagMC.hpp"
#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <fstream>
#include <vector>
#include <cmath>
#include "EquilData.h"
#include "Integrator.h"
#include "Source.h"
#include "CoordTransform.h"
#include "SimpleLogger.h"
#include "mpi.h"
#include "Particle.h"
 
using namespace moab;

using moab::DagMC;

moab::DagMC* DAG;
double * vecNorm(double vector[3]);

static const char input_file[] = "triple_block.h5m";
static const char scaled_fixed_end[] = "scaled_fixed_end.h5m";
static const char smardda_intersect[] ="smardda-intersect.txt";
static const char ray_qry_exps[] = "exps00000200.qry";
static const char sduct[] = "sduct.h5m";
static const char hashtag_mesh[] = "hashtag_mesh.h5m";

class aegisUnitTest: public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};



// ray_fire test with triple block geometry 
// simply tests ray_fire functionality
TEST_F(aegisUnitTest, Triple_Block_rayfire) {
  double max_tol = 1.0e-6;
  DAG = new DagMC();
  DAG->load_file(input_file); // open test dag file 
  DAG->init_OBBTree();
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};       // ray launch direction
  double origin[3] = {0.0, 0.0, 0.0};  // ray launch origin
  double next_surf_dist;
  EntityHandle next_surf;
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  double expected_next_surf_dist = 5.0;
  EXPECT_NEAR(expected_next_surf_dist, next_surf_dist, max_tol);
}


TEST_F(aegisUnitTest, Ray_Propagation_Pipe) {
  DAG = new DagMC();
  DAG->load_file(scaled_fixed_end); // open big pipe file 
  DAG->init_OBBTree();
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {-1.0, 1.0, 3.0};       // ray launch direction
  double origin[3] = {0.0, -15.0, -25.0};  // ray launch origin
  double prev_dir[3];
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; 
  int surfaces_intersected;
  int expected_surfaced_intersected;
  EntityHandle next_surf;
  int nrayfire=0; // No. of ray_fire calls 
  EntityHandle prev_surf; // previous surface id 
  double dir_mag; // magnitude of ray direction vector
  double normal[3]; // vector of surface facet normal intersected by last ray_fire
  double reflect_dot; // dot product of ray dir vector and normal vector  

  history.reset(); // reset history before launching any rays

  // launch ray

  dir_mag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  for (int i=0; i<3; ++i) { // calculate next ray launch point
    dir[i] = dir[i]/dir_mag;
  }
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  nrayfire +=1;
  for (int i=0; i<3; ++i) { // calculate next ray launch point
    origin[i] = origin[i] + (next_surf_dist * dir[i]);
  }
  DAG->next_vol(next_surf, vol_h, vol_h); // move to next volume id (determined from next_surf)
  // launch
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  nrayfire +=1;
  for (int i=0; i<3; ++i) { // calculate next ray launch point
    origin[i] = origin[i] + (next_surf_dist * dir[i]);
  }
  while (next_surf !=0){
    prev_dir[0] = dir[0];
    prev_dir[1] = dir[1];
    prev_dir[2] = dir[2];

    DAG->get_angle(next_surf, NULL, normal, &history);
    reflect_dot = dir[0]*normal[0] + dir[1]*normal[1] + dir[2]*normal[2];
    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      dir[i] = dir[i] - 2*reflect_dot*normal[i];
    }

    prev_surf = next_surf;
    prev_surf_dist = next_surf_dist;

    if (next_surf_dist < 1e-3 || next_surf_dist == 0){ // fix distance if too small 
      dir[0] = prev_dir[0];
      dir[1] = prev_dir[1];
      dir[2] = prev_dir[2];
    } 

    if (dir != prev_dir){
      history.reset(); // reset if direction changed
    }

    DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
    nrayfire +=1;

    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      origin[i] = origin[i] + (next_surf_dist * dir[i]);
    }
    if (prev_surf == next_surf){
      DAG->next_vol(next_surf, vol_h, vol_h);
    }
  }
EXPECT_EQ(nrayfire, 15);
}

TEST_F(aegisUnitTest, SMARDDA_comparison_test) {
  DAG = new DagMC();
  DAG->load_file(sduct); // open big pipe file 
  DAG->init_OBBTree();
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {-1.0, 1.0, 3.0};       // ray launch direction
  double origin[3] = {0.0, -15.0, -25.0};  // ray launch origin
  double prev_dir[3];
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; 
  EntityHandle next_surf;
  int nrayfire=0; // No. of ray_fire calls 
  EntityHandle prev_surf; // previous surface id 
  double dir_mag; // magnitude of ray direction vector
  double normal[3]; // vector of surface facet normal intersected by last ray_fire

  std::vector<std::vector<double>> rayqry; // xyz data from rayqry file
  std::ifstream ray_input(ray_qry_exps);
  std::string line;
  double word;

  if (ray_input) 
  {
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

    std::ifstream smardda_pts(smardda_intersect);
    std::vector<std::vector<double>> smardda_qry;

    if (smardda_pts) {
      while(getline(smardda_pts, line, '\n'))        
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
        smardda_qry.emplace_back(tempVec);
      }    
    }
    else 
    {
        std::cout<<"file cannot be opened"<<std::endl;
    }
    smardda_pts.close();

    //now you can do the whatever processing you want on the vector
    int j=0;
    int k;
    int rayqrymax = rayqry.size();
    double dir_array[rayqry.size()][3]; 

    //check out the elements of the 2D vector so the we can confirm if it contains all the right elements(rows and columns)
    for(std::vector<double> &newvec: rayqry)
    {
      k=0;
        for(const double &elem: newvec)
        {
          dir_array[j][k]=elem;
          //std::cout << elem<< std::endl;
          k+=1;
        }
        //std::cout<<std::endl;
        j+=1;
    }
    double raydirs[rayqrymax][3];
    for (int j=0; j<rayqrymax; j+=2)
     {
      for (int k=0; k<3; k++)
      {
        raydirs[j][k] = dir_array[j+1][k] - dir_array[j][k]; 
      }
    }

    // define qrydata origin
    double qryorigin[3];
    double intersect[3];
    qryorigin[0] = dir_array[0][0];
    qryorigin[1] = dir_array[0][1];
    qryorigin[2] = dir_array[0][2];
    double *normdir;
    double xdiff, ydiff, zdiff;
    double max_xdiff=0, max_ydiff=0, max_zdiff=0;
    int smardda_count=0; 

    //loop through qrydata
    for (int qrycount=0; qrycount < rayqrymax; qrycount+=2){
      normdir = vecNorm(raydirs[qrycount]);
      history.reset(); // reset history before launching ray

      // launch
      DAG->ray_fire(vol_h, qryorigin, normdir, next_surf, next_surf_dist, &history, 0, 1);

      // calculate next intersection point and write out to textfile
      for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
        intersect[i] = qryorigin[i] + (next_surf_dist * normdir[i]);
        }
      // compare element by element difference between smardda_pts and rayqry
      xdiff = fabs(intersect[0] - smardda_qry[smardda_count][0]);
      ydiff = fabs(intersect[1] - smardda_qry[smardda_count][1]);
      zdiff = fabs(intersect[2] - smardda_qry[smardda_count][2]);
      if (xdiff > max_xdiff) {
        max_xdiff = xdiff;
      }
      if (ydiff > max_ydiff) {
        max_ydiff = ydiff;
      }
      if (zdiff > max_zdiff) {
        max_zdiff = zdiff;
      }

      smardda_count += 1;
    }
    
    double max_tol = 0.6;    

    std::cout << "max_xdiff = " << max_xdiff << std::endl; 
    std::cout << "max_ydiff = " << max_ydiff << std::endl; 
    std::cout << "max_zdiff = " << max_zdiff << std::endl; 

    EXPECT_TRUE(max_xdiff < max_tol);
    EXPECT_TRUE(max_ydiff < max_tol);
    EXPECT_TRUE(max_zdiff < max_tol);

  
}


TEST_F(aegisUnitTest, eqdsk_read) {
  MPI_Init(NULL, NULL);
  auto equilibrium = std::make_shared<EquilData>();
  equilibrium->read_eqdsk("test.eqdsk");
  eqdskData eqdsk = equilibrium->get_eqdsk_struct();

  std::string header = " disr     610.00 msec    nw      nh      vde       0  65 129";
  // Test if header is correctly read
  EXPECT_TRUE(eqdsk.header == header);

  // Test if integer values are correctly read in

  EXPECT_EQ(equilibrium->nw, 65);
  EXPECT_EQ(equilibrium->nh, 129);
  EXPECT_EQ(eqdsk.nbdry, 89);
  EXPECT_EQ(eqdsk.nlim, 57);

  // Test size of arrays 
  EXPECT_EQ(eqdsk.fpol.size(), 65);
  EXPECT_EQ(eqdsk.pres.size(), 65);
  EXPECT_EQ(eqdsk.ffprime.size(), 65);
  EXPECT_EQ(eqdsk.pprime.size(), 65);
  EXPECT_EQ(eqdsk.qpsi.size(), 65);
  EXPECT_EQ(eqdsk.rbdry.size(), 89);
  EXPECT_EQ(eqdsk.zbdry.size(), 89);
  EXPECT_EQ(eqdsk.rlim.size(), 57);
  EXPECT_EQ(eqdsk.zlim.size(), 57);


  EXPECT_EQ(eqdsk.psi.size()*eqdsk.psi[0].size(), 8385);

  // Test random elements of array
  EXPECT_FLOAT_EQ(eqdsk.fpol[35], 33.2359842);
  EXPECT_FLOAT_EQ(eqdsk.pres[13], 457077.309);
  EXPECT_FLOAT_EQ(eqdsk.ffprime[54], 252.582128);
  EXPECT_FLOAT_EQ(eqdsk.pprime[11], 1742.17328);  
  EXPECT_FLOAT_EQ(eqdsk.psi[56][1], 1.39236796);
  EXPECT_FLOAT_EQ(eqdsk.psi[1][24], -6.80347259);
  EXPECT_FLOAT_EQ(eqdsk.psi[12][102], -2.22815321);
  EXPECT_FLOAT_EQ(eqdsk.qpsi[52], 1.64743552);
  EXPECT_FLOAT_EQ(eqdsk.rbdry[4], 7.9430907);  
  EXPECT_FLOAT_EQ(eqdsk.zbdry[39], 0.986207476);
  EXPECT_FLOAT_EQ(eqdsk.rlim[27], 557.200115);
  EXPECT_FLOAT_EQ(eqdsk.zlim[1], -149.995359);
  MPI_Finalize();
}


TEST_F(aegisUnitTest, cart_to_polar_transform){
  std::vector<double> input = {2.4, 5.3, -2.0};
  std::vector<double> output;
  output.reserve(3);

  output = CoordTransform::cart_to_polar(input);

  EXPECT_FLOAT_EQ(output[0], 5.8180752);
  EXPECT_FLOAT_EQ(output[1], -2);
  EXPECT_FLOAT_EQ(output[2], -1.1455913);

}

TEST_F(aegisUnitTest, polar_to_cart_transform){
  std::vector<double> input = {5.8180752, -2.0, -1.1455913};
  std::vector<double> output;
  output.reserve(3);

  output = CoordTransform::polar_to_cart(input);

  EXPECT_FLOAT_EQ(output[0], 2.4);
  EXPECT_FLOAT_EQ(output[1], 5.3);
  EXPECT_FLOAT_EQ(output[2], -2.0);

}

TEST_F(aegisUnitTest, polar_to_flux_transform){
  MPI_Init(NULL, NULL);
  auto equilibrium = std::make_shared<EquilData>();
  equilibrium->read_eqdsk("test.eqdsk");
  equilibrium->init_interp_splines();
  equilibrium->centre(1);

  std::vector<double> input = {2.4, 5.3, -2.0};
  std::vector<double> output;
  std::string direction;
  output.reserve(3);

  output = CoordTransform::polar_to_flux(input, equilibrium);

  EXPECT_FLOAT_EQ(output[0], -2.5001357);
  EXPECT_FLOAT_EQ(output[1], 2.1932724);
  EXPECT_FLOAT_EQ(output[2], -2);
  MPI_Finalize();
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



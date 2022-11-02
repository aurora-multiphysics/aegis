#include <gtest/gtest.h>
#include "DagMC.hpp"
#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"

using namespace moab;

using moab::DagMC;

moab::DagMC* DAG;

static const char input_file[] = "triple_block.h5m";
double eps = 1.0e-6;

class DagmcSimpleTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};


// Static DAGMC version call
TEST(DagmcCallTest, VersionCall) {
  moab::DagMC dagmc_instance;
  ErrorCode rval = dagmc_instance.load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  } 


TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc) {
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  // make new dagmc into that moab
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);

  ErrorCode rval;

  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
}


// ray_fire test with triple block geometry 
// simply tests ray_fire functionality
TEST_F(DagmcSimpleTest, Triple_Block_rayfire) {
  DAG = new DagMC();
  DAG->load_file(input_file); // open test dag file 
  DAG->init_OBBTree();
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};       // ray along x direction
  double origin[3] = {-5.0, 0.0, 0.0};  // origin at -5 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  double expected_next_surf_dist = 10.0;
  EXPECT_NEAR(expected_next_surf_dist, next_surf_dist, eps);
}


TEST_F(DagmcSimpleTest, Ray_Propagation_Pipe) {
  DAG = new DagMC();
  DAG->load_file("big_pipe.h5m"); // open test dag file 
  DAG->init_OBBTree();
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 1.0};       // ray launch direction
  double origin[3] = {-6.0, 0.0, -6.0};  // ray launch origin
  double next_surf_dist;
  int surfaces_intersected;
  int expected_surfaced_intersected;
  EntityHandle next_surf;
  int nrayfire=0; // No. of ray_fire calls 
  EntityHandle prev_surf; // previous surface id
  double prev_surf_dist; // previous next_surf_dist before the next ray_fire call

  history.reset(); // reset history before launching any rays

  // launch ray
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
    dir[0]=-dir[0]; // flip Z-component to keep ray in pipe
    prev_surf = next_surf;
    prev_surf_dist = next_surf_dist;
    DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
    nrayfire +=1;
    if (next_surf_dist < 1e-3){ // fix distance if too small 
      next_surf_dist = prev_surf_dist;
    }  
    for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
      origin[i] = origin[i] + (next_surf_dist * dir[i]);
    }
    if (prev_surf == next_surf){
      DAG->next_vol(next_surf, vol_h, vol_h);
    }
  EXPECT_EQ(nrayfire, 8);
  }
}
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


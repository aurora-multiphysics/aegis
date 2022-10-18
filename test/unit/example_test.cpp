#include "DataClass.hpp"
#include <gtest/gtest.h>
#include "DagMC.hpp"
#include <iostream>

TEST(DagmcCallTest, VersionCall) {
  moab::DagMC dagmc_instance;
  float version = dagmc_instance.version();
  EXPECT_FLOAT_EQ(version, 3.2) << "DAGMC version verified as 3.2";
  } 


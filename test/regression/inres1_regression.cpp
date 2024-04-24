#include <gtest/gtest.h>
#include <stdio.h>
#include <bits/stdc++.h>
#include <assert.h>
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
#include <algorithm>
#include <unordered_map>
#include <time.h>
#include <any>
#include <mpi.h>


#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "Inputs.h"
#include "SimpleLogger.h"
#include "EquilData.h"
#include "Source.h"
#include "Integrator.h"
#include "CoordTransform.h"
#include "alglib/interpolation.h"
#include "Particle.h"
#include "ParticleSimulation.h"
 
 using namespace moab;

 using moab::DagMC;


double dot_product(std::vector<double> vector_a, std::vector<double> vector_b);

 class aegisRegressionTests: public ::testing::Test {
  protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
 };



 TEST_F(aegisRegressionTests, inres1)
 {
  MPI_Init(NULL, NULL);

  std::string smarddaFile = "smardda_inres1_no_shadowQ.txt";
  double qvalue;
  std::vector<std::pair<double,double>> smardda_qValues;
  std::string line;
  std::ifstream smardda_qValues_in(smarddaFile);
  if (smardda_qValues_in)
  {
    while (std::getline(smardda_qValues_in, line, '\n'))
    {
      std::istringstream ss(line);
      std::vector<double> tempVec;
      while(ss >> qvalue)
      {
        tempVec.push_back(qvalue);  
      }
      smardda_qValues.push_back(std::make_pair(tempVec[1], tempVec[0]));
    }
  }
  else 
  {
    FAIL() << "Cannot find '" << smarddaFile << "'";
  }
//-------------- RUN AEGIS -------------
  std::string configFilename = "aegis_settings.json";
  auto configFile = std::make_shared<JsonHandler>(configFilename);

  auto equilibrium = std::make_shared<EquilData>();
  equilibrium->setup(configFile);
  equilibrium->move();
  equilibrium->psiref_override();
  equilibrium->init_interp_splines();
  equilibrium->centre(1);
  equilibrium->write_bfield();

  ParticleSimulation aegis(configFile, equilibrium);

  if (std::filesystem::exists(configFilename)){
    aegis.Execute_serial();
  }

  else {
    FAIL() << "Cannot find '" << configFilename << "'";
  }
//-------------- COMPARE AEGIS AGAINST SMARDDA -------------

  std::sort(aegis.psiQ_values.begin(), aegis.psiQ_values.end());
  std::sort(smardda_qValues.begin(), smardda_qValues.end());

  std::cout << std::endl << "Printing first 10 elements..." << std::endl;
  for (int i=0; i<10; i++)
  {
    std::cout << "Element - " << i << std::endl;
    std::cout << "PSI = " << aegis.psiQ_values[i].first << " Q = " << aegis.psiQ_values[i].second << " --> AEGIS" << std::endl;
    std::cout << "PSI = " << smardda_qValues[i].first << " Q = " << smardda_qValues[i].second << " --> SMARDDA" << std::endl;
    std::cout << std::endl;
  }  

  double Q_rel_sum = 0.0;

  std::cout << aegis.target_num_facets();

  for (int i=0; i<aegis.target_num_facets(); i++){
    Q_rel_sum += std::pow((aegis.psiQ_values[i].second - smardda_qValues[i].second), 2);
  }


  double L2_NORM = std::sqrt( (1.0/aegis.target_num_facets())*Q_rel_sum );
  const double EXPECTED_L2_NORM = 349791;

  const auto AEGIS_MAX = *std::max_element(aegis.psiQ_values.begin(),aegis.psiQ_values.end(),[](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
  const auto SMARDDA_MAX = *std::max_element(smardda_qValues.begin(),smardda_qValues.end(),[](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });

  double percentTol = 7; // 5% tolerance
  double MAX_REL_ERROR = fabs(std::fabs((AEGIS_MAX.second - SMARDDA_MAX.second)/SMARDDA_MAX.second)*100);
  double L2_NORM_ERROR = fabs(std::fabs((L2_NORM - EXPECTED_L2_NORM)/EXPECTED_L2_NORM)*100);
  std::cout << "---------------------------" << std::endl;
  std::cout << "MAX Q AEGIS = " << AEGIS_MAX.second << std::endl;
  std::cout << "MAX Q SMARDDA = " << SMARDDA_MAX.second << std::endl;
  std::cout << "MAX Q %ERROR = " << MAX_REL_ERROR << std::endl;
  std::cout << "L2_NORM = " << L2_NORM << std::endl;
  std::cout << "L2_NORM %Error = " << L2_NORM_ERROR << std::endl;
  std::cout << "---------------------------" << std::endl << std::endl;
  
  EXPECT_LT(MAX_REL_ERROR, percentTol);
  ASSERT_LT(L2_NORM_ERROR, percentTol);
  MPI_Finalize();

 }


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
#include <limits>  // for min/max values
#include <set>
#include <vector>
//#include <gtest/gtest.h>

#include "DagMC.hpp"

int main() { 
    moab::DagMC dagmc_instance;
    std::cout << dagmc_instance.version() << std::endl;
//    TEST(DagmcCallTest, VersionCall) {
//        EXPECT_STRCASEEQ(version_string,"3.2");
//	} 
//     int Factorial(int n);
//     TEST(FactorialTest, HandlesZeroInput) {
//       EXPECT_EQ(Factorial(0), 1);
// }
    return 0;
}

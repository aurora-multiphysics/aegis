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

#include "DagMC.hpp"

int main() { 

    moab::DagMC dagmc_instance;
    std:: cout << "DAGMC Version installed -> " <<  moab::DagMC::version() << std::endl ;
    return 0;
}

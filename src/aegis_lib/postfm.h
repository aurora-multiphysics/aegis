#ifndef postfm__
#define postfm__

#include <vector>
#include <variant>
#include <cmath>
#include <string>
#include "equData.h"


std::vector<double> cartesian_to_polar_toroidal(std::vector<double> inputVector,
                                                   std::string direction);

// Coordinate system transforms

std::vector<double> polar_to_flux(std::vector<double> inputVector, std::string direction,
                                  equData& EquData);



#endif
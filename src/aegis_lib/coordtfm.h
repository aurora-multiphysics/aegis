#ifndef coordtfm__
#define coordtfm__

#include <vector>
#include <variant>
#include <cmath>
#include <string>
#include "equData.h"

namespace coordTfm
{

std::vector<double> cart_to_polar(std::vector<double> inputVector,
                                                   std::string direction);

// Coordinate system transforms

std::vector<double> polar_to_flux(std::vector<double> inputVector, std::string direction,
                                  equData& EquData);

}

#endif
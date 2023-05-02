#ifndef coordtfm__
#define coordtfm__

#include <vector>
#include <variant>
#include <cmath>
#include <string>
#include "equData.h"

namespace coordTfm
{

// Cartesian to polar toroidal (direction = "backwards" for polar->cartesian)
std::vector<double> cart_to_polar(std::vector<double> inputVector,
                                  std::string direction);

// Polar toroidal to flux coordinates defined by psi spline
std::vector<double> polar_to_flux(std::vector<double> inputVector, std::string direction,
                                  equData& EquData);

}

#endif
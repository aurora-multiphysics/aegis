#ifndef CoordTransform__
#define CoordTransform__

#include <vector>
#include <variant>
#include <cmath>
#include <string>
#include "EquilData.h"

namespace CoordTransform
{


// Cartesian to polar toroidal (direction = "backwards" for polar->cartesian)
std::vector<double> cart_to_polar(std::vector<double> inputVector,
                                  std::string direction);
std::vector<double> cart_to_polar(double e0, double e1, double e2, std::string direction);

// Polar toroidal to flux coordinates defined by psi spline
std::vector<double> polar_to_flux(std::vector<double> inputVector, std::string direction,
                                  const std::shared_ptr<EquilData>& equilibrium);
std::vector<double> polar_to_flux(double e0, double e1, double e2, std::string direction, const std::shared_ptr<EquilData>& equilibrium);
};

#endif
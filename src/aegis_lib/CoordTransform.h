#ifndef CoordTransform__
#define CoordTransform__

#include <vector>
#include <variant>
#include <cmath>
#include <string>
#include "EquilData.h"
#include "alglib/interpolation.h"


namespace CoordTransform
{


// Cartesian to polar toroidal (direction = "backwards" for polar->cartesian)
std::vector<double> cart_to_polar(std::vector<double> inputVector);
std::vector<double> cart_to_polar_xy(std::vector<double> inputVector);
std::vector<double> cart_to_polar(double e0, double e1, double e2);

std::vector<double> polar_to_cart(std::vector<double> inputVector);
std::vector<double> polar_to_cart(double e0, double e1, double e2);

// Polar toroidal to flux coordinates defined by psi spline
std::vector<double> polar_to_flux(std::vector<double> inputVector,
                                  const std::shared_ptr<EquilData>& equilibrium);


std::vector<double> polar_to_flux(double e0, double e1, double e2, const std::shared_ptr<EquilData>& equilibrium);

};

#endif
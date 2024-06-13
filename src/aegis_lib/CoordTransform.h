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
std::vector<double> cart_to_polar(const std::vector<double> &xyz);

// polar (R,Z) coords from (X,Y,Z)
inline std::vector<double> cart_to_polar_xy(const std::vector<double> &xyz)
{
  double r = sqrt(pow(xyz[0], 2) + pow(xyz[1], 2));
  std::vector<double> rz = {r, xyz[2]};
  return rz;
};

// get polar angle (PHI) from (X,Y)
inline double cart_to_polar_phi(const std::vector<double> &xyz)
{
  double phi = atan2(-xyz[1], xyz[0]);
  return phi;
};


std::vector<double> cart_to_polar(const double &e0, const double &e1, const double &e2);

std::vector<double> polar_to_cart(const std::vector<double> &inputVector);
std::vector<double> polar_to_cart(const double &e0, const double &e1, const double &e2);

// Polar toroidal to flux coordinates defined by psi spline
std::vector<double> polar_to_flux(const std::vector<double> &inputVector,
                                  const std::shared_ptr<EquilData>& equilibrium);


std::vector<double> polar_to_flux(const double &e0, const double &e1, const double &e2, const std::shared_ptr<EquilData>& equilibrium);

};

#endif
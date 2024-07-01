#include <vector>
#include <variant>
#include <cmath>
#include <string>

#include "CoordTransform.h"
#include "EquilData.h"
#include "alglib/interpolation.h"
#include "SimpleLogger.h"

std::vector<double>
CoordTransform::cart_to_polar(const std::vector<double> & inputVector)
{
  double x = inputVector[0];
  double y = inputVector[1];
  double z = inputVector[2];
  double r = sqrt(pow(x, 2) + pow(y, 2)); // calculate
  double phi = atan2(-y, x);              // calculate phi
  std::vector<double> outputVector = {r, z, phi};
  return outputVector;
}

std::vector<double>
CoordTransform::polar_to_cart(const std::vector<double> & inputVector)
{
  double r = inputVector[0];
  double z = inputVector[1];
  double phi = inputVector[2];
  double x = r * cos(phi);  // calculate x
  double y = -r * sin(phi); // calculate y
  std::vector<double> outputVector = {x, y, z};
  return outputVector;
}

std::vector<double>
CoordTransform::cart_to_polar(const double & e0, const double & e1, const double & e2)
{
  double x = e0;
  double y = e1;
  double z = e2;
  double r = sqrt(pow(x, 2) + pow(y, 2)); // calculate
  double phi = atan2(-y, x);              // calculate phi
  std::vector<double> outputVector = {r, z, phi};
  return outputVector;
}

std::vector<double>
CoordTransform::polar_to_cart(const double & e0, const double & e1, const double & e2)
{
  double r = e0;
  double z = e1;
  double phi = e2;
  double x = r * cos(phi);  // calculate x
  double y = -r * sin(phi); // calculate y
  std::vector<double> outputVector = {x, y, z};
  return outputVector;
}

std::vector<double>
CoordTransform::polar_to_flux(const std::vector<double> & inputVector,
                              const std::shared_ptr<EquilData> & equilibrium)
{
  double r = inputVector[0];
  double z = inputVector[1];
  double phi = inputVector[2];
  double psi =
      alglib::spline2dcalc(equilibrium->psiSpline, r, z); // spline interpolation of psi(R,Z)
  double theta = atan2(z - equilibrium->zcen, r - equilibrium->rcen);
  if (theta < -M_PI_2)
  {
    theta = 2 * M_PI + theta;
  }
  std::vector<double> outputVector = {-psi, theta, phi};
  return outputVector;
}

std::vector<double>
CoordTransform::polar_to_flux(const double & e0, const double & e1, const double & e2,
                              const std::shared_ptr<EquilData> & equilibrium)
{
  double r = e0;
  double z = e1;
  double phi = e2;
  double psi =
      alglib::spline2dcalc(equilibrium->psiSpline, r, z); // spline interpolation of psi(R,Z)
  double theta = atan2(z - equilibrium->zcen, r - equilibrium->rcen);
  if (theta < -M_PI_2)
  {
    theta = 2 * M_PI + theta;
  }
  std::vector<double> outputVector = {-psi, theta, phi};
  return outputVector;
}

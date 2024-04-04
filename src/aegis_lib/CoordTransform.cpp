#include <vector>
#include <variant>
#include <cmath>
#include <string>

#include "CoordTransform.h"
#include "EquilData.h"
#include "alglib/interpolation.h"
#include "SimpleLogger.h"


std::vector<double> CoordTransform::cart_to_polar(std::vector<double> inputVector,
                                                   std::string direction)
{
  std::vector<double> outputVector(3);
  double r; // polar r
  double phi; // polar phi
  double x; // cart x
  double y; // cart y
  double z; // cart z

  if (direction == "backwards")
  {
    r = inputVector[0];
    z = inputVector[1];
    phi = inputVector[2];

    x = r*cos(phi); // calculate x

    y = -r*sin(phi); // calculate y

    outputVector[0] = x;
    outputVector[1] = y;
    outputVector[2] = z;
  }
  else
  {
    x = inputVector[0];
    y = inputVector[1];
    z = inputVector[2];

    r = sqrt(pow(x,2) + pow(y,2)); // calculate
    phi = atan2(-y,x); // calculate phi

    outputVector[0] = r;
    outputVector[1] = z;
    outputVector[2] = phi;
  }
  return outputVector;
}

std::vector<double> CoordTransform::cart_to_polar(double e0, double e1, double e2, std::string direction)
{
  std::vector<double> outputVector(3);
  double r; // polar r
  double phi; // polar phi
  double x; // cart x
  double y; // cart y
  double z; // cart z

  if (direction == "backwards")
  {
    r = e0;
    z = e1;
    phi = e2;

    x = r*cos(phi); // calculate x

    y = -r*sin(phi); // calculate y

    outputVector[0] = x;
    outputVector[1] = y;
    outputVector[2] = z;
  }
  else
  {
    x = e0;
    y = e1;
    z = e2;

    r = sqrt(pow(x,2) + pow(y,2)); // calculate
    phi = atan2(-y,x); // calculate phi

    outputVector[0] = r;
    outputVector[1] = z;
    outputVector[2] = phi;
  }
  return outputVector;
}

std::vector<double> CoordTransform::polar_to_flux(std::vector<double> inputVector,
                                            std::string direction, const std::shared_ptr<EquilData>& equilibrium)
{
  std::vector<double> outputVector(3);
  double r; // local polar r
  double z; // local polar z
  double phi; // local polar phi
  double psi; // local psi
  double theta; // local theta

  if (direction == "backwards") // backwards transform (flux -> polar coords) TODO
  {
    psi = inputVector[0];
    theta = inputVector[1];
    phi = inputVector[2];


  }
  else // fowards transform (polar -> flux coords)
  {
    r = inputVector[0];
    z = inputVector[1];
    phi = inputVector[2];

    psi = alglib::spline2dcalc(equilibrium->psiSpline, r, z); // spline interpolation of psi(R,Z)

    theta = atan2(z-equilibrium->zcen, r-equilibrium->rcen);
    if (theta < -M_PI_2)
    {
      theta = 2*M_PI+theta;
    }
  }


  outputVector[0] = -psi;
  outputVector[1] = theta;
  outputVector[2] = phi;

  return outputVector;
}

std::vector<double> CoordTransform::polar_to_flux(double e0, double e1, double e2,
                                            std::string direction, const std::shared_ptr<EquilData>& equilibrium)
{
  std::vector<double> outputVector(3);
  double r; // local polar r
  double z; // local polar z
  double phi; // local polar phi
  double psi; // local psi
  double theta; // local theta

  if (direction == "backwards") // backwards transform (flux -> polar coords) TODO
  {
    psi = e0;
    theta = e1;
    phi = e2;


  }
  else // fowards transform (polar -> flux coords)
  {
    r = e0;
    z = e1;
    phi = e2;

    psi = alglib::spline2dcalc(equilibrium->psiSpline, r, z); // spline interpolation of psi(R,Z)

    theta = atan2(z-equilibrium->zcen, r-equilibrium->rcen);
    if (theta < -M_PI_2)
    {
      theta = 2*M_PI+theta;
    }
  }


  outputVector[0] = -psi;
  outputVector[1] = theta;
  outputVector[2] = phi;

  return outputVector;
}

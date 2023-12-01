#include <vector>
#include <variant>
#include <cmath>
#include <string>

#include "coordtfm.h"
#include "equData.h"
#include "alglib/interpolation.h"
#include "simpleLogger.h"


std::vector<double> coordTfm::cart_to_polar(std::vector<double> inputVector,
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

std::vector<double> coordTfm::polar_to_flux(std::vector<double> inputVector,
                                            std::string direction, equData& EquData)
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

    psi = alglib::spline2dcalc(EquData.psiSpline, r, z); // spline interpolation of psi(R,Z)

    theta = atan2(z-EquData.zcen, r-EquData.rcen);
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

// coordSystem =
// 1 - Cartesian
// 2 - Polar-Toroidal
// 3 - Flux Coordinates
vecTfm::vecTfm(std::vector<double> vector, int coordSystem, equData &EquData)
{
  switch (coordSystem)
  {
    case 1:
      cartCoord = vector;
      polarCoord = coordTfm::cart_to_polar(cartCoord, "forwards");
      fluxCoord = coordTfm::polar_to_flux(polarCoord, "forwards", EquData);
    case 2:
      polarCoord = vector;
      fluxCoord = coordTfm::polar_to_flux(polarCoord, "forwards", EquData);
      cartCoord = coordTfm::cart_to_polar(cartCoord, "backwards");
    // case 3:
    //   flux = vector;
    //   // TODO backwards tfm from flux
    default:
      LOG_INFO << "No coordinate system provided for vecTfm constructor. Assume Cartesian";
      cartCoord = vector;
  }
}


// return cartesian vector
std::vector<double> vecTfm::cart()
{
  return cartCoord;
}

std::vector<double> vecTfm::polar()
{
  // not sure why this is broken atm
  // cart and flux returns work but this is broken for some reason
  std::cout << polarCoord[0] << polarCoord[1] << polarCoord[2] << std::endl;
  return polarCoord;
}

std::vector<double> vecTfm::flux()
{
  return fluxCoord;
}


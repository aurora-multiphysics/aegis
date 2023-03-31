#include <vector>
#include <variant>
#include <cmath>
#include <string>

#include "coordtfm.h"
#include "equData.h" 
#include "alglib/interpolation.h"


std::vector<double> coordTfm::cartesian_to_polar_toroidal(std::vector<double> inputVector,
                                                   std::string direction)
{
  std::vector<double> outputVector;
  outputVector.reserve(3);
  double r; // polar r 
  double zeta; // polar zeta
  double x; // cart x 
  double y; // cart y
  double z; // cart z

  if (direction == "backwards")
  {
    r = inputVector[0]; 
    z = inputVector[1]; 
    zeta = inputVector[2];

    x = r*cos(zeta); // calculate x
    y = -r*sin(zeta); // calculate y

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
    zeta = atan2(-y,x); // calculate zeta 

    outputVector[0] = r;
    outputVector[1] = z; 
    outputVector[2] = zeta;
  } 
  return outputVector;
}

std::vector<double> coordTfm::polar_to_flux(std::vector<double> inputVector, std::string direction,
                                  equData& EquData)
{
  std::vector<double> outputVector;
  outputVector.reserve(3);
  double r; // local polar r
  double z; // local polar z
  double zeta; // local polar zeta
  double psi; // local psi  
  double theta; // local theta 

  if (direction == "backwards") // backwards transform (flux -> polar coords) TODO
  {
    psi = inputVector[0];
    theta = inputVector[1];
    zeta = inputVector[2];

    
  }
  else // fowards transform (polar -> flux coords)
  {
    r = inputVector[0];
    z = inputVector[1];
    zeta = inputVector[2];

    psi = alglib::spline2dcalc(EquData.psiSpline, r, z); // spline interpolation of psi(R,Z)
    theta = atan2(z-EquData.zcen, r-EquData.rcen);
    if (theta < -M_PI_2) 
    {
      theta = 2*M_PI+theta;
    }
  }

  outputVector[0] = psi;
  outputVector[1] = theta;
  outputVector[2] = zeta;

  return outputVector;
}

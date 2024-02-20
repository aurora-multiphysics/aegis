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

// Polar toroidal to flux coordinates defined by psi spline
std::vector<double> polar_to_flux(std::vector<double> inputVector, std::string direction,
                                  EquilData& EquData);

};


  class vecTfm
  {
    private:
    std::vector<double> cartCoord;
    std::vector<double> polarCoord;
    std::vector<double> fluxCoord;

    public: 
    vecTfm(std::vector<double> vector, int coordSystem, EquilData &EquData);
    std::vector<double> cart();
    std::vector<double> polar();
    std::vector<double> flux();
  };

#endif
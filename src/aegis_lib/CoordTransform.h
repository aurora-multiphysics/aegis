#ifndef CoordTransform__
#define CoordTransform__

#include <vector>
#include <variant>
#include <cmath>
#include <string>
#include "EquilData.h"
#include "alglib/interpolation.h"

class Position3D : public AegisBase
{
  public:
  Position3D() = default;
  template <typename T>
  Position3D(const T coords, coordinateSystem coordSys) : _e0(coords[0]), _e1(coords[1]), _e2(coords[2]), _coordSys(coordSys) {};
  Position3D(double e0, double e1, double e2, coordinateSystem coordSys) : _e0(e0), _e1(e1), _e2(e2), _coordSys(coordSys) {};


  // return current coordinate system 
  coordinateSystem coordSys() { return _coordSys; } 

  // return STL array of coords
  std::array<double, 3> coords() { return {_e0, _e1, _e2}; } 

  // return STL vector of coords
  std::vector<double> coords_vec() { return {_e0, _e1, _e2}; }

  // update coordinates without changing system 
  template <typename T> void update_coords(T newCoords)
  {
    _e0 = newCoords[0]; 
    _e1 = newCoords[1];
    _e2 = newCoords[2]; 
  }

  // update coordinates without changing system
  void update_coords(double e0New, double e1New, double e2New);
  
  // update coordinates with a new coordinate system. Provide coordinates of new system.
  template <typename T> void update_coords(T newCoords, coordinateSystem newCoordSys)
  {
    _e0 = newCoords[0]; 
    _e1 = newCoords[1];
    _e2 = newCoords;
    _coordSys = newCoordSys; 
  }

  // update coordinates with a new coordinate system. Provide coordinates of new system.
  void update_coords(double e0New, double e1New, double e2New, coordinateSystem newCoordSys);

  // generalised coordinates
  double _e0 = 0.0; // first general coordinate
  double _e1 = 0.0 ; // second general coordinate 
  double _e2 = 0.0; // third general coordinate


  // operator overloading
  // template <typename T> 
  // Position3D& operator=(T other)
  // {
  //   _e0 = other[0];
  //   _e1 = other[1];
  //   _e2 = other[2];
  //   return *this;
  // }

  Position3D& operator+=(Position3D);
  Position3D& operator+=(double v);
  Position3D& operator-=(Position3D);
  Position3D& operator-=(double v);
  Position3D& operator*=(double v);
  Position3D& operator/=(double v);

  const double& operator[](int i) const
  {
    switch (i) {
    case 0:
      return _e0;
    case 1:
      return _e1;
    case 2:
      return _e2;
    default:
      throw std::out_of_range {"Index in Position must be between 0 and 2."};
    }
  }
  double& operator[](int i)
  {
    switch (i) {
    case 0:
      return _e0;
    case 1:
      return _e1;
    case 2:
      return _e2;
    default:
      throw std::out_of_range {"Index in Position must be between 0 and 2."};
    }
  }


  protected:

  private:
  coordinateSystem _coordSys = coordinateSystem::CARTESIAN;
};


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

inline Position3D cart_to_polar_xy(const Position3D &pos)
{
  double r = sqrt(pow(pos[0], 2) + pow(pos[1], 2));
  double z = pos[2];
  Position3D newPos(r, z, 0.0, coordinateSystem::POLAR);
  return newPos;
}

// get polar angle (PHI) from (X,Y)
inline double cart_to_polar_phi(const std::vector<double> &xyz)
{
  double phi = atan2(-xyz[1], xyz[0]);
  return phi;
};

Position3D cart_to_polar(const Position3D & pos);
Position3D polar_to_cart(const Position3D & pos);
Position3D polar_to_flux(Position3D & pos, const std::shared_ptr<EquilData> & equilibrium);

std::vector<double> cart_to_polar(const double &e0, const double &e1, const double &e2);

std::vector<double> polar_to_cart(const std::vector<double> &inputVector);
std::vector<double> polar_to_cart(const double &e0, const double &e1, const double &e2);

// Polar toroidal to flux coordinates defined by psi spline
std::vector<double> polar_to_flux(const std::vector<double> &inputVector,
                                  const std::shared_ptr<EquilData>& equilibrium);

std::vector<double> polar_to_flux(const double &e0, const double &e1, const double &e2, const std::shared_ptr<EquilData>& equilibrium);
};

#endif
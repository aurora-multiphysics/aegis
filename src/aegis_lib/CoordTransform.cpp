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

// Position 3D methods

Position3D &
Position3D::operator+=(Position3D other)
{
  _e0 += other._e0;
  _e1 += other._e1;
  _e2 += other._e2;
  return *this;
}

Position3D &
Position3D::operator+=(double v)
{
  _e0 += v;
  _e1 += v;
  _e2 += v;
  return *this;
}

Position3D &
Position3D::operator-=(Position3D other)
{
  _e0 -= other._e0;
  _e1 -= other._e1;
  _e2 -= other._e2;
  return *this;
}

Position3D &
Position3D::operator-=(double v)
{
  _e0 -= v;
  _e1 -= v;
  _e2 -= v;
  return *this;
}

Position3D &
Position3D::operator*=(double v)
{
  _e0 *= v;
  _e1 *= v;
  _e2 *= v;
  return *this;
}

Position3D &
Position3D::operator/=(double v)
{
  _e0 /= v;
  _e1 /= v;
  _e2 /= v;
  return *this;
}

// update coordinates without changing system
void
Position3D::update_coords(double e0New, double e1New, double e2New)
{
  _e0 = e0New;
  _e1 = e1New;
  _e2 = e2New;
}

// update coordinates with a new coordinate system. Provide coordinates of new system.
void
Position3D::update_coords(double e0New, double e1New, double e2New, coordinateSystem newCoordSys)
{
  _e0 = e0New;
  _e1 = e1New;
  _e2 = e2New;
  _coordSys = newCoordSys;
}

// Position 3D coordTransforms

Position3D
CoordTransform::cart_to_polar(const Position3D & pos)
{
  double x = pos[0];
  double y = pos[1];
  double z = pos[2];
  double r = sqrt(pow(x, 2) + pow(y, 2)); // calculate
  double phi = atan2(-y, x);              // calculate phi
  Position3D newPos(r, z, phi, coordinateSystem::POLAR);
  return newPos;
}

Position3D
CoordTransform::polar_to_cart(const Position3D & pos)
{
  double r = pos[0];
  double z = pos[1];
  double phi = pos[2];
  double x = r * cos(phi);  // calculate x
  double y = -r * sin(phi); // calculate y
  Position3D newPos(x, y, z, coordinateSystem::CARTESIAN);
  return newPos;
}

Position3D
CoordTransform::polar_to_flux(Position3D & pos, const std::shared_ptr<EquilData> & equilibrium)
{
  // transform to polar before flux
  if (pos.coordSys() == coordinateSystem::CARTESIAN)
  {
    pos = CoordTransform::cart_to_polar(pos);
  }

  double r = pos[0];
  double z = pos[1];
  double phi = pos[2];
  double psi =
      alglib::spline2dcalc(equilibrium->psiSpline, r, z); // spline interpolation of psi(R,Z)
  double theta = atan2(z - equilibrium->zcen, r - equilibrium->rcen);
  if (theta < -M_PI_2)
  {
    theta = 2 * M_PI + theta;
  }
  Position3D newPos(-psi, theta, phi, coordinateSystem::FLUX);
  return newPos;
}

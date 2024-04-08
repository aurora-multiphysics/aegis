#include "Source.h"
#include "CoordTransform.h"
#include <cctype>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>

// define plane via cartesian plane equation ax + by + cz + d = 0
TriangleSource::TriangleSource(std::vector<double> xyz1, std::vector<double> xyz2,
                               std::vector<double> xyz3, moab::EntityHandle handle,
                               std::string launchType)
{
  xyzA = xyz1;
  xyzB = xyz2;
  xyzC = xyz3;

  std::vector<double> line1(3), line2(3);
  for (int i = 0; i < 3; i++)
  {
    line1[i] = xyzB[i] - xyzA[i];
    line2[i] = xyzC[i] - xyzA[i];
  }

  std::vector<double> normalVec(3);

  // recover vector normal to plane (Cross product of two edges in the triangle)
  normalVec[0] = line1[1] * line2[2] - line1[2] * line2[1];
  normalVec[1] = line1[2] * line2[0] - line1[0] * line2[2];
  normalVec[2] = line1[0] * line2[1] - line1[1] * line2[0];

  // recover constant D in plane equation
  D = -(normalVec[0] * xyzA[0] + normalVec[1] * xyzA[1] + normalVec[2] * xyzA[2]);
  normal = normalVec;
  double magnitude = sqrt(pow(normalVec[0], 2) + pow(normalVec[1], 2) + pow(normalVec[2], 2));
  normalVec[0] = normalVec[0] / magnitude;
  normalVec[1] = normalVec[1] / magnitude;
  normalVec[2] = normalVec[2] / magnitude;
  unitNormal = normalVec;

  entityHandle = handle;

  if (launchType == "random")
  {
    launchPos = random_pt();
  }
  else
  {
    launchPos = centroid();
  }
}

std::vector<double>
TriangleSource::random_pt()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0, 1);
  double randomA = dist(rng);
  double randomB = dist(rng);

  if ((randomA + randomB) > 1)
  {
    randomA = 1 - randomA;
    randomB = 1 - randomB;
  }

  std::vector<double> randomPt(3);

  for (int i = 0; i < 3; i++)
  {
    randomPt[i] = xyzA[i] + randomA * (xyzB[i] - xyzA[i]) + randomB * (xyzC[i] - xyzA[i]);
  }
  return randomPt;
}

std::vector<double>
TriangleSource::centroid()
{
  std::vector<double> centroid(3);
  double x, y, z; // xyz coords of centroid
  x = (xyzA[0] + xyzB[0] + xyzC[0]) / 3;
  y = (xyzA[1] + xyzB[1] + xyzC[1]) / 3;
  z = (xyzA[2] + xyzB[2] + xyzC[2]) / 3;

  centroid[0] = x;
  centroid[1] = y;
  centroid[2] = z;

  return centroid;
}

// set heatflux and B.n on surface element at particle launch position
void
Sources::set_heatflux_params(const std::shared_ptr<EquilData> & equilibrium,
                             const std::string formula)
{
  double psid;
  std::vector<double> polarPos, fluxPos;
  std::vector<double> bField;
  std::vector<double> bFieldXYZ;

  polarPos = CoordTransform::cart_to_polar(launchPos, "forwards");
  fluxPos = CoordTransform::polar_to_flux(polarPos, "forwards", equilibrium);

  bField = equilibrium->b_field(polarPos, "forwards");
  bField = equilibrium->b_field_cart(bField, polarPos[2], 0);
  double product = 0;
  for (int i = 0; i < 3; i++)
  {
    product = product + bField[i] * unitNormal[i];
  }

  Bn = product; // store B.n

  psi = fluxPos[0]; // store psi at particle start
  psid = psi + equilibrium->psibdry;
  Q = equilibrium->omp_power_dep(psid, Bn, "exp"); // store Q at particle start
}

void
Sources::update_heatflux(double newHeatflux)
{
  Q = newHeatflux;
}

std::vector<double>
Sources::launch_pos()
{
  return launchPos;
}

double
Sources::BdotN()
{
  return Bn;
}

double
Sources::heatflux()
{
  return Q;
}

moab::EntityHandle
Sources::entity_handle()
{
  return entityHandle;
}

double
Sources::get_psi()
{
  return psi;
}

std::vector<double>
Sources::get_normal()
{
  return normal;
}
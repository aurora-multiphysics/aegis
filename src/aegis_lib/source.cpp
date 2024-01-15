#include <stdio.h>
#include <iostream>
#include <cctype>
#include <random>
#include <math.h>
#include "source.h"


pointSource::pointSource(std::vector<double> xyz)
{
  for (int i=0; i<xyz.size(); i++)
  {
    r[i] = xyz[i];
  }
}

void pointSource::set_dir(double newDir[3])
{
  dir[0] = newDir[0];
  dir[1] = newDir[1];
  dir[2] = newDir[2];
}


// sample a random direction isotropically from point 
void pointSource::get_isotropic_dir()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0,1.0);
  double theta = 2 * M_PI * dist(rng);
  double phi = acos(1 - 2 * dist(rng));
  dir[0] = sin(phi) * cos(theta);
  dir[1] = sin(phi) * sin(theta);
  dir[2] = cos(phi);
}

void pointSource::get_hemisphere_surface_dir(double surfaceNormal[3])
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0,1.0);
  double theta = 2 * M_PI * dist(rng);
  double phi = acos(1 - 2 * dist(rng));
  dir[0] = sin(phi) * cos(theta);
  dir[1] = sin(phi) * sin(theta);
  dir[2] = cos(phi);
}


// boxSource constructor, a square plane can be described by providing position vectors
// with one dimension that is constant 
boxSource::boxSource(double xyz1[3], double xyz2[3])
{
  for (int i=0; i<3; i++)
  {  
    pA[i] = xyz1[i];
    pB[i] = xyz2[i];
  }  
}



void boxSource::get_pt()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0,1);
  double random[3]={dist(rng), dist(rng), dist(rng)};

  for (int i=0; i<3; i++){
    pR[i] = pA[i] + random[i] * (pB[i]-pA[i]);
  }
}

void boxSource::get_dir()
{
  std::vector<double> test(3);
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0.0,1.0);
  double theta = 2 * M_PI * dist(rng);
  double phi = acos(1 - 2 * dist(rng));
  dir[0] = sin(phi) * cos(theta);
  dir[1] = sin(phi) * sin(theta);
  dir[2] = cos(phi);
}

/////////

// define plane via cartesian plane equation ax + by + cz + d = 0
triSource::triSource(std::vector<double> xyz1, std::vector<double> xyz2, std::vector<double> xyz3, moab::EntityHandle handle)
{
  xyzA = xyz1;
  xyzB = xyz2;
  xyzC = xyz3;

  std::vector<double> line1(3), line2(3);
  for (int i=0; i<3; i++)
  {
    line1[i] = xyzB[i] - xyzA[i];
    line2[i] = xyzC[i] - xyzA[i];
  }

  std::vector<double> normalVec(3);


  // recover vector normal to plane (Cross product of two edges in the triangle)
  normalVec[0] = line1[1]*line2[2] - line1[2]*line2[1];
  normalVec[1] = line1[2]*line2[0] - line1[0]*line2[2];
  normalVec[2] = line1[0]*line2[1] - line1[1]*line2[0];

  // recover constant D in plane equation
  D = -(normalVec[0]*xyzA[0] + normalVec[1]*xyzA[1] + normalVec[2]*xyzA[2]);
  this->normal = normalVec;
  double magnitude = sqrt(pow(normalVec[0],2) + pow(normalVec[1],2) + pow(normalVec[2],2));
  normalVec[0] = normalVec[0]/magnitude;
  normalVec[1] = normalVec[1]/magnitude;
  normalVec[2] = normalVec[2]/magnitude; 
  this->unitNormal = normalVec;

  this->entityHandle = handle;

}



std::vector<double> triSource::random_pt()
{
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_real_distribution<double> dist(0,1);
  double randomA = dist(rng);
  double randomB = dist(rng);

  if ((randomA + randomB) > 1) 
  {
    randomA = 1 - randomA;
    randomB = 1 - randomB;
  }

  std::vector<double> randomPt(3);

  for(int i=0; i<3; i++)
  {
    randomPt[i] = xyzA[i] + randomA*(xyzB[i] - xyzA[i]) + randomB*(xyzC[i] - xyzA[i]);
  }
  return randomPt;
}

std::vector<double> triSource::centroid()
{
  std::vector<double> centroid(3);
  double x, y, z; // xyz coords of centroid
  x = (xyzA[0] + xyzB[0] + xyzC[0])/3; 
  y = (xyzA[1] + xyzB[1] + xyzC[1])/3; 
  z = (xyzA[2] + xyzB[2] + xyzC[2])/3; 

  centroid[0] = x;
  centroid[1] = y;
  centroid[2] = z;

  return centroid;
}



void triSource::dagmcInstance(moab::DagMC* DAG)
{
  DAGInstance = DAG;
  DAGInstance->write_mesh("dag.out", 1);
}

double triSource::dot_product(std::vector<double> externalVector) 
{
  double product = 0;
  for (int i; i<3; i++)
  {
    product = product + externalVector[i]*this->unitNormal[i];
  }
  return product;
}
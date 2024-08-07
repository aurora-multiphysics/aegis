#include <gtest/gtest.h>
#include <stdio.h>
#include "DagMC.hpp"
#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <fstream>
#include <vector>
#include <cmath>
#include "EquilData.h"
#include "Integrator.h"
#include "Source.h"
#include "CoordTransform.h"
#include "SimpleLogger.h"
#include "mpi.h"
 
using namespace moab;

using moab::DagMC;

class convexHullUnitTest: public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};

struct Point
{
  double r,z;
};

int orientation(Point p, Point q, Point r);
void jarvis_march(std::vector<std::vector<double>> &polarCoords);



TEST_F(convexHullUnitTest, recoverSurfaces) {
  moab::ErrorCode rval;

  // initialise DAGMC
  auto DAG = std::make_unique<DagMC>();
  rval = DAG->load_file("watertight_MASTU.h5m"); EXPECT_EQ(0, rval);
  rval = DAG->init_OBBTree(); EXPECT_EQ(0, rval);
  
  // get list of volumes and surfaces
  moab::Range surfaces, volumes;
  rval = DAG->setup_geometry(surfaces, volumes); EXPECT_EQ(0, rval);

  // get implicit complement volume
  auto GTT = DAG->geom_tool();
  auto MOAB = DAG->moab_instance();
  EntityHandle implicitComplement;
  GTT->get_implicit_complement(implicitComplement);
  MOAB->list_entity(implicitComplement);


  // rval = DAG->moab_instance()->list_entity(implicitComplement); EXPECT_EQ(0, rval);

  // get children meshsets of implicit complement
  std::vector<EntityHandle> children;
  rval = DAG->moab_instance()->get_child_meshsets(implicitComplement, children, 1); EXPECT_EQ(0, rval);
  
  // print information about childsets
  // DAG->moab_instance()->list_entities(children.data(), children.size());

  // get vertices in implicit complement
  std::vector<EntityHandle> vertices;
  int sense = 0;
  for(const auto & child:children)
  {
    rval = DAG->moab_instance()->get_entities_by_type(child, MBVERTEX, vertices, false); EXPECT_EQ(0, rval);
  }
  std::vector<double> vertexCoords(vertices.size()*3);
  rval = DAG->moab_instance()->get_coords(&vertices[0], vertices.size(), vertexCoords.data()); EXPECT_EQ(0, rval);

  std::exit(1);

  // write out implicit complement coords
  double f = 0.95;
  double x = 0;
  double y = 0;
  double z = 0.485034;

  std::ofstream implicitComplementCoordsCart("implicit_complement_coords_xyz.txt");
  std::ofstream implicitComplementCoordsPolar("implicit_complement_coords_rz.txt");
  std::vector<std::vector<double>> coords(vertices.size());
  std::vector<double> polarCoords;
  std::vector<std::vector<double>> points;
  int index = 0;
  for (auto &i:coords)
  {
    i.push_back(vertexCoords[index]);
    i.push_back(vertexCoords[index+1]);
    i.push_back(vertexCoords[index+2]);
    index +=3;

    // i[0] = f*(i[0]-x)+x;
    // i[1] = f*(i[1]-y)+y;
    // i[2] = f*(i[2]-z)+z;

    implicitComplementCoordsCart << i[0] << " " << i[1] << " " << i[2] << std::endl;

    polarCoords = CoordTransform::cart_to_polar_xy(i);
    implicitComplementCoordsPolar << polarCoords[0] << " " << polarCoords[1] << std::endl;
    points.push_back(polarCoords);
  }

  jarvis_march(points);

}


int orientation(Point p, Point q, Point r)
{
  double val = (q.z - p.z) * (r.r - q.r) -
               (q.r - p.r) * (r.z - q.z);
  if (val == 0) return 0; // collinear
  return (val > 0)? 1:2; // clock or anticlock wise
}

void jarvis_march(std::vector<std::vector<double>> &polarCoords)
{
  int n = polarCoords.size();

  std::vector<Point> points;
  for (const auto &i:polarCoords)
  {
    Point temp; 
    temp.r = i[0];
    temp.z = i[1];
    points.push_back(temp);
  }


  std::vector<Point> convexHull; // initialise results
  
  // find the leftmost point
  int l = 0;
  for (int i = 1; i < points.size(); i++)
  {
    if (points[i].r < points[l].r) { l = i; }
  }

  // start from leftmost point moving counter clockwise until we reach start again
  int p = l, q;

  do 
  {
    // Add current point to our convex hull result
    convexHull.push_back(points[p]);

    // search for a point q such that the orientation(p, q, x) is counterclockwise for all points x.
    q = (p+1)%n;
    for (int i = 0; i < n; i++)
    {
      // if i more counterclockwise than current q then update q
      if (orientation(points[p], points[i], points[q]) == 2) { q = i; }
    }

    // now q is the most counterclockwise with respect to p.
    // set p as q fr next iteration so that q is added to the convex hull
    p = q;
    
  } while (p != l); // while we dont come back to the first point

  // print our result to a file
  std::ofstream hull("convex_hull_pts.txt");
  for (const auto &i:convexHull)
  {
    hull << i.r << " " << i.z << std::endl;
  }

}

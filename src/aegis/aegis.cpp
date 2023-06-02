#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>  // for min/max values
#include <set>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>


#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "settings.hpp"
#include "simpleLogger.h"
#include "equData.h"
#include "source.h"
#include "integrator.h"
#include "coordtfm.h"
#include "alglib/interpolation.h"



using namespace moab;
using namespace coordTfm;

using moab::DagMC;
using moab::OrientedBoxTreeTool;
moab::DagMC* DAG;

void next_pt(double prev_pt[3], double origin[3], double next_surf_dist,
                          double dir[3], std::ofstream &ray_intersect);
double dot_product(double vector_a[], double vector_b[]);
double dot_product(std::vector<double> vector_a, std::vector<double> vector_b);
void reflect(double dir[3], double prev_dir[3], EntityHandle next_surf);
double * vecNorm(double vector[3]);


// LOG macros

//LOG_TRACE << "this is a trace message";             ***WRITTEN OUT TO LOGFILE***
//LOG_DEBUG << "this is a debug message";             ***WRITTEN OUT TO LOGFILE***
//LOG_WARNING << "this is a warning message";         ***WRITTEN OUT TO CONSOLE AND LOGFILE***
//LOG_ERROR << "this is an error message";            ***WRITTEN OUT TO CONSOLE AND LOGFILE***
//LOG_FATAL << "this is a fatal error message";       ***WRITTEN OUT TO CONSOLE AND LOGFILE***

int main() {
  settings settings;
  settings.load_settings();
  LOG_WARNING << "h5m Faceted Geometry file to be used = " << settings.geo_input;


  static const char* input_file = settings.geo_input.c_str();
  static const char* ray_qry_exps = settings.ray_qry.c_str();
  std::string eqdsk_file = settings.eqdsk_file;
  moab::Range Surfs, Vols, Facets, Facet_vertices;
  DAG = new DagMC(); // New DAGMC instance
  DAG->load_file(input_file); // open test dag file
  DAG->init_OBBTree(); // initialise OBBTree
  DAG->setup_geometry(Surfs, Vols);
  //DAG->create_graveyard();
  DAG->remove_graveyard();
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, Facets);
  LOG_WARNING << "No of triangles in geometry " << Facets.size();

  //DAG->moab_instance()->get_entities_by_dimension()

  moab::EntityHandle triangle_set, vertex_set;
  DAG->moab_instance()->create_meshset( moab::MESHSET_SET, triangle_set );
  DAG->moab_instance()->add_entities(triangle_set, Facets);
  int num_facets;
  DAG->moab_instance()->get_number_entities_by_handle(triangle_set, num_facets);

  DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, Facet_vertices);
  LOG_WARNING << "Number of vertices in geometry " << Facet_vertices.size();
  std::cout << std::endl;
  DAG->moab_instance()->create_meshset( moab::MESHSET_SET, vertex_set );
  DAG->moab_instance()->add_entities(vertex_set, Facet_vertices);
  DAG->moab_instance()->get_number_entities_by_dimension(vertex_set, 1 ,num_facets);
  moab::Range vertAdjs;

  moab::Range ents;
  DAG->moab_instance()->get_entities_by_handle(0, ents);
  for (auto i:Facets)
  {
    moab::Range verts;
    DAG->moab_instance()->get_adjacencies(&i, 1, 0, false, verts);
    //std::cout << "Tri " << DAG->moab_instance()->id_from_handle(i) << " vertex adjacencies:" << std::endl;
    std::vector<double> coords(3*verts.size());
    DAG->moab_instance()->get_coords(verts, &coords[0]);
    for (int j=0; j<verts.size(); j++)
    {
      //std::cout << "Vertex " << DAG->moab_instance()->id_from_handle(verts[j]) << "(" <<  
      //coords[0+j] << ", " << coords[1+j] << ", "  << coords[2+j] << ")" << std::endl;
    }
      //std::cout << std::endl;

    // for (auto j:verts)
    // {
    //   double coords[3];
    //   DAG->moab_instance()->get_coords(&j, 3);
    //   std::cout << "Coords of each vertex" << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
    // }
  }



  //DAG->write_mesh("dag.out", 1);
  EntityHandle prev_surf; // previous surface id
  EntityHandle next_surf; // surface id
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; // previous next_surf_dist before the next ray_fire call
  DagMC::RayHistory history; // initialise RayHistory object
  EntityHandle vol_h = DAG->entity_by_index(3, 1);



  double dir_mag; // magnitude of ray direction vector
  double reflect_dot; // dot product of ray dir vector and normal vector
  double *normdir;
  int lost_rays=0;

  // --------------------------------------------------------------------------------------------------------
  if (settings.runcase=="rayqry") // rayqry test case
  {
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "--------------------RAYQRY CASE-----------------------" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

    std::vector<std::vector<double>> rayqry; // xyz data from rayqry file
    // Read in qry data
    std::ifstream ray_input(ray_qry_exps);
    double word;
    std::string line;
    if (ray_input) {
          while(getline(ray_input, line, '\n'))
          {
              //create a temporary vector that will contain all the columns
              std::vector<double> tempVec;
              std::istringstream ss(line);
              //read word by word(or int by int)
              while(ss >> word)
              {
                  //LOG_WARNING<<"word:"<<word;
                  //add the word to the temporary vector
                  tempVec.push_back(word);
              }
              //now all the words from the current line has been added to the temporary vector
              rayqry.emplace_back(tempVec);
          }
      }
      else
      {
          LOG_FATAL << "rayqry file cannot be opened";
      }
      ray_input.close();

      //now you can do the whatever processing you want on the vector
      int j=0;DAG->write_mesh("dag.out", 1);
      int k;
      int qrymax = rayqry.size();
      double dir_array[rayqry.size()][3];

      //lets check out the elements of the 2D vector so the we can confirm if it contains all the right elements(rows and columns)
      for(std::vector<double> &newvec: rayqry)
      {
        k=0;
          for(const double &elem: newvec)
          {
            dir_array[j][k]=elem;
            //LOG_TRACE << elem;
            k+=1;
          }
          j+=1;
      }

      double raydirs[qrymax][3];
      for (int j=0; j<qrymax; j+=2)
      {
        for (int k=0; k<3; k++)
        {
          raydirs[j][k] = dir_array[j+1][k] - dir_array[j][k];
        }
        //LOG_TRACE << raydirs[j][0] << ", " << raydirs[j][1] << ", " << raydirs[j][2]; // print out all ray directions from ray_qry

      }
      // define qrydata origin
      double qryorigin[3];
      double intersect[3];
      qryorigin[0] = dir_array[0][0];
      qryorigin[1] = dir_array[0][1];
      qryorigin[2] = dir_array[0][2];
      std::ofstream ray_intersect("ray_intersect.txt"); // stream for ray-surface intersection points
      ray_intersect << std::setprecision(5) << std::fixed;
      //ray_intersect << qryorigin[0] << ' ' << qryorigin[1] << ' ' << qryorigin[2] << ' ' ; // store first ray origin

      //loop through qrydata
      for (int qrycount=0; qrycount < qrymax; qrycount+=2){
        normdir = vecNorm(raydirs[qrycount]);
        // print normalised ray directions

        LOG_TRACE << normdir[0] << ", " << normdir[1] << ", " << normdir[2]; // print out all ray directions from ray_qry

        history.reset(); // reset history before launching ray
        // launch
        DAG->ray_fire(vol_h, qryorigin, normdir, next_surf, next_surf_dist, &history, 0, 1);
        // Count number of lost rays
        if (next_surf == 0)
        {
          lost_rays += 1;
        }

        // calculate next intersection point and write out to textfile
        for (int i=0; i<3; ++i)
        { // Calculate ray intersect point
          intersect[i] = qryorigin[i] + (next_surf_dist * normdir[i]);
          ray_intersect << intersect[i] << ' ';
        }
        ray_intersect << std::endl;
      }
      LOG_ERROR << "Lost ray count = " << lost_rays;
  }

  // --------------------------------------------------------------------------------------------------------
  else if (settings.runcase == "specific") // specific ray case
  {
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "--------------------SOURCE DEFINED CASE---------------------" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    std::ofstream ray_coords1; // define stream ray_coords1
    std::ofstream ray_coords2; // define stream ray_coords2
    std::ofstream ray_intersect; // stream for ray-surface intersection points
    int reflection = 1;
    int nSample = 10000; // Number of rays sampled
    ray_coords1.open("ray_coords1.txt"); // write out stream ray_coords1 to "ray_coords1.txt" file
    ray_coords2.open("ray_coords2.txt"); // write out stream ray_coords2 to "ray_coords1.txt" file
    ray_intersect.open("ray_intersect.txt"); // write out stream to "ray_intersect.txt" file


    surfaceIntegrator integrator(Facets);
    //std::map<EntityHandle, int>::iterator it = integrator.nRays.begin();
    int lostRays=0;

    double intersect_pt[3];
    std::vector<double> pSource(3); 
    pSource = {0, 0, -30};
    EntityHandle meshElement;
    double distanceMesh;
    pointSource spatialSource(pSource);
    double reflected_dir[3];
    double reflect_dot;
    double surface_normal[3];
    EntityHandle obb_root;
    DAG->next_vol(Surfs[0], vol_h, vol_h);
    DAG->get_root(vol_h, obb_root);
    //DAG->obb_tree()->print(obb_root, std::cout, false);
    EntityHandle facet_hit;
    for (int i=0; i<nSample; i++)
    {
      spatialSource.get_isotropic_dir();
      history.reset();
      DAG->ray_fire(vol_h, spatialSource.r, spatialSource.dir, next_surf, next_surf_dist, &history, 0, 1);
      history.get_last_intersection(facet_hit);

      // std::cout << "obb_tree() " << DAG->obb_tree() << std::endl;
      // std::cout << "get_geom_tag() " << DAG->geom_tag() << std::endl;
      if (next_surf == 0)
      {
        next_surf_dist = 0;
        lostRays +=1;
      }
      else
      {
        for (int j=0; j<3; j++)
        {
          intersect_pt[j] = spatialSource.r[j] + next_surf_dist*spatialSource.dir[j];
        }
        integrator.count_hit(facet_hit);


        ray_intersect << spatialSource.r[0] << ' ' << spatialSource.r[1] << ' ' << spatialSource.r[2] << std::endl;
        ray_intersect << intersect_pt[0] << ' ' << intersect_pt[1] << ' ' << intersect_pt[2] << std::endl;

      ray_coords1 << spatialSource.dir[0] << ' ' << spatialSource.dir[1] << ' ' << spatialSource.dir[2] << std::endl;



      // REFLECTION CODE
        // if (reflection==1)
        // {
          // DAG->get_angle(next_surf, NULL, surface_normal, &history);
          // reflect_dot = dot_product(spatialSource.dir, surface_normal);

          // reflected_dir[0] = spatialSource.dir[0] - 2*reflect_dot*surface_normal[0];
          // reflected_dir[1] = spatialSource.dir[1] - 2*reflect_dot*surface_normal[1];
          // reflected_dir[2] = spatialSource.dir[2] - 2*reflect_dot*surface_normal[2];

          // DAG->ray_fire(vol_h, intersect_pt, reflected_dir, next_surf, next_surf_dist, &history, 0, 1);
          // history.get_last_intersection(facet_hit);
          //   if (reflected_dir != spatialSource.dir)
          //   {
          //     history.reset();
          //   }
          //   for (int j=0; j<3; j++)
          //   {
          //     intersect_pt[j] = intersect_pt[j] + next_surf_dist*reflected_dir[j];
          //   }
          //   DAG->closest_to_location(vol_h, intersect_pt, distanceMesh, &meshElement);
          //   integrator.count_hit(facet_hit);

          //   ray_intersect << intersect_pt[0] << ' ' << intersect_pt[1] << ' ' << intersect_pt[2] << std::endl;
          // //}
      }
    }


    LOG_WARNING << "Number of rays lost = " << lostRays;
    LOG_WARNING << "Number of rays launched = " << nSample;
    LOG_WARNING << "Number of ray-facet intersections = " << integrator.raysHit;

    integrator.facet_values(integrator.nRays);






  }

  // --------------------------------------------------------------------------------------------------------

  else if (settings.runcase=="eqdsk")
  {
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "--------------------READING EQDSK---------------------" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    equData EquData;


    EquData.read_eqdsk(eqdsk_file);
    //EquData.write_eqdsk_out();

    EquData.init_interp_splines();
    EquData.gnuplot_out();
    EquData.centre(1);

    // attempting to create a trace through magnetic field 
    // Currently not working

    
    bool plotRZ = true;
    bool plotXYZ = true;
    EquData.write_bfield(plotRZ, plotXYZ);
    std::vector<double> cartPosSource(3);
    double s, ds;
    double normB[3];
    double norm;
    std::ofstream trace1("trace1.txt");
    std::ofstream trace2("trace2.txt");
    std::ofstream trace3("trace3.txt");

    double newptA[3];
    double phi;
    std::vector<double> Bfield;
    std::vector<double> polarPos(3); 
    std::vector<double> newPt(3);

    ////////// Particle tracking

    int nS;
    // Plasma facing surface of HCLL first wall structure
    moab::Range targetFacets; // range containing all of the triangles in the surface of interest
     
    // can specify particular surfaces of interest 
    EntityHandle targetSurf; // surface of interest
    targetSurf = DAG->entity_by_id(2, 2143); // front facing surface of HCLL


    DAG->moab_instance()->list_entity(targetSurf); // list geometric information about surface
    DAG->moab_instance()->get_entities_by_type(targetSurf, MBTRI, targetFacets);
    surfaceIntegrator integrator(targetFacets);
    EntityHandle hit;

    std::vector<double> triA(3), triB(3), triC(3);
    std::vector<double> randTri;
    double fieldDir[3];
    double Bn; // B.n at surface of geometry 
    ds = 0.01;
    nS = 100000;
    int iteration_count = 0; 
    int trace_count = 0;

    int zSign;
    int terminationType;
    for (auto i:targetFacets)
    {
      DAG->next_vol(Surfs[0], vol_h, vol_h);
      iteration_count +=1;
      moab::Range HCLLverts;
      DAG->moab_instance()->get_adjacencies(&i, 1, 0, false, HCLLverts);
      std::vector<double> HCLLcoords(9);
      DAG->moab_instance()->get_coords(HCLLverts, &HCLLcoords[0]);
     
      for (int j=0; j<3; j++)
      {
        triA[j] = HCLLcoords[j];
        triB[j] = HCLLcoords[j+3];
        triC[j] = HCLLcoords[j+6];
      }

      triSource Tri(triA, triB, triC);
      randTri = Tri.random_pt();

      double triStart[3];

      
      

      trace1 << randTri[0] << " " << randTri[1] << " " << randTri[2] << std::endl;
      trace3 << randTri[0] << " " << randTri[1] << " " << randTri[2] << std::endl;

      triStart[0] = randTri[0];
      triStart[1] = randTri[1];
      triStart[2] = randTri[2];

      if (triStart[2] > 0) 
      {
        zSign = 1;
        terminationType = 2;
      }
      else 
      {
        zSign = -1;
        terminationType = 2;
      }

      Bfield = EquData.b_field(randTri, "cart");

      if (Bfield[0] == 0 && Bfield[1] == 0 && Bfield[2] == 0)
      {
        LOG_INFO << "Position of triangle start not in magnetic field. Skipping to next triangle";
        continue;
      }

      polarPos = coordTfm::cart_to_polar(randTri,"forwards");
      Bfield = EquData.b_field_cart(Bfield, polarPos[2], 0);
      norm = sqrt(pow(Bfield[0],2) + pow(Bfield[1],2) + pow(Bfield[2],2));
      normB[0] = Bfield[0]/norm;
      normB[1] = Bfield[1]/norm; 
      normB[2] = Bfield[2]/norm;

      Bn = dot_product(Bfield,Tri.normal);
      if (Bn < 0)
      {
        fieldDir[0] = -normB[0];
        fieldDir[1] = -normB[1];
        fieldDir[2] = -normB[2];

      }
      else if (Bn > 0)
      {
        fieldDir[0] = normB[0];
        fieldDir[1] = normB[1];
        fieldDir[2] = normB[2];
      }

      // DAG->ray_fire(vol_h, triStart, normB, next_surf, next_surf_dist, &history, ds, 1);
      // std::cout << Bfield << std::endl;

      DAG->ray_fire(vol_h, triStart, fieldDir, next_surf, next_surf_dist, &history, ds, 1);
      if (next_surf != 0)
        {
          DAG->next_vol(next_surf, vol_h, vol_h);
          history.get_last_intersection(hit);
          integrator.count_hit(hit);
          //history.reset();
          continue;
        }
      for (int i=0; i<3; ++i)
      {
        newPt[i] = triStart[i] + fieldDir[i]*ds;
      }
      trace1 << newPt[0] << " " << newPt[1] << " " << newPt[2] << std::endl;
      history.rollback_last_intersection();

      for (int i=0; i < nS; ++i, s+ds)
      {
        newptA[0] = newPt[0]; 
        newptA[1] = newPt[1];
        newptA[2] = newPt[2];
        DAG->ray_fire(vol_h, newptA, fieldDir, next_surf, next_surf_dist, &history, ds, 1);
        for (int i=0; i<3; ++i)
        {
          newPt[i] = newPt[i] + fieldDir[i]*ds;
        }

        trace1 << newPt[0] << " " << newPt[1] << " " << newPt[2] << std::endl;

        if (next_surf != 0)
        {
          DAG->next_vol(next_surf, vol_h, vol_h);
          history.get_last_intersection(hit);
          integrator.count_hit(hit);
          //history.reset();
          break;
        }
        history.rollback_last_intersection();

        
        Bfield = EquData.b_field(newPt, "cart");

        if (Bfield[0] == 0)
        {
          LOG_INFO << "TRACE STOPPED BECAUSE LEAVING MAGNETIC FIELD";
          break;
        }
        else
        {
          polarPos = coordTfm::cart_to_polar(newPt,"forwards");
          Bfield = EquData.b_field_cart(Bfield, polarPos[2], 0);
          norm = sqrt(pow(Bfield[0],2) + pow(Bfield[1],2) + pow(Bfield[2],2)); 
          normB[0] = Bfield[0]/norm;
          normB[1] = Bfield[1]/norm; 
          normB[2] = Bfield[2]/norm;

          if (Bn < 0)
          {
            fieldDir[0] = -normB[0];
            fieldDir[1] = -normB[1];
            fieldDir[2] = -normB[2];

          }
          else if (Bn > 0)
          {
            fieldDir[0] = normB[0];
            fieldDir[1] = normB[1];
            fieldDir[2] = normB[2];
          }    
        }


        if (zSign == 1)
        {
          if (newPt[2] < EquData.zcen) // break when crossing zcen
          {
            break;
          } 
        }
        else if (zSign == -1)
        {
          if (newPt[2] > EquData.zcen)
          {
            break;
          }
        }
      }
            


    }

    LOG_WARNING << "Number of rays launched = " << targetFacets.size();
    LOG_WARNING << "Number of ray-facet intersections = " << integrator.raysHit;
    integrator.facet_values(integrator.nRays);


//////////



  }
  else // No runcase specified
  {
    LOG_FATAL << "No runcase specified - please choose from 'specific', 'rayqry' or 'eqdsk'";
  }

  return 0;
}

// Updates origin array for use in ray_fire
void next_pt(double prev_pt[3], double origin[3], double next_surf_dist,
            double dir[3], std::ofstream &ray_intersect){
  // prev_pt -> Array of previous ray launch point
  // next_surf_dist -> Input distance to next surface (from ray_fire)
  // dir -> Array of direction vector of ray
  // write_stream -> filestream to write out ray_intersection points



  for (int i=0; i<3; ++i) { // loop to calculate next ray launch point
    origin[i] = prev_pt[i] + (next_surf_dist * dir[i]);
    ray_intersect << origin[i] << ' ';
  }
  ray_intersect ;
  return;
}

double * vecNorm(double vector[3]){
  static double normalised_vector[3];
  double vector_mag;

  vector_mag = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
  for (int i=0; i<3; i++){
    normalised_vector[i] = vector[i]/vector_mag;
  }
  return normalised_vector;
}


// void next_dir(double raydir[3], const int qrymax, double dir_array[][3]){
//   for (int j=0; j<qrymax; j++)
//   {
//     for (int k=0; k<3; k++)
//   {
//     raydir[k] = dir_array[j+1][k] - dir_array[j][k];
//     //LOG_WARNING << raydir[k];
//   }
//   //LOG_WARNING ;
//   }
//   return;
// }


double dot_product(double vector_a[], double vector_b[]){
   double product = 0;
   for (int i = 0; i < 3; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}

double dot_product(std::vector<double> vector_a, std::vector<double> vector_b){
   double product = 0;
   for (int i = 0; i < 3; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}


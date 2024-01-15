#include <stdio.h>
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
#include <time.h>
#include <any>

#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyLine.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkCompositeDataSet.h>
#include <vtkInformation.h>
#include <vtkSTLReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkAppendFilter.h>

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
#include "vtkAegis.h"
#include "particle.h"


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



  clock_t start = clock();
  settings settings;
  settings.load_params();
  settings.print_params();
 
  LOG_WARNING << "h5m Faceted Geometry file to be used = " << settings.sValues["DAGMC_input"];

  static const char* dagmc_input_file = settings.sValues["DAGMC_input"].c_str();
  static const char* vtk_input_file = settings.sValues["VTK_input"].c_str();

  static const char* ray_qry_exps = settings.sValues["ray_qry"].c_str();
  std::string eqdsk_file = settings.sValues["eqdsk_file"];
  moab::Range Surfs, Vols, Facets, Facet_vertices, Vertices;
  DAG = new DagMC(); // New DAGMC instance
  DAG->load_file(dagmc_input_file); // open test dag file
  DAG->init_OBBTree(); // initialise OBBTree
  DAG->setup_geometry(Surfs, Vols);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, Facets);
  LOG_WARNING << "No of triangles in geometry " << Facets.size();

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
    std::vector<double> coords(3*verts.size());
    DAG->moab_instance()->get_coords(verts, &coords[0]);
  }

  EntityHandle prev_surf; // previous surface id
  EntityHandle next_surf; // surface id
  double next_surf_dist=0.0; // distance to the next surface ray will intersect
  double prev_surf_dist; // previous next_surf_dist before the next ray_fire call
  DagMC::RayHistory history; // initialise RayHistory object
  EntityHandle vol_h = DAG->entity_by_index(3, 1);




  // --------------------------------------------------------------------------------------------------------

  if (settings.sValues["runcase"]=="eqdsk") // structure in place in case I want to have completely different run options
  {
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "--------------------READING EQDSK---------------------" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    equData EquData;
    std::cout << "RMOVE = "  << settings.dValues["rmove"] << std::endl; 
    EquData.read_eqdsk(eqdsk_file);
    EquData.move(settings.dValues["rmove"], settings.dValues["zmove"], settings.dValues["fscale"]); // same values smardda uses for EQ3. TODO - add these parameters to input config file
    EquData.psibdry = settings.dValues["psiref"]; // psiref in geoq.ctl = this 
    EquData.init_interp_splines();
    EquData.gnuplot_out();
    EquData.centre(1);

    double psol = settings.dValues["Psol"];
    double lambda_q = settings.dValues["lambda_q"];

    std::vector<double> vertexCoordinates;
    DAG->moab_instance()->get_vertex_coordinates(vertexCoordinates);
    LOG_WARNING << "NUMBER OF VERTICES = " <<  vertexCoordinates.size() << std::endl;
    int numVertex = vertexCoordinates.size()/3;
    std::vector<std::vector<double>> vertexList(numVertex, std::vector<double> (3));

    // create list of vertexes as XYZ vectors
    for (int i = 0; i<numVertex; i++)
    {
      vertexList[i][0] = vertexCoordinates[i];
      vertexList[i][1] = vertexCoordinates[i + numVertex];
      vertexList[i][2] = vertexCoordinates[i + numVertex*2];
    }

    EquData.psi_limiter(vertexList);
    
    bool plotRZ = true;
    bool plotXYZ = true;
    EquData.write_bfield(plotRZ, plotXYZ);
    std::vector<double> cartPosSource(3);
    std::ofstream launchPosTxt("launch_positions.txt");
    LOG_WARNING << "EquData.psibdry = " << EquData.psibdry;

    double phi;
    std::vector<double> Bfield;
    std::vector<double> polarPos(3);
    std::vector<double> newPt(3);

    // Get triangles from the surface(s) of interest
    moab::Range targetFacets, targetFacetsMeshset; // range containing all of the triangles in the surface of interest
    moab::Range tempFacets;
    // can specify particular surfaces of interest
    EntityHandle targetSurf; // surface of interest
    targetSurf = Surfs[0];

    std::unordered_map<EntityHandle, moab::Range> surfFacets;
    moab::Range surfFacetsKeys;
    int numFacets = 0;
    if (settings.vValues["surfs"][0] != 0)
    {
      int keyctr = 0;
      for (auto &s:settings.vValues["surfs"])
      {
        targetSurf = DAG->entity_by_id(2,s); // surface of interest
        DAG->moab_instance()->get_entities_by_type(targetSurf, MBTRI, surfFacets[targetSurf]);
        std::cout << "Surface ID [" << DAG->get_entity_id(s) << "] " << "No. of elements " << surfFacets[s].size() << std::endl;
        surfFacetsKeys.insert(targetSurf);
        numFacets += surfFacets[targetSurf].size();
        targetFacets.merge(surfFacets[targetSurf]);         
      }
      LOG_WARNING << "Surface IDs provided. Launching from surfaces given by global IDs:";
      for (auto i:settings.vValues["surfs"])
      {LOG_WARNING << "Surface ID - " << i << " ";}
    }
    else
    {
      targetFacets = Facets;
      for (auto &s:Surfs)
      {
        DAG->moab_instance()->get_entities_by_type(s, MBTRI, surfFacets[s]);
        surfFacetsKeys.insert(s);
        std::cout << "Surface ID [" << DAG->get_entity_id(s) << "] " << "No. of elements " << surfFacets[s].size() << std::endl;
        numFacets += surfFacets[s].size();             
        targetFacets.merge(surfFacets[s]);  
      }
      LOG_WARNING << "No surface ID provided. Launching from all surfaces by default. WARNING - Will take a significant amount of time for larger geometry sets";
    }
    LOG_WARNING << "Total Number of Triangles rays launched from = "
                << numFacets;

    surfaceIntegrator integrator(targetFacets);
    EntityHandle hit;

    std::vector<double> triA(3), triB(3), triC(3); // triangle nodes
    double Bn; // B.n at surface of geometry
    double s; 
    double ds = settings.dValues["dsTrack"];
    int nS = settings.iValues["nTrack"];

    int iteration_count = 0;
    int trace_count = 0;
    double coords[9];
    int zSign;

    std::string launchType;
    if (settings.sValues["launchPos"] == "fixed")
    {
      launchType = "fixed";
      LOG_WARNING << "Launching from triangle centroid/barycentre";
    }
    else
    {
      launchType = "random";
      LOG_WARNING << "Launching from random positions in triangles";
    }
    // INITIALISE VTK STUFF -------------------------------------------------

    vtkAegis aegisVTK;

    vtkNew<vtkMultiBlockDataSet> multiBlockRoot, multiBlockBranch;
    std::map<std::string, vtkNew<vtkMultiBlockDataSet>> vtkParticleTracks;

    const char* branchShadowedPart = "Shadowed Particles";
    const char* branchLostPart = "Lost Particles";
    const char* branchDepositingPart = "Depositing Particles";
    const char* branchMaxLengthPart = "Max Length Particles";

    multiBlockRoot->SetBlock(0, multiBlockBranch); // set block 
    multiBlockRoot->GetMetaData(static_cast<int>(0)) // name block
                  ->Set(vtkCompositeDataSet::NAME(), "Particle Tracks");
    LOG_INFO << "Initialising particle_tracks root ";

    // create arrays for various cell data  
    aegisVTK.new_vtkArray("Q", 1);
    aegisVTK.new_vtkArray("B.n_direction", 1);
    aegisVTK.new_vtkArray("Normal", 3);
    aegisVTK.new_vtkArray("B_field", 3);
    aegisVTK.new_vtkArray("Psi_Start", 1);
    aegisVTK.new_vtkArray("B.n", 1);

    int facetCounter=0;

    std::ofstream psi_values("psi_values.txt");
    std::string particleTrace = settings.sValues["trace"];
    double rOutrBdry = settings.dValues["rOutrBdry"];
    LOG_WARNING << "rOutrBdry (from input file) = " << rOutrBdry;
    // LOOP OVER FACETS OF INTEREST IN SURFACE s -----------------------------------------------------------    
    for (auto i:targetFacets)
    {
      //DAG->moab_instance()->list_entity(s);
      facetCounter +=1;
      iteration_count +=1;
      moab::Range HCLLverts;
      std::vector<EntityHandle> verts;
      DAG->moab_instance()->get_adjacencies(&i, 1, 0, false, verts);
      DAG->moab_instance()->get_coords(&verts[0], verts.size(), coords);

      particleBase particle;

      for (int j=0; j<3; j++)
      {
        triA[j] = coords[j];
        triB[j] = coords[j+3];
        triC[j] = coords[j+6];
      }
      triSource Tri(triA, triB, triC, i);

      aegisVTK.arrays["Normal"]->InsertNextTuple3(Tri.unitNormal[0], Tri.unitNormal[1], Tri.unitNormal[2]);
      
      if (launchType == "fixed") // get launchpos on triangle
      {
        particle.set_pos(Tri.centroid());
      }
      else
      {
        particle.set_pos(Tri.random_pt());
      }
      integrator.launchPositions[i] = particle.launchPos;

      particle.check_if_in_bfield(EquData);
      if (particle.outOfBounds)
      {
        LOG_INFO << "Particle start is out of magnetic field bounds. Skipping to next triangle. Check correct eqdsk is being used for the given geometry";
        aegisVTK.arrays["B_field"]->InsertNextTuple3(0,0,0);
        aegisVTK.arrays["Psi_Start"]->InsertNextTuple1(0);
        aegisVTK.arrays["B.n"]->InsertNextTuple1(0);
        aegisVTK.arrays["B.n_direction"]->InsertNextTuple1(0);
        aegisVTK.arrays["Q"]->InsertNextTuple1(0.0);
        continue;
      }

      particle.set_dir(EquData); // Set unit direction vector along cartesian magnetic field vector

      vtkNew<vtkPoints> vtkpoints;
      int vtkPointCounter = 0;

      if (particleTrace == "yes")
      {
        vtkpoints->InsertNextPoint(particle.launchPos[0], particle.launchPos[1], particle.launchPos[2]);
        vtkPointCounter +=1;
      }

      std::string forwards = "forwards";
      polarPos = coordTfm::cart_to_polar(particle.launchPos, forwards);

      aegisVTK.arrays["B_field"]->InsertNextTuple3(particle.BfieldXYZ[0], particle.BfieldXYZ[1], particle.BfieldXYZ[2]);
      Bn = dot_product(particle.BfieldXYZ,Tri.unitNormal);
      
      // CALCULATING Q at surface 
      double psi;

      std::vector<double> temp;
      temp = coordTfm::cart_to_polar(particle.launchPos, "forwards");
      temp = coordTfm::polar_to_flux(temp, "forwards", EquData);
      psi = temp[0];
      double psid = psi + EquData.psibdry;
      double Q;
      aegisVTK.arrays["Psi_Start"]->InsertNextTuple1(psi);
      aegisVTK.arrays["B.n"]->InsertNextTuple1(Bn);
      // --------------------------- TODO move Bn check and Q calculation into a class (maybe triangle source class?)
      particle.align_dir_to_surf(Bn);
      if (Bn < 0)
      {
        aegisVTK.arrays["B.n_direction"]->InsertNextTuple1(-1.0);
      }
      else if (Bn > 0)
      {
        aegisVTK.arrays["B.n_direction"]->InsertNextTuple1(1.0);
        
      }
      Q = EquData.omp_power_dep(psid, psol, lambda_q, Bn, "exp");
      Q = std::fabs(Q);
      psi_values << Q << " " << psi << " " << psid << " " << particle.BfieldXYZ[0] << " " << particle.BfieldXYZ[1] 
                 << " " << particle.BfieldXYZ[2] << " " << Bn << std::endl;


      // Don't intersect on initial launch 
      int ray_orientation = 1;
      history.reset();
      DAG->ray_fire(vol_h, particle.launchPos.data(), particle.dir.data(), next_surf, next_surf_dist, &history,ds,ray_orientation);
      if (next_surf != 0) 
      {
        history.get_last_intersection(hit);
        history.rollback_last_intersection();
        DAG->next_vol(next_surf, vol_h, vol_h);
        LOG_INFO << "---- RAY HIT ON LAUNCH [" << facetCounter << "] ----";
      }

      for (int j=0; j<3; ++j)
      {
        newPt[j] = particle.launchPos[j] + particle.dir[j]*ds;
        particle.pos[j] = particle.pos[j] + particle.dir[j]*ds;
      }

      if (particleTrace == "yes")
      {
        vtkpoints->InsertNextPoint(newPt[0], newPt[1], newPt[2]);
        vtkPointCounter +=1;
      }
      //history.rollback_last_intersection();
      history.reset();
      s = ds;
      bool traceEnded = false; 

      // --------------------- LOOP OVER PARTICLE TRACK --------------------- 
      for (int j=0; j < nS; ++j) // loop along fieldline
      {
        history.reset();
        traceEnded = false;
        s += ds;
        DAG->ray_fire(vol_h, particle.pos.data(), particle.dir.data(), next_surf, next_surf_dist, &history, ds, ray_orientation);

        if (next_surf != 0) // Terminate fieldline trace since shadowing surface hit
        {
          newPt[0] = newPt[0] +particle.dir[0]*next_surf_dist;
          newPt[1] = newPt[1] +particle.dir[1]*next_surf_dist;
          newPt[2] = newPt[2] +particle.dir[2]*next_surf_dist;

          //DAG->next_vol(next_surf, vol_h, vol_h);
          history.get_last_intersection(hit);
          integrator.count_hit(hit);
          LOG_INFO << "Surface " << next_surf << " hit after travelling " << s << " units";
          history.rollback_last_intersection();
          history.reset();

          integrator.store_heat_flux(i, 0.0);
          if (particleTrace == "yes")
          {
            vtkpoints->InsertNextPoint(newPt[0], newPt[1], newPt[2]);
            vtkPointCounter +=1;
            // if block does not exist, create it
            if (vtkParticleTracks.find(branchShadowedPart) == vtkParticleTracks.end())
            {
              int staticCast = aegisVTK.multiBlockCounters.size();
              multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchShadowedPart]); // set block 
              multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                              ->Set(vtkCompositeDataSet::NAME(), branchShadowedPart); 
              std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchShadowedPart << std::endl;
              aegisVTK.multiBlockCounters[branchShadowedPart] = 0;
            }  
            vtkSmartPointer<vtkPolyData> polydataTrack;
            polydataTrack = aegisVTK.new_track(branchShadowedPart, vtkpoints, 0.0);
            vtkParticleTracks[branchShadowedPart]->SetBlock(aegisVTK.multiBlockCounters[branchShadowedPart], polydataTrack);
            
          }
          aegisVTK.arrays["Q"]->InsertNextTuple1(0.0);
          traceEnded = true;
          break; // break loop over ray if surface hit
        }
        else 
        {
          newPt[0] = newPt[0] +particle.dir[0]*ds;
          newPt[1] = newPt[1] +particle.dir[1]*ds; // ds because a surface isnt reached so next_surf_dist does not have a value
          newPt[2] = newPt[2] +particle.dir[2]*ds;
        }
        
        if (particleTrace == "yes")
        {
          vtkpoints->InsertNextPoint(newPt[0], newPt[1], newPt[2]);
          vtkPointCounter +=1;    
        }

        particle.set_pos(newPt);
        particle.check_if_in_bfield(EquData);
        if (particle.outOfBounds) // Terminate fieldline trace since particle has left magnetic domain
        {
          LOG_INFO << "TRACE STOPPED BECAUSE LEAVING MAGNETIC FIELD";
          integrator.count_lost_ray();
          integrator.store_heat_flux(i, 0.0);
          if (particleTrace == "yes")
          {
            if (vtkParticleTracks.find(branchLostPart) == vtkParticleTracks.end())
            {
              int staticCast = aegisVTK.multiBlockCounters.size();
              multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchLostPart]); // set block 
              multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                              ->Set(vtkCompositeDataSet::NAME(), branchLostPart); 
              std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchLostPart << std::endl;
              aegisVTK.multiBlockCounters[branchLostPart] = 0;
            }  
            vtkSmartPointer<vtkPolyData> polydataTrack;
            polydataTrack = aegisVTK.new_track(branchLostPart, vtkpoints, 0.0);
            vtkParticleTracks[branchLostPart]->SetBlock(aegisVTK.multiBlockCounters[branchLostPart], polydataTrack);
          }
          aegisVTK.arrays["Q"]->InsertNextTuple1(0.0);
          traceEnded = true;
          break; // break if ray leaves magnetic field
        }
        else
        {
          particle.set_dir(EquData);
          particle.align_dir_to_surf(Bn);
        }

        particle.check_if_midplane_reached(EquData.zcen, EquData.rbdry, rOutrBdry);

        if (particle.atMidplane !=0) // if particle is 1 or 2 then deposit heat
        {
          polarPos = coordTfm::cart_to_polar(newPt, "forwards");
          std::vector<double> fluxPos = coordTfm::polar_to_flux(particle.get_pos("polar"), "forwards", EquData);
          integrator.store_heat_flux(i, Q);
          if (particleTrace == "yes")
          {
            if (vtkParticleTracks.find(branchDepositingPart) == vtkParticleTracks.end())
            {
              int staticCast = aegisVTK.multiBlockCounters.size();
              multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchDepositingPart]); // set block 
              multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                              ->Set(vtkCompositeDataSet::NAME(), branchDepositingPart); 
              std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchDepositingPart << std::endl;
              aegisVTK.multiBlockCounters[branchDepositingPart] = 0;
            }  
            vtkSmartPointer<vtkPolyData> polydataTrack;
            polydataTrack = aegisVTK.new_track(branchDepositingPart, vtkpoints, Q);
            vtkParticleTracks[branchDepositingPart]->SetBlock(aegisVTK.multiBlockCounters[branchDepositingPart], polydataTrack);
          }
          aegisVTK.arrays["Q"]->InsertNextTuple1(Q);
          traceEnded = true;
          break; // break if ray hits omp
        }

      }
      
      if (traceEnded == false)
      {
        if (vtkParticleTracks.find(branchMaxLengthPart) == vtkParticleTracks.end())
        {
          int staticCast = aegisVTK.multiBlockCounters.size();
          multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchMaxLengthPart]); // set block 
          multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                          ->Set(vtkCompositeDataSet::NAME(), branchMaxLengthPart); 
          std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchMaxLengthPart << std::endl;
          aegisVTK.multiBlockCounters[branchMaxLengthPart] = 0;
        }  
        vtkSmartPointer<vtkPolyData> polydataTrack;
        polydataTrack = aegisVTK.new_track(branchMaxLengthPart, vtkpoints, 0.0);
        vtkParticleTracks[branchMaxLengthPart]->SetBlock(aegisVTK.multiBlockCounters[branchMaxLengthPart], polydataTrack);
        aegisVTK.arrays["Q"]->InsertNextTuple1(0.0);
        LOG_INFO << "Fieldline trace reached maximum length before intersection";
        traceEnded = true;
      }
    }

    LOG_WARNING << "Number of rays launched = " << facetCounter;
    LOG_WARNING << "Number of shadowed ray intersections = " << integrator.raysHit;
    LOG_WARNING << "Number of rays depositing power from omp = " << integrator.raysHeatDep;
    LOG_WARNING << "Number of rays lost from magnetic domain = " << integrator.raysLost;
    LOG_WARNING << "Number of rays that reached the maximum allowed length = " 
                << aegisVTK.multiBlockCounters[branchMaxLengthPart];

    integrator.csv_out(integrator.powFac);

    integrator.piecewise_multilinear_out(integrator.powFac);

    aegisVTK.write_unstructuredGrid(vtk_input_file, "out.vtk");

    if (particleTrace == "yes")
    {  
      vtkNew<vtkXMLMultiBlockDataWriter> vtkMBWriter;
      vtkMBWriter->SetFileName("particle_tracks.vtm");
      vtkMBWriter->SetInputData(multiBlockRoot);
      vtkMBWriter->Write();
    }

  }

  else // No runcase specified
  {
    LOG_FATAL << "No runcase specified - please choose from 'eqdsk'";
  }


  clock_t end = clock();
  double elapsed = double(end - start)/CLOCKS_PER_SEC;

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "Elapsed Aegis run time = " << elapsed << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

  return 0;
}



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


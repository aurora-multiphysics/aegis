#include "aegisClass.h"


// void Aegis::setup_dagmc(
  
// )


// void Aegis::setup_equillibrium(){
  
// }

void AegisClass::Execute(){

  clock_t start = clock(); // start clock time

  // set parameters
  runSettings.load_params(); 
  runSettings.print_params();
  dagmcInputFile = runSettings.sValues["DAGMC_input"];
  vtkInputFile = runSettings.sValues["VTK_input"];
  eqdskInputFile = runSettings.sValues["eqdsk_file"];
  powerSOL = runSettings.dValues["Psol"];
  lambdaQ = runSettings.dValues["lambda_q"];
  trackStepSize = runSettings.dValues["dsTrack"];
  maxTrackSteps = runSettings.iValues["nTrack"];
  particleLaunchPos = runSettings.sValues["launchPos"];

  if (particleLaunchPos == "fixed")
  {
    LOG_WARNING << "Launching from triangle centroid/barycentre";
  }
  else
  {
    LOG_WARNING << "Launching from random positions in triangles";
  }

  userROutrBdry = runSettings.dValues["rOutrBdry"];
  LOG_WARNING << "User set R value at Outer midplane = " << userROutrBdry << std::endl;

  drawParticleTracks = runSettings.sValues["trace"];
  rmove = runSettings.dValues["rmove"];
  zmove = runSettings.dValues["zmove"];
  fscale = runSettings.dValues["fscale"];
  psiref = runSettings.dValues["psiref"];

  // setup dagmc instance
  DAG = std::make_unique<moab::DagMC>();
  DAG->load_file(dagmcInputFile.c_str());
  DAG->init_OBBTree();
  DAG->setup_geometry(surfsList, volsList);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, facetsList);
  numFacets = facetsList.size();
  std::cout << "Number of Triangles in Geometry " << numFacets << std::endl;
  volID = DAG->entity_by_index(3,1);

  // setup B Field data
  bFieldData.read_eqdsk(eqdskInputFile);
  bFieldData.move(rmove, zmove, fscale);
  bFieldData.psibdry = psiref; // abstract out to function
  bFieldData.init_interp_splines();
  bFieldData.centre(1);

  std::vector<double> vertexCoordinates;
  DAG->moab_instance()->get_vertex_coordinates(vertexCoordinates);
  LOG_WARNING << "NUMBER OF NODES = " <<  vertexCoordinates.size() << std::endl;
  int numNodes = vertexCoordinates.size()/3;
  std::vector<std::vector<double>> vertexList(numNodes, std::vector<double> (3));

  for (int i = 0; i<numNodes; i++)
  {
    vertexList[i][0] = vertexCoordinates[i];
    vertexList[i][1] = vertexCoordinates[i + numNodes];
    vertexList[i][2] = vertexCoordinates[i + numNodes*2];
  }

  bFieldData.psi_limiter(vertexList);
  bFieldData.write_bfield(plotBFieldRZ, plotBFieldXYZ);


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


  // Initialise VTK  
  vtkInterface.new_vtkArray("Q", 1);
  vtkInterface.new_vtkArray("B.n_direction", 1);
  vtkInterface.new_vtkArray("Normal", 3);
  vtkInterface.new_vtkArray("B_field", 3);
  vtkInterface.new_vtkArray("Psi_Start", 1);
  vtkInterface.new_vtkArray("B.n", 1);


  double phi = 0.0;
  std::vector<double> Bfield; 
  std::vector<double> polarPos(3);
  std::vector<double> newPt(3);
  
  std::vector<double> triCoords(9);
  std::vector<double> triA(3), triB(3), triC(3);
  integrator = std::make_unique<surfaceIntegrator>(facetsList);
  // loop over all facets
  int facetCounter = 0;
  std::ofstream aegisOUT("aegis_new_OUT.txt");

  traceEnded = false;

  for (auto facet:facetsList){
    facetCounter +=1;
    particleBase particle;

    std::vector<moab::EntityHandle> triNodes;
    DAG->moab_instance()->get_adjacencies(&facet, 1, 0, false, triNodes);
    DAG->moab_instance()->get_coords(&triNodes[0], triNodes.size(), triCoords.data());

    for (int j=0; j<3; j++){
      triA[j] = triCoords[j];
      triB[j] = triCoords[j+3];
      triC[j] = triCoords[j+6];
    }
    triSource Tri(triA, triB, triC, facet); 
    vtkInterface.arrays["Normal"]->InsertNextTuple3(Tri.unitNormal[0], Tri.unitNormal[1], Tri.unitNormal[2]);

    if (particleLaunchPos == "fixed"){
      particle.set_pos(Tri.centroid());
    }
    else{
      particle.set_pos(Tri.random_pt());
    }
    integrator->launchPositions[facet] = particle.get_pos("cart");
    
    particle.check_if_in_bfield(bFieldData); // if out of bounds skip to next triangle
    if (particle.outOfBounds){
      LOG_INFO << "Particle start is out of magnetic field bounds. Skipping to next triangle. Check correct eqdsk is being used for the given geometry";
      vtkInterface.arrays["B_field"]->InsertNextTuple3(0,0,0);
      vtkInterface.arrays["Psi_Start"]->InsertNextTuple1(0);
      vtkInterface.arrays["B.n"]->InsertNextTuple1(0);
      vtkInterface.arrays["B.n_direction"]->InsertNextTuple1(0);
      vtkInterface.arrays["Q"]->InsertNextTuple1(0.0);
      continue;
    }

    vtkNew<vtkPoints> vtkpoints;
    int vtkPointCounter = 0;

    vtkpoints->InsertNextPoint(particle.launchPos[0], particle.launchPos[1], particle.launchPos[2]);
    vtkPointCounter +=1;

    particle.set_dir(bFieldData);
    polarPos = coordTfm::cart_to_polar(particle.launchPos, "forwards");

    vtkInterface.arrays["B_field"]->InsertNextTuple3(particle.BfieldXYZ[0], particle.BfieldXYZ[1], particle.BfieldXYZ[2]);
    BdotN = Tri.dot_product(particle.BfieldXYZ);
    psi = particle.get_psi(bFieldData); 
    psid = psi + bFieldData.psibdry; 

    vtkInterface.arrays["Psi_Start"]->InsertNextTuple1(psi);
    vtkInterface.arrays["B.n"]->InsertNextTuple1(BdotN);
    particle.align_dir_to_surf(BdotN);

    if (BdotN < 0){
      vtkInterface.arrays["B.n_direction"]->InsertNextTuple1(-1.0);
    }
    else if (BdotN > 0){
      vtkInterface.arrays["B.n_direction"]->InsertNextTuple1(1.0);
    }
    Q = bFieldData.omp_power_dep(psid, powerSOL, lambdaQ, BdotN, "exp");
    double psiOnSurface = psi;
    
    
    // Start ray tracing
    DAG->ray_fire(volID, particle.launchPos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history,trackStepSize,rayOrientation);
    if (nextSurf != 0) 
    {
      history.get_last_intersection(intersectedFacet);
      history.rollback_last_intersection();
      DAG->next_vol(nextSurf, volID, volID);
      LOG_INFO << "---- RAY HIT ON LAUNCH [" << facetCounter << "] ----";
    }
    history.reset();
    double trackLength = trackStepSize;


    for (int j=0; j<3; ++j)
    {
      newPt[j] = particle.launchPos[j] + particle.dir[j]*trackStepSize;
    }

    vtkpoints->InsertNextPoint(newPt[0], newPt[1], newPt[2]);
    vtkPointCounter +=1;


    // loop along particle track
    for (int j=0; j<maxTrackSteps; ++j){
      history.reset();
      trackLength += trackStepSize;
      DAG->ray_fire(volID, particle.pos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history, trackStepSize, rayOrientation);

      if (nextSurf != 0){
        newPt[0] = particle.pos[0] + particle.dir[0] * nextSurfDist;
        newPt[1] = particle.pos[1] + particle.dir[1] * nextSurfDist;
        newPt[2] = particle.pos[2] + particle.dir[2] * nextSurfDist;

        history.get_last_intersection(intersectedFacet);
        integrator->count_hit(intersectedFacet);
        LOG_INFO << "Surface " << nextSurf << " hit after travelling " << trackLength << " units";
        history.rollback_last_intersection();
        history.reset();

        vtkpoints->InsertNextPoint(newPt[0], newPt[1], newPt[2]);
        vtkPointCounter +=1;
        // if block does not exist, create it
        if (vtkParticleTracks.find(branchShadowedPart) == vtkParticleTracks.end())
        {
          int staticCast = vtkInterface.multiBlockCounters.size();
          multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchShadowedPart]); // set block 
          multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                          ->Set(vtkCompositeDataSet::NAME(), branchShadowedPart); 
          std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchShadowedPart << std::endl;
          vtkInterface.multiBlockCounters[branchShadowedPart] = 0;
        }  
        vtkNew<vtkPolyData> polydataTrack;
        polydataTrack = vtkInterface.new_track(branchShadowedPart, vtkpoints, 0.0);
        vtkParticleTracks[branchShadowedPart]->SetBlock(vtkInterface.multiBlockCounters[branchShadowedPart], polydataTrack);


        vtkInterface.arrays["Q"]->InsertNextTuple1(0.0);
        integrator->store_heat_flux(facet,0.0);
        psiQ_values.push_back(std::make_pair(psiOnSurface,Q));
        traceEnded = true;
        break;
      }
      else{
        newPt[0] = particle.pos[0] + particle.dir[0] * trackStepSize;
        newPt[1] = particle.pos[1] + particle.dir[1] * trackStepSize;
        newPt[2] = particle.pos[2] + particle.dir[2] * trackStepSize;
      }
      
      vtkpoints->InsertNextPoint(newPt[0], newPt[1], newPt[2]);
      vtkPointCounter +=1; 

      particle.set_pos(newPt);
      particle.check_if_in_bfield(bFieldData);
      if (particle.outOfBounds){
        LOG_INFO << "TRACE STOPPED BECAUSE LEAVING MAGNETIC FIELD";
        integrator->count_lost_ray();

        if (vtkParticleTracks.find(branchLostPart) == vtkParticleTracks.end())
        {
          int staticCast = vtkInterface.multiBlockCounters.size();
          multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchLostPart]); // set block 
          multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                          ->Set(vtkCompositeDataSet::NAME(), branchLostPart); 
          std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchLostPart << std::endl;
          vtkInterface.multiBlockCounters[branchLostPart] = 0;
        }  
        vtkNew<vtkPolyData> polydataTrack;
        polydataTrack = vtkInterface.new_track(branchLostPart, vtkpoints, 0.0);
        vtkParticleTracks[branchLostPart]->SetBlock(vtkInterface.multiBlockCounters[branchLostPart], polydataTrack);

        vtkInterface.arrays["Q"]->InsertNextTuple1(0.0);
        integrator->store_heat_flux(facet,0.0);
        psiQ_values.push_back(std::make_pair(psiOnSurface,Q));
        traceEnded = true;
        break;
      }
      else{
        particle.set_dir(bFieldData);
        particle.align_dir_to_surf(BdotN);
      }

      particle.check_if_midplane_reached(bFieldData.zcen, bFieldData.rbdry, userROutrBdry);
      if (particle.atMidplane != 0){
        
        if (vtkParticleTracks.find(branchDepositingPart) == vtkParticleTracks.end())
        {
          int staticCast = vtkInterface.multiBlockCounters.size();
          multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchDepositingPart]); // set block 
          multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                          ->Set(vtkCompositeDataSet::NAME(), branchDepositingPart); 
          std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchDepositingPart << std::endl;
          vtkInterface.multiBlockCounters[branchDepositingPart] = 0;
        }  
        vtkNew<vtkPolyData> polydataTrack;
        polydataTrack = vtkInterface.new_track(branchDepositingPart, vtkpoints, Q);
        vtkParticleTracks[branchDepositingPart]->SetBlock(vtkInterface.multiBlockCounters[branchDepositingPart], polydataTrack);

        vtkInterface.arrays["Q"]->InsertNextTuple1(Q);
        integrator->store_heat_flux(facet,Q);
        psiQ_values.push_back(std::make_pair(psiOnSurface,Q));
        traceEnded = true;
        break;
      }

    } 

    if (traceEnded == false){
      if (vtkParticleTracks.find(branchMaxLengthPart) == vtkParticleTracks.end())
      {
        int staticCast = vtkInterface.multiBlockCounters.size();
        multiBlockBranch->SetBlock(staticCast, vtkParticleTracks[branchMaxLengthPart]); // set block 
        multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                        ->Set(vtkCompositeDataSet::NAME(), branchMaxLengthPart); 
        std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchMaxLengthPart << std::endl;
        vtkInterface.multiBlockCounters[branchMaxLengthPart] = 0;
      }  
      vtkNew<vtkPolyData> polydataTrack;
      polydataTrack = vtkInterface.new_track(branchMaxLengthPart, vtkpoints, 0.0);
      vtkParticleTracks[branchMaxLengthPart]->SetBlock(vtkInterface.multiBlockCounters[branchMaxLengthPart], polydataTrack);
      
      vtkInterface.arrays["Q"]->InsertNextTuple1(0.0);
      integrator->store_heat_flux(facet,0.0);
      psiQ_values.push_back(std::make_pair(psiOnSurface,Q));
      LOG_INFO << "Fieldline trace reached maximum length before intersection";
      traceEnded = true;
    }

  }



  vtkInterface.write_unstructuredGrid(vtkInputFile.c_str(), "out.vtk");

  vtkNew<vtkXMLMultiBlockDataWriter> vtkMBWriter;
  vtkMBWriter->SetFileName("particle_tracks.vtm");
  vtkMBWriter->SetInputData(multiBlockRoot);
  vtkMBWriter->Write();

  integrator->print_particle_stats();

  clock_t end = clock();
  double elapsed = double(end - start)/CLOCKS_PER_SEC;

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "Elapsed Aegis run time = " << elapsed << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

}

// void AegisClass::particle_is_shadowed(EntityHandle facet, particleBase particle) {

//   history.get_last_intersection(intersectedFacet);
//   integrator->count_hit(intersectedFacet);
//   LOG_INFO << "Surface " << nextSurf << " hit after travelling " << trackLength << " units";
//   history.rollback_last_intersection();
//   history.reset();

//   vtkInterface.arrays["Q"]->InsertNextTuple1(0.0);
//   integrator->store_heat_flux(facet,0.0);
//   psiQ_values.push_back(std::make_pair(psiOnSurface,Q));
//   traceEnded = true;
// }

int AegisClass::num_facets(){
  return numFacets;
}

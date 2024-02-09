#include "aegisClass.h"
#include <mpi.h>

void AegisClass::Execute(std::string settingsFile){

  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  settingsFileName = settingsFile;
  init_solve();
  init_geometry();

  vtkInterface->init();  
//  polarDAG = std::make_unique<moab::DagMC>(DAG->moab_instance());

  std::vector<double> Bfield; 
  std::vector<double> polarPos(3);
  std::vector<double> newPt(3);  
  std::vector<double> triCoords(9);
  std::vector<double> triA(3), triB(3), triC(3);
  moab::Range targetSurfaceList = select_target_surface();
  integrator = std::make_unique<surfaceIntegrator>(targetSurfaceList);

  int totalNumberOfFacets = targetSurfaceList.size();
  int numberofFacets = totalNumberOfFacets / nprocs;
  int startFacet = rank * numberofFacets;
  int endFacet = startFacet + numberofFacets;
  nFacets = 0;
  
  for (int i=startFacet; i<endFacet; ++i){ // loop over all facets
    const auto facet = targetSurfaceList[i];
    ++nFacets;
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
    // vtkInterface->insert_next_uStrGrid("Normal", Tri.unitNormal);

    if (particleLaunchPos == "fixed"){
      particle.set_pos(Tri.centroid());
    }
    else{
      particle.set_pos(Tri.random_pt());
    }
    integrator->set_launch_position(facet, particle.get_pos("cart"));
    
    particle.check_if_in_bfield(bFieldData); // if out of bounds skip to next triangle
    if (particle.outOfBounds){
      LOG_INFO << "Particle start is out of magnetic field bounds. Skipping to next triangle. Check correct eqdsk is being used for the given geometry";
      // vtkInterface->insert_next_uStrGrid("B_field", {0.0, 0.0, 0.0});
      // vtkInterface->insert_next_uStrGrid("Psi_Start", 0.0);
      // vtkInterface->insert_next_uStrGrid("B.n", 0.0);
      // vtkInterface->insert_next_uStrGrid("B.n_direction", 0.0);
    //  vtkInterface->insert_next_uStrGrid("Q", 0.0);
    qValues.push_back(0.0);
      continue;
    }
    vtkInterface->init_new_vtkPoints();
    // vtkInterface->insert_next_point_in_track(particle.launchPos);

    particle.set_dir(bFieldData);
    polarPos = coordTfm::cart_to_polar(particle.launchPos, "forwards");

    // vtkInterface->insert_next_uStrGrid("B_field", particle.BfieldXYZ);
    BdotN = Tri.dot_product(particle.BfieldXYZ); 
    psi = particle.get_psi(bFieldData); 

    psiValues.push_back(psi);
    psid = psi + bFieldData.psibdry; 

    // vtkInterface->insert_next_uStrGrid("Psi_Start", psi);
    // vtkInterface->insert_next_uStrGrid("B.n", BdotN);
    particle.align_dir_to_surf(BdotN);

    if (BdotN < 0){
      // vtkInterface->insert_next_uStrGrid("B.n_direction", -1.0);
    }
    else if (BdotN > 0){
      // vtkInterface->insert_next_uStrGrid("B.n_direction", 1.0);
    }
    Q = bFieldData.omp_power_dep(psid, powerSOL, lambdaQ, BdotN, "exp");
    
    psiOnSurface = psi;
    
    // Start ray tracing
    DagMC::RayHistory history;
    ray_hit_on_launch(particle, history);
    trackLength = trackStepSize;

    particle.update_vectors(trackStepSize, bFieldData);
    vtkInterface->insert_next_point_in_track(particle.pos);

    // loop along particle track
    for (int j=0; j<maxTrackSteps; ++j){
      trackLength += trackStepSize;
      DAG->ray_fire(volID, particle.pos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history, trackStepSize, rayOrientation);

      if (nextSurf != 0){
        particle.update_vectors(nextSurfDist); // update position to surface intersection point
        terminate_particle(facet, history, terminationState::SHADOWED);
        break;
      }
      else{
        particle.update_vectors(trackStepSize); // update position by stepsize
      }
      
      vtkInterface->insert_next_point_in_track(particle.pos);
      particle.check_if_in_bfield(bFieldData);
      if (particle.outOfBounds){
        terminate_particle(facet, history, terminationState::LOST);
        break;
      }
      else{
        particle.set_dir(bFieldData);
        particle.align_dir_to_surf(BdotN);
      }
      particle.check_if_midplane_reached(bFieldData.zcen, bFieldData.rbdry, userROutrBdry);
      
      if (particle.atMidplane != 0 && !noMidplaneTermination){ 
        terminate_particle(facet, history, terminationState::DEPOSITING);
        break;
      }

      if (j == (maxTrackSteps-1))
      {
        terminate_particle(facet, history, terminationState::MAXLENGTH);
        break;
      }

    } 

  }
  // write out data and print final

  std::vector<double> allQValues;
  std::vector<double> allPsiValues;
  std::array<int, 4> particleStats = integrator->particle_stats(); 
  std::array<int, 4> totalParticleStats;


  if (rank != 0){
    MPI_Send(psiValues.data(), psiValues.size(), MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
    MPI_Send(qValues.data(), qValues.size(), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);

  }
  else {
    allPsiValues.insert(allPsiValues.end(), psiValues.begin(), psiValues.end());
    allQValues.insert(allQValues.end(), qValues.begin(), qValues.end());
    totalParticleStats = integrator->particle_stats();
  }


  for (int i=1; i<nprocs; ++i){
    if (rank == 0){
      MPI_Recv(psiValues.data(), psiValues.size(), MPI_DOUBLE, i, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      allPsiValues.insert(allPsiValues.end(), psiValues.begin(), psiValues.end());

      MPI_Recv(qValues.data(), qValues.size(), MPI_DOUBLE, i, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      allQValues.insert(allQValues.end(), qValues.begin(), qValues.end());

      MPI_Recv(particleStats.data(), particleStats.size(), MPI_INT, i, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int j=0; j<particleStats.size(); ++j){
        totalParticleStats[j] += particleStats[j]; 
      }
    }
  }

  if (rank == 0){

    for (int i=0; i<num_facets(); ++i){
      vtkInterface->insert_next_uStrGrid("Q", allQValues[i]);
      vtkInterface->insert_next_uStrGrid("Psi_Start", allPsiValues[i]);

    }

    vtkInterface->write_unstructuredGrid(vtkInputFile, "out.vtk");
    if (drawParticleTracks == "yes"){
      vtkInterface->write_multiBlockData("particle_tracks.vtm");
    }
    print_particle_stats(totalParticleStats);
  }

  mpi_particle_stats();

  }


void AegisClass::init_solve(){
  
  clock_t start = clock(); // start clock time

  // set runtime parameters
  runSettings.load_params(settingsFileName); 
  //runSettings.print_params();
  dagmcInputFile = runSettings.sValues["DAGMC_input"];
  vtkInputFile = "target_facets.stl";
  eqdskInputFile = runSettings.sValues["eqdsk_file"];
  powerSOL = runSettings.dValues["Psol"];
  lambdaQ = runSettings.dValues["lambda_q"];
  trackStepSize = runSettings.dValues["dsTrack"];
  maxTrackSteps = runSettings.iValues["nTrack"];
  particleLaunchPos = runSettings.sValues["launchPos"];
  drawParticleTracks = runSettings.sValues["trace"];
  std::string noDep = runSettings.sValues["no_deposition"];
  if (noDep == "yes") {noMidplaneTermination = true;}

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

  rmove = runSettings.dValues["rmove"];
  zmove = runSettings.dValues["zmove"];
  fscale = runSettings.dValues["fscale"];
  psiref = runSettings.dValues["psiref"];


  //initialise memeory for smart pointers
  DAG = std::make_unique<moab::DagMC>();
  vtkInterface = std::make_unique<VtkInterface>(drawParticleTracks);
}

void AegisClass::init_geometry(){
  // setup dagmc instance
  DAG->load_file(dagmcInputFile.c_str());
  DAG->init_OBBTree();
  DAG->setup_geometry(surfsList, volsList);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, facetsList);
  std::cout << "Number of Triangles in Geometry " << facetsList.size() << std::endl;
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
}


void AegisClass::terminate_particle(const moab::EntityHandle &facet, DagMC::RayHistory &history, terminationState termination){
  double heatflux;

  switch(termination){
    case terminationState::DEPOSITING:
      heatflux = Q;
      vtkInterface->write_particle_track(branchDepositingPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);
      LOG_INFO << "Midplane reached. Depositing power";
      break;

    case terminationState::SHADOWED:
      heatflux = 0.0;
      history.get_last_intersection(intersectedFacet);
      vtkInterface->write_particle_track(branchShadowedPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);
      LOG_INFO << "Surface " << nextSurf << " hit after travelling " << trackLength << " units";
      break;

    case terminationState::LOST:
      heatflux = 0.0;   
      vtkInterface->write_particle_track(branchLostPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);    
      LOG_INFO << "TRACE STOPPED BECAUSE LEAVING MAGNETIC FIELD";
      break;

    case terminationState::MAXLENGTH:
      heatflux = 0.0;
      vtkInterface->write_particle_track(branchMaxLengthPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);        
      LOG_INFO << "Fieldline trace reached maximum length before intersection";
  }

//  vtkInterface->insert_next_uStrGrid("Q", heatflux);
  qValues.push_back(heatflux);
  psiQ_values.push_back(std::make_pair(psiOnSurface, Q));
}

void AegisClass::ray_hit_on_launch(particleBase &particle, DagMC::RayHistory &history){
  DAG->ray_fire(volID, particle.launchPos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history,trackStepSize,rayOrientation);
  if (nextSurf != 0) 
  {
    history.get_last_intersection(intersectedFacet);
    DAG->next_vol(nextSurf, volID, volID);
    LOG_INFO << "---- RAY HIT ON LAUNCH [" << facetCounter << "] ----";
  }
}

// Get triangles from the surface(s) of interest
moab::Range AegisClass::select_target_surface(){ 

  moab::Range targetFacets; // range containing all of the triangles in the surface of interest
  // can specify particular surfaces of interest
  EntityHandle targetSurf; // surface of interest
  targetSurf = surfsList[0];

  std::unordered_map<EntityHandle, moab::Range> surfFacets;
  moab::Range surfFacetsKeys;
  int numTargetFacets = 0;
  if (runSettings.vValues["surfs"][0] != 0)
  {
    for (auto &s:runSettings.vValues["surfs"])
    {
      targetSurf = DAG->entity_by_id(2,s); // surface of interest
      DAG->moab_instance()->get_entities_by_type(targetSurf, MBTRI, surfFacets[targetSurf]);
      std::cout << "Surface ID [" << s << "] " << "No. of elements " << surfFacets[targetSurf].size() << std::endl;
      surfFacetsKeys.insert(targetSurf);
      numTargetFacets += surfFacets[targetSurf].size();
      targetFacets.merge(surfFacets[targetSurf]);         
    }

    LOG_WARNING << "Surface IDs provided. Launching from surfaces given by global IDs:";
  }
  else
  {
    targetFacets = facetsList;
    for (auto &s:surfsList)
    {
      DAG->moab_instance()->get_entities_by_type(s, MBTRI, surfFacets[s]);
      surfFacetsKeys.insert(s);
      std::cout << "Surface ID [" << s << "] " << "No. of elements " << surfFacets[s].size() << std::endl;
      numTargetFacets += surfFacets[s].size();             
      targetFacets.merge(surfFacets[s]);  
    }
    LOG_WARNING << "No surface ID provided. Launching from all surfaces by default. WARNING - Will take a significant amount of time for larger geometry sets";
  }
  LOG_WARNING << "Total Number of Triangles rays launched from = "
              << numTargetFacets;
  moab::EntityHandle targetFacetsSet;
  moab::Range target_facets;
  DAG->moab_instance()->create_meshset(MESHSET_SET, targetFacetsSet);
  DAG->moab_instance()->add_entities(targetFacetsSet, targetFacets);
  target_facets.insert(targetFacetsSet);
  DAG->moab_instance()->write_file("target_facets.stl", NULL, NULL, target_facets);

  numFacets = targetFacets.size();
  return targetFacets;
}

int AegisClass::num_facets(){
  return numFacets;
}

void AegisClass::print_particle_stats(std::array<int, 4> particleStats){
  int particlesCounted = 0;
  for (const auto i:particleStats){
    particlesCounted += i;
  }

  LOG_WARNING << "Number of particles launched = " << particlesCounted<< std::endl;
  LOG_WARNING << "Number of particles depositing power from omp = " << particleStats[0] << std::endl;
  LOG_WARNING << "Number of shadowed particle intersections = " << particleStats[1] << std::endl;
  LOG_WARNING << "Number of particles lost from magnetic domain = " << particleStats[2] << std::endl;
  LOG_WARNING << "Number of particles terminated upon reaching max tracking length = " << particleStats[3] << std::endl; 
  LOG_WARNING << "Number of particles not accounted for = " << (num_facets() - particlesCounted) << std::endl;
}

void AegisClass::mpi_particle_stats(){
 
  for (int i=0; i<nprocs; ++i){
    if (rank == i){
      std::cout << std::endl << "process " << i << " has the following particle stats:" << std::endl;
      for (auto k: integrator->particle_stats())
      {
        std::cout << k << std::endl;
      } 
    }
  }
}




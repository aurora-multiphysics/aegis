#include "ParticleSimulation.h"
#include <mpi.h>
#include <memory>

ParticleSimulation::ParticleSimulation(std::string filename)
{
  JSONsettings = std::make_shared<InputJSON>(filename);

  read_params(JSONsettings);

  equilibrium.setup(JSONsettings);
  DAG = std::make_unique<moab::DagMC>();
  vtkInterface = std::make_unique<VtkInterface>(JSONsettings);

  init_geometry();
}

void ParticleSimulation::Execute(){

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  vtkInterface->init();  
  moab::Range targetSurfaceList = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetSurfaceList);

  std::vector<double> Bfield; 
  std::vector<double> polarPos(3);
  std::vector<double> newPt(3);  
  std::vector<double> triCoords(9);
  std::vector<double> triA(3), triB(3), triC(3);

  int totalNumberOfFacets = targetSurfaceList.size();
  int numberofFacets = totalNumberOfFacets / nprocs;
  int startFacet = rank * numberofFacets;
  int endFacet = startFacet + numberofFacets;
  nFacets = 0;
  
  for (int i=startFacet; i<endFacet; ++i){ // loop over all facets
    const auto facet = targetSurfaceList[i];
    ++nFacets;
    ParticleBase particle;
    std::vector<moab::EntityHandle> triNodes;
    DAG->moab_instance()->get_adjacencies(&facet, 1, 0, false, triNodes);
    DAG->moab_instance()->get_coords(&triNodes[0], triNodes.size(), triCoords.data());

    for (int j=0; j<3; j++){
      triA[j] = triCoords[j];
      triB[j] = triCoords[j+3];
      triC[j] = triCoords[j+6];
    }
    TriSource Tri(triA, triB, triC, facet); 

    if (particleLaunchPos == "fixed"){
      particle.set_pos(Tri.centroid());
    }
    else{
      particle.set_pos(Tri.random_pt());
    }
    integrator->set_launch_position(facet, particle.get_pos("cart"));
    
    particle.check_if_in_bfield(equilibrium); // if out of bounds skip to next triangle
    if (particle.outOfBounds){
      LOG_INFO << "Particle start is out of magnetic field bounds. Skipping to next triangle. Check correct eqdsk is being used for the given geometry";
      //vtkInterface->insert_zero_uStrGrid();
    qValues.push_back(0.0);
      continue;
    }
    vtkInterface->init_new_vtkPoints();
    // vtkInterface->insert_next_point_in_track(particle.launchPos);

    particle.set_dir(equilibrium);
    polarPos = CoordTransform::cart_to_polar(particle.launchPos, "forwards");
    BdotN = Tri.dot_product(particle.BfieldXYZ);

    psi = particle.get_psi(equilibrium); 
    psiOnSurface = psi;
    psiValues.push_back(psi);
    psid = psi + equilibrium.psibdry; 
    
    particle.align_dir_to_surf(BdotN);
    Q = equilibrium.omp_power_dep(psid, BdotN, "exp");
    
    // Start ray tracing
    DagMC::RayHistory history;
    ray_hit_on_launch(particle, history);
    trackLength = trackStepSize;

    particle.update_vectors(trackStepSize, equilibrium);
    vtkInterface->insert_next_point_in_track(particle.pos);

    // vtkInterface->insert_next_uStrGrid("Normal", Tri.unitNormal);
    // vtkInterface->insert_next_uStrGrid("B_field", particle.BfieldXYZ);
    // vtkInterface->insert_next_uStrGrid("Psi_Start", psi);
    // vtkInterface->insert_next_uStrGrid("B.n", BdotN);
    if (BdotN < 0){
      // vtkInterface->insert_next_uStrGrid("B.n_direction", -1.0);
    }
    else if (BdotN > 0){
      // vtkInterface->insert_next_uStrGrid("B.n_direction", 1.0);
    }

    loop_over_particle_track(facet, particle, history);

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

    vtkInterface->write_unstructuredGrid("out.vtk");
    vtkInterface->write_multiBlockData("particle_tracks.vtm");
    print_particle_stats(totalParticleStats);
  }

  mpi_particle_stats();

  }


void ParticleSimulation::read_params(const std::shared_ptr<InputJSON> &inputs){

  json aegisNamelist;
  if (inputs->data.contains("aegis_params"))
  {
    aegisNamelist = inputs->data["aegis_params"];
    dagmcInputFile = aegisNamelist["DAGMC"];
    trackStepSize = aegisNamelist["step_size"];
    maxTrackSteps = aegisNamelist["max_steps"];
    particleLaunchPos = aegisNamelist["launch_pos"];
    noMidplaneTermination = aegisNamelist["force_no_deposition"];

    if (aegisNamelist.contains("target_surfs"))
    {
      for (auto i:aegisNamelist["target_surfs"])
      {
        vectorOfTargetSurfs.push_back(i);
      }
    }
  }
  
  if (particleLaunchPos == "fixed")
  {
    LOG_WARNING << "Launching from triangle centroid/barycentre";
  }
  else
  {
    LOG_WARNING << "Launching from random positions in triangles";
  }

}

void ParticleSimulation::init_geometry(){
  // setup dagmc instance
  DAG->load_file(dagmcInputFile.c_str());
  DAG->init_OBBTree();
  DAG->setup_geometry(surfsList, volsList);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, facetsList);
  std::cout << "Number of Triangles in Geometry " << facetsList.size() << std::endl;
  volID = DAG->entity_by_index(3,1);

  // setup B Field data
  
  equilibrium.move();
  equilibrium.psiref_override();
  equilibrium.init_interp_splines();
  equilibrium.centre(1);

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

  equilibrium.psi_limiter(vertexList);
  equilibrium.write_bfield(plotBFieldRZ, plotBFieldXYZ);
}

void ParticleSimulation::loop_over_particle_track(const moab::EntityHandle &facet, ParticleBase &particle, DagMC::RayHistory &history)
{

  for (int step=0; step<maxTrackSteps; ++step)
  {
    trackLength += trackStepSize;
    DAG->ray_fire(volID, particle.pos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history, trackStepSize, rayOrientation);

    if (nextSurf != 0)
    {
      particle.update_vectors(nextSurfDist); // update position to surface intersection point
      terminate_particle(facet, history, terminationState::SHADOWED);
      break;
    }
    else
    {
      particle.update_vectors(trackStepSize); // update position by stepsize
    }
    
    vtkInterface->insert_next_point_in_track(particle.pos);
    particle.check_if_in_bfield(equilibrium);
    if (particle.outOfBounds){
      terminate_particle(facet, history, terminationState::LOST);
      break;
    }

    else
    {
      particle.set_dir(equilibrium);
      particle.align_dir_to_surf(BdotN);
    }
    particle.check_if_midplane_reached(equilibrium.get_midplane_params());
    
    if (particle.atMidplane != 0 && !noMidplaneTermination){ 
      terminate_particle(facet, history, terminationState::DEPOSITING);
      break;
    }

    if (step == (maxTrackSteps-1))
    {
      terminate_particle(facet, history, terminationState::MAXLENGTH);
      break;
    }

  } 

}

void ParticleSimulation::terminate_particle(const moab::EntityHandle &facet, DagMC::RayHistory &history, terminationState termination){
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

void ParticleSimulation::ray_hit_on_launch(ParticleBase &particle, DagMC::RayHistory &history){
  DAG->ray_fire(volID, particle.launchPos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history,trackStepSize,rayOrientation);
  if (nextSurf != 0) 
  {
    history.get_last_intersection(intersectedFacet);
    DAG->next_vol(nextSurf, volID, volID);
    LOG_INFO << "---- RAY HIT ON LAUNCH [" << facetCounter << "] ----";
  }
}

// Get triangles from the surface(s) of interest
moab::Range ParticleSimulation::select_target_surface(){ 

  moab::Range targetFacets; // range containing all of the triangles in the surface of interest
  // can specify particular surfaces of interest
  EntityHandle targetSurf; // surface of interest
  targetSurf = surfsList[0];

  std::unordered_map<EntityHandle, moab::Range> surfFacets;
  moab::Range surfFacetsKeys;
  int numTargetFacets = 0;


  if (!vectorOfTargetSurfs.empty())
  {
    for (auto &surfID:vectorOfTargetSurfs)
    {
      targetSurf = DAG->entity_by_id(2,surfID); // surface of interest
      DAG->moab_instance()->get_entities_by_type(targetSurf, MBTRI, surfFacets[targetSurf]);
      if (rank == 0)
      {
        std::cout << "Surface ID [" << surfID << "] " << "No. of elements " << surfFacets[targetSurf].size() << std::endl;
      }
      surfFacetsKeys.insert(targetSurf);
      numTargetFacets += surfFacets[targetSurf].size();
      targetFacets.merge(surfFacets[targetSurf]);  
    }

    LOG_WARNING << "Surface IDs provided. Launching from surfaces given by global IDs:";
  }
  else
  {
    targetFacets = facetsList;
    for (auto &surfEH:surfsList)
    {
      DAG->moab_instance()->get_entities_by_type(surfEH, MBTRI, surfFacets[surfEH]);
      surfFacetsKeys.insert(surfEH);
      if (rank == 0)
      {
        std::cout << "Surface ID [" << surfEH << "] " << "No. of elements " << surfFacets[surfEH].size() << std::endl;
      }
      numTargetFacets += surfFacets[surfEH].size();             
      targetFacets.merge(surfFacets[surfEH]); 
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

int ParticleSimulation::num_facets(){
  return numFacets;
}

void ParticleSimulation::print_particle_stats(std::array<int, 4> particleStats){
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

void ParticleSimulation::mpi_particle_stats(){
 
  for (int i=0; i<nprocs; ++i){
    if (rank == i){
      std::cout << std::endl << "process " << i << " has the following particle stats:" << std::endl;
      std::array localRankParticleStats = integrator->particle_stats();

      std::cout << "Depositing - " << localRankParticleStats[0] << std::endl;
      std::cout << "SHADOWED - " << localRankParticleStats[1] << std::endl;
      std::cout << "LOST - " << localRankParticleStats[2] << std::endl;
      std::cout << "MAX LENGTH - " << localRankParticleStats[3] << std::endl;
    }
  }
}
#include "ParticleSimulation.h"
#include <mpi.h>

// setup AEGIS simulation 
ParticleSimulation::ParticleSimulation(std::string filename)
{
  set_mpi_params();
  setup_profiling();

  JSONsettings = std::make_shared<InputJSON>(filename);

  read_params(JSONsettings);

  equilibrium.setup(JSONsettings);
  DAG = std::make_unique<moab::DagMC>();
  vtkInterface = std::make_unique<VtkInterface>(JSONsettings);

  init_geometry();
}

// perform AEGIS simulation solve with dynamic mpi load balancing
void ParticleSimulation::Execute_dynamic_mpi(){

  MPI_Status mpiStatus;

  vtkInterface->init();
  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  int totalFacets = targetFacets.size();
  int nFacetsPerProc = totalFacets / nprocs;

  const int noMoreWork = -1;

  std::vector<double> handlerQVals(num_facets());

  int workerStartIndex;
  MPI_Request request;

  std::ofstream fg("");


  if(rank == 0)  
  { // handler process...
    std::cout << "Dynamic task scheduling with " << (nprocs-1) << " processes, each handling "<< dynamicTaskSize << " facets" << std::endl;
    int facetsHandled = 0;
    for(int procID = 1; procID < nprocs; procID++)
    { // Send the initial indexes for each process statically
      MPI_Send(&facetsHandled, 1, MPI_INT, procID, procID, MPI_COMM_WORLD); 
      facetsHandled += dynamicTaskSize;
    }  


    std::vector<double> bufferQVals(dynamicTaskSize);

    int activeWorkers = nprocs - 1;

    do{
      int avaialbleProcess = 0;
      // call recieve from any source to get an avaialble process
      MPI_Recv(&avaialbleProcess, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &mpiStatus);
      if(facetsHandled < totalFacets)
      { // send that process the next index to start on in the total list
        MPI_Send(&facetsHandled, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD);  

        facetsHandled += dynamicTaskSize;
        MPI_Irecv(bufferQVals.data(), bufferQVals.size(), MPI_DOUBLE, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &request);
        MPI_Irecv(&workerStartIndex, 1, MPI_INT, MPI_ANY_SOURCE, 101, MPI_COMM_WORLD, &request);

        if (rank == 0) {std::cout << "Facets Handled = " << workerStartIndex  << std::endl;}

        // int workerEndIndex = workerStartIndex + dynamicTaskSize;
        // std::vector<double>::iterator itr = handlerQVals.begin() + workerStartIndex;
        // handlerQVals.insert(itr, bufferQVals.begin(), bufferQVals.end());
        
      }
      else
      { // send message to workers to tell them no more work available
        MPI_Send(&noMoreWork, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD);
        activeWorkers--;
      }

    } while(activeWorkers > 0); // while some workers are still active






    // for (int ctr=0; ctr<bufferQVals.size(); ++ctr)
    // {
    //   printf("Qvals on rank [%d] = %f[%d] \n", rank, bufferQVals[ctr], ctr+1);     
    // }
  

  }

  else
  { // worker processes...
    std::vector<double> bufferQVals;
    int startFacetIndex = 0;
    // recieve the initial set of indexes
    MPI_Recv(&startFacetIndex, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &mpiStatus);    
    // Work with the iterations startFacetIndex

    int start = startFacetIndex;
    int end = start + dynamicTaskSize;

    // processing the initial set of facets. Call loop over facets here
    bufferQVals = loop_over_facets(start, end);
    // send local (to worker) array back to handler process
    // MPI_Send(bufferQVals.data(), bufferQVals.size(), MPI_DOUBLE, 0, rank*100, MPI_COMM_WORLD);
    // MPI_Send(&start, 1, MPI_INT, 0, rank*101, MPI_COMM_WORLD);
    
    // for (int ctr=0; ctr<bufferQVals.size(); ++ctr)
    // {
    //   printf("Qvals on rank [%d] = %f[%d] \n", rank, bufferQVals[ctr], ctr+1);     
    // }

    // for (int i=0; i<dynamicTaskSize; ++i)
    // {
    //   initialQValsRankn << "Qvals on rank[" << rank << "] = " << bufferQVals[i] << "[" << i << "]" << std::endl; 
    // }

    do
    {
      MPI_Send(&rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Recv(&startFacetIndex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &mpiStatus);
      
      if (startFacetIndex + dynamicTaskSize > num_facets())
      {
        dynamicTaskSize = num_facets() - startFacetIndex;
      }

      start = startFacetIndex;
      end = start + dynamicTaskSize;

      if(startFacetIndex != noMoreWork)
      { // processing the next set of available work that has been dynamically allocated 
        bufferQVals = loop_over_facets(start, end);
        // send local (to worker) array back to handler process
        
        MPI_Igatherv(bufferQVals.data(), bufferQVals.size(), MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD, &request);
        // MPI_Send(bufferQVals.data(), bufferQVals.size(), MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
        MPI_Send(&start, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
 
        // for (int ctr=0; ctr<bufferQVals.size(); ++ctr)
        // {
        //   printf("Qvals on rank [%d] = %f[%d] \n", rank, bufferQVals[ctr], ctr+1); 
        // }
 
      }
    } while(startFacetIndex != noMoreWork); // while there is work available
  }

  // if (rank == 0)
  // {
  //   MPI_Recv(bufferQVals.data(), bufferQVals.size(), MPI_DOUBLE, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &mpiStatus);
  //   MPI_Recv(&workerStartIndex, 1, MPI_INT, MPI_ANY_SOURCE, 101, MPI_COMM_WORLD, &mpiStatus);

  //   int workerEndIndex = workerStartIndex + dynamicTaskSize;
  //   for (int i=workerStartIndex; i<workerEndIndex; ++i)
  //   {
  //   }
  // }

  // write out data and print final

  std::array<int, 5> particleStats = integrator->particle_stats(); 
  std::array<int, 5> totalParticleStats;


   if (rank != 0){
  //   MPI_Send(psiValues.data(), psiValues.size(), MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
  //   MPI_Send(qValues.data(), qValues.size(), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);

   }
   else 
   {
  //   allPsiValues.insert(allPsiValues.end(), psiValues.begin(), psiValues.end());
  //   allQValues.insert(allQValues.end(), qValues.begin(), qValues.end());
    totalParticleStats = integrator->particle_stats();
  }


  for (int i=1; i<nprocs; ++i){
    if (rank == 0){
  //     MPI_Recv(psiValues.data(), psiValues.size(), MPI_DOUBLE, i, 9, MPI_COMM_WORLD, &mpiStatus);
  //     allPsiValues.insert(allPsiValues.end(), psiValues.begin(), psiValues.end());

  //     MPI_Recv(qValues.data(), qValues.size(), MPI_DOUBLE, i, 10, MPI_COMM_WORLD, &mpiStatus);
  //     allQValues.insert(allQValues.end(), qValues.begin(), qValues.end());

      MPI_Recv(particleStats.data(), particleStats.size(), MPI_INT, i, 11, MPI_COMM_WORLD, &mpiStatus);
      for (int j=0; j<particleStats.size(); ++j){
        totalParticleStats[j] += particleStats[j]; 
      }
    }
  }


  //write_aegis_out()
  mpi_particle_stats();

  print_particle_stats(totalParticleStats);

  }

// Read parameters from aegis_settings.json config file
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
    dynamicTaskSize = aegisNamelist["dynamic_task_size"];


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
    log_string(LogLevel::WARNING, "Launching from triangle centroid/barycentre");
  }
  else
  {
    log_string(LogLevel::WARNING, "Launching from random positions in triangles");
  }

}

// initialise CAD geometry for AEGIS and magnetic field equilibrium for CAD Geometry
void ParticleSimulation::init_geometry(){
  // setup dagmc instance
  DAG->load_file(dagmcInputFile.c_str());
  DAG->init_OBBTree();
  DAG->setup_geometry(surfsList, volsList);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, facetsList);

  if (rank == 0) {std::cout << "Number of Triangles in Geometry " << facetsList.size() << std::endl;}
  implComplementVol = DAG->entity_by_index(3, volsList.size());
  
  if (!DAG->is_implicit_complement(implComplementVol))
    {
      std::cout << "Particle not in implicit complement. Check volumes." << std::endl;
    }

  // setup B Field data
  
  equilibrium.move();
  equilibrium.psiref_override();
  equilibrium.init_interp_splines();
  equilibrium.centre(1);

  std::vector<double> vertexCoordinates;
  DAG->moab_instance()->get_vertex_coordinates(vertexCoordinates);
  int numNodes = vertexCoordinates.size()/3;
  if (rank == 0) {std::cout << "NUMBER OF NODES = " <<  numNodes << std::endl;}

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

// loop over facets in target surfaces
std::vector<double> ParticleSimulation::loop_over_facets(int startFacet, int endFacet)
{ 

  int size = endFacet - startFacet;  
  std::vector<double> heatfluxVals;

  std::vector<double> Bfield; 
  std::vector<double> polarPos(3);
  std::vector<double> newPt(3);  
  std::vector<double> triCoords(9);
  std::vector<double> triA(3), triB(3), triC(3);

  for (int i=startFacet; i<endFacet; ++i)
  { // loop over all facets
    
    if (i >= targetFacets.size()) // handle padded values
    {
      heatfluxVals.push_back(-1);
      integrator->count_particle(0, terminationState::PADDED, 0.0);
      continue;
    }
    const auto facet = targetFacets[i];
    ParticleBase particle;
    std::vector<moab::EntityHandle> triNodes;
    DAG->moab_instance()->get_adjacencies(&facet, 1, 0, false, triNodes);
    DAG->moab_instance()->get_coords(&triNodes[0], triNodes.size(), triCoords.data());

    for (int j=0; j<3; j++)
    {
      triA[j] = triCoords[j];
      triB[j] = triCoords[j+3];
      triC[j] = triCoords[j+6];
    }
    TriSource Tri(triA, triB, triC, facet); 

    if (particleLaunchPos == "fixed")
    {
      particle.set_pos(Tri.centroid());
    }
    else
    {
      particle.set_pos(Tri.random_pt());
    }
    integrator->set_launch_position(facet, particle.get_pos("cart"));
    
    particle.check_if_in_bfield(equilibrium); // if out of bounds skip to next triangle
    if (particle.outOfBounds)
    {
      if (rank == 0) {LOG_INFO << "Particle start is out of magnetic field bounds. Skipping to next triangle. Check correct eqdsk is being used for the given geometry";}
      continue;
    }


    vtkInterface->init_new_vtkPoints();
    vtkInterface->insert_next_point_in_track(particle.launchPos);

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
    //ray_hit_on_launch(particle, history);
    trackLength = trackStepSize;

    particle.update_vectors(trackStepSize, equilibrium);

    vtkInterface->insert_next_point_in_track(particle.pos);


    terminationState state = loop_over_particle_track(facet, particle, history);

  if (state == terminationState::DEPOSITING) {heatfluxVals.push_back(Q);}
  else {heatfluxVals.push_back(0.0);}

  }


  return heatfluxVals;
}

// loop over single particle track from launch to termination
terminationState ParticleSimulation::loop_over_particle_track(const moab::EntityHandle &facet, ParticleBase &particle, DagMC::RayHistory &history)
{
  double euclidDistToNextSurf = 0.0;
  nextSurf = 0;

  for (int step=0; step<maxTrackSteps; ++step)
  {
    trackLength += trackStepSize;
    iterationCounter++;

    DAG->ray_fire(implComplementVol, particle.pos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history, trackStepSize, rayOrientation);
    numberOfRayFireCalls++; 


    if (nextSurf != 0)
    {
      particle.update_vectors(nextSurfDist); // update position to surface intersection point
      terminate_particle(facet, history, terminationState::SHADOWED);
      return terminationState::SHADOWED;
    }
    else
    {
      particle.update_vectors(trackStepSize); // update position by stepsize
    }
    
    vtkInterface->insert_next_point_in_track(particle.pos);
    particle.check_if_in_bfield(equilibrium);
    if (particle.outOfBounds){
      terminate_particle(facet, history, terminationState::LOST);
      return terminationState::LOST;
    }

    else
    {
      particle.set_dir(equilibrium);
      particle.align_dir_to_surf(BdotN);
    }
    particle.check_if_midplane_reached(equilibrium.get_midplane_params());
    
    if (particle.atMidplane != 0 && !noMidplaneTermination){ 
      terminate_particle(facet, history, terminationState::DEPOSITING);
      return terminationState::DEPOSITING;
    }

    if (step == (maxTrackSteps-1))
    {
      terminate_particle(facet, history, terminationState::MAXLENGTH);
      return terminationState::MAXLENGTH;
    }

  } 

  return terminationState::LOST; // default return lost particle

}

// terminate particle depending on 1 of 4 termination states:
//
// DEPOSITING - Particle reaches outer midplane and deposits power on facet. Heatflux = Q
// SHADOWED - Particle hits another piece of geometry. Heatflux = 0 
// LOST - Particle leaves magnetic field so trace stops. Heatflux = 0
// MAX LENGTH - Particle reaches maximum user set length before anything else. Heatflux = 0
void ParticleSimulation::terminate_particle(const moab::EntityHandle &facet, DagMC::RayHistory &history, terminationState termination){
  double heatflux;

  switch(termination){
    case terminationState::DEPOSITING:
      heatflux = Q;
      vtkInterface->write_particle_track(branchDepositingPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);
      if (rank == 0) {LOG_INFO << "Midplane reached. Depositing power";}
      break;

    case terminationState::SHADOWED:
      heatflux = 0.0;
      history.get_last_intersection(intersectedFacet);
      vtkInterface->write_particle_track(branchShadowedPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);
      if (rank == 0) {LOG_INFO << "Surface " << nextSurf << " hit after travelling " << trackLength << " units";}
      break;

    case terminationState::LOST:
      heatflux = 0.0;   
      vtkInterface->write_particle_track(branchLostPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);    
      if (rank == 0) {LOG_INFO << "TRACE STOPPED BECAUSE LEAVING MAGNETIC FIELD";}
      break;

    case terminationState::MAXLENGTH:
      heatflux = 0.0;
      vtkInterface->write_particle_track(branchMaxLengthPart, heatflux);
      integrator->count_particle(facet, termination, heatflux);        
      if (rank == 0) {LOG_INFO << "Fieldline trace reached maximum length before intersection";}
  }

  qValues.push_back(heatflux);
  psiQ_values.push_back(std::make_pair(psiOnSurface, Q));
}

void ParticleSimulation::ray_hit_on_launch(ParticleBase &particle, DagMC::RayHistory &history){
  DAG->ray_fire(implComplementVol, particle.launchPos.data(), particle.dir.data(), nextSurf, nextSurfDist, &history,trackStepSize,rayOrientation);
  if (nextSurf != 0) 
  {
    history.get_last_intersection(intersectedFacet);
    EntityHandle nextVolEH;
    DAG->next_vol(nextSurf, implComplementVol, nextVolEH);
    if (rank == 0) {LOG_INFO << "---- RAY HIT ON LAUNCH [" << facetCounter << "] ----";}
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
    std::stringstream surfaceIDsOut;
    surfaceIDsOut << "Surface IDs provided. Number of facets in each surface with global IDs: ";
    for (auto &surfID:vectorOfTargetSurfs)
    {
      targetSurf = DAG->entity_by_id(2,surfID); // surface of interest
      DAG->moab_instance()->get_entities_by_type(targetSurf, MBTRI, surfFacets[targetSurf]);
      if (rank == 0)
      {
        surfaceIDsOut << "[" << surfID << "] - " << surfFacets[targetSurf].size() << ", ";
      }
      surfFacetsKeys.insert(targetSurf);
      numTargetFacets += surfFacets[targetSurf].size();
      targetFacets.merge(surfFacets[targetSurf]);  
    }
    surfaceIDsOut << std::endl;
    log_string(LogLevel::WARNING, surfaceIDsOut.str());

  }
  else
  {
    targetFacets = facetsList;
    std::stringstream surfaceIDsOut;
    surfaceIDsOut << "No surface IDs provided. Launching from all surfaces by default. Facets per surface ID: ";

    for (auto &surfEH:surfsList)
    {
      DAG->moab_instance()->get_entities_by_type(surfEH, MBTRI, surfFacets[surfEH]);
      surfFacetsKeys.insert(surfEH);
      if (rank == 0)
      {
        surfaceIDsOut << "[" << surfEH << "] " << " - " << surfFacets[surfEH].size() << ",";
      }
      numTargetFacets += surfFacets[surfEH].size();             
      targetFacets.merge(surfFacets[surfEH]); 
    }
  log_string(LogLevel::WARNING, surfaceIDsOut.str());
  log_string(LogLevel::WARNING, "WARNING - Will take a significant amount of time for larger geometry sets");
  }
  if (rank == 0) {LOG_WARNING << "Total Number of Triangles rays launched from = " << numTargetFacets;}

  numFacets = targetFacets.size();
  return targetFacets;
}

// return total number of facets in target surface(s) of interest 
int ParticleSimulation::num_facets(){
  return numFacets;
}

// print particle stats for the entire run
void ParticleSimulation::print_particle_stats(std::array<int, 5> particleStats){
  int particlesCounted = 0;
  //MPI_Barrier(MPI_COMM_WORLD); // barrier to ensure this is always printed last
  for (const auto i:particleStats){
    particlesCounted += i;
  }
  particlesCounted -=  particleStats[4]; // remove padded particles
  if (rank == 0)
  {
    LOG_WARNING << "Number of particles launched = " << particlesCounted << std::endl;
    LOG_WARNING << "Number of particles depositing power from omp = " << particleStats[0] << std::endl;
    LOG_WARNING << "Number of shadowed particle intersections = " << particleStats[1] << std::endl;
    LOG_WARNING << "Number of particles lost from magnetic domain = " << particleStats[2] << std::endl;
    LOG_WARNING << "Number of particles terminated upon reaching max tracking length = " << particleStats[3] << std::endl; 
    LOG_WARNING << "Number of padded null particles = " << particleStats[4] << std::endl;
    LOG_WARNING << "Number of particles not accounted for = " << (num_facets() - particlesCounted) << std::endl;
  }
}

// print individual particle stats for each MPI rank
void ParticleSimulation::mpi_particle_stats(){
 
  for (int i=1; i<nprocs; ++i){
    if (rank == i){
      std::cout << std::endl << "process " << i << " has the following particle stats:" << std::endl;
      std::array localRankParticleStats = integrator->particle_stats();

      std::cout << "DEPOSITING - " << localRankParticleStats[0] << std::endl;
      std::cout << "SHADOWED - " << localRankParticleStats[1] << std::endl;
      std::cout << "LOST - " << localRankParticleStats[2] << std::endl;
      std::cout << "MAX LENGTH - " << localRankParticleStats[3] << std::endl;

      int totalParticlesHandled = localRankParticleStats[0] + localRankParticleStats[1]
                                + localRankParticleStats[2] + localRankParticleStats[3];

      std::cout << "TOTAL - " << totalParticlesHandled << std::endl;

    }
  }
}

// run AEGIS simulation on single core
void ParticleSimulation::Execute(){

  vtkInterface->init();  

  moab::ErrorCode rval;

  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  implicit_complement_testing();

  int start = 0;
  int end = targetFacets.size();

  std::vector<double> qvalues;
  qvalues = loop_over_facets(start, end);
  if (qvalues.empty())
  {
    log_string(LogLevel::ERROR, "Error - loop over facets returned no heatfluxes, please check logfile. Exiting...");
    std::exit(EXIT_FAILURE);
  }
  attach_mesh_attribute("Heatflux", targetFacets, qvalues);
  write_out_mesh(meshWriteOptions::BOTH); 


  // write out data and print final

  std::array<int, 5> particleStats = integrator->particle_stats(); 

    vtkInterface->write_multiBlockData("particle_tracks.vtm");
    print_particle_stats(particleStats);

    std::cout << "Number of Ray fire calls = " << numberOfRayFireCalls << std::endl;
    std::cout << "Number of Closest Location calls = " << numberOfClosestLocCalls << std::endl;
    std::cout << "Number of iterations performed = " << iterationCounter << std::endl;
  }

void ParticleSimulation::Execute_padded_mpi(){

  vtkInterface->init();  
  MPI_Status mpiStatus;

  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  int totalFacets = num_facets();
  int remainder = totalFacets % nprocs;
  int nPadded = 0;
  if (remainder != 0)
  {
    int paddedTotalFacets = num_facets() + (nprocs-remainder);
    nPadded = paddedTotalFacets - totalFacets;
    totalFacets = paddedTotalFacets;
  }
  int nFacetsPerProc = totalFacets / nprocs;
  int startFacet = rank * nFacetsPerProc;
  int endFacet = startFacet + nFacetsPerProc;
  int root_rank = 0;


  std::vector<double> qvalues; // qvalues buffer local to each processor
  std::vector<double> rootQvalues(totalFacets); // total qvalues buffer on root process for IO  

  qvalues = loop_over_facets(startFacet, endFacet); // perform main loop

  std::array<int, 5> particleStats = integrator->particle_stats(); 
  std::array<int, 5> totalParticleStats;
  std::vector<double> nonRootQvals(nFacetsPerProc);
  if (rank != 0)
  { 
    MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);
  }
  else
  {
    totalParticleStats = integrator->particle_stats();
    MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, rootQvalues.data(), qvalues.size(), MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    if (rootQvalues.size() > num_facets())
    {
      rootQvalues.resize(rootQvalues.size() - nPadded);
    }

    for(int i=1; i<nprocs; ++i)
    {
      MPI_Recv(particleStats.data(), particleStats.size(), MPI_INT, i, 11, MPI_COMM_WORLD, &mpiStatus);
      for (int j=0; j<particleStats.size(); ++j){
        totalParticleStats[j] += particleStats[j]; 
      }
    }
    attach_mesh_attribute("Heatflux", targetFacets, rootQvalues);
    write_out_mesh(meshWriteOptions::BOTH, targetFacets); 
  }


  mpi_particle_stats();
  print_particle_stats(totalParticleStats);


}

void ParticleSimulation::Execute_mpi(){
  startTime = MPI_Wtime();
  vtkInterface->init();  
  MPI_Status mpiStatus;

  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  int totalFacets = num_facets();
  int nFacetsPerProc = totalFacets / nprocs;
  int remainder = totalFacets % nprocs;

  int root_rank = 0;


  std::vector<double> qvalues; // qvalues buffer local to each processor
  std::vector<double> rootQvalues(totalFacets); // total qvalues buffer on root process for IO  

  std::vector<int> recieveCounts(nprocs, nFacetsPerProc); 
  std::vector<int>::iterator recvItr = recieveCounts.begin();

  int remaindersHandled = 0;  
  while (remaindersHandled < remainder) // distribute remainders across processes
  {
    *recvItr += 1;
    std::advance(recvItr, 1);
    remaindersHandled++;
    if (recvItr == recieveCounts.end()) { recvItr = recieveCounts.begin(); }
  }

  std::vector<int> displacements(nprocs, 0);
  int displ = 0;
  
  for (int i=0; i<nprocs; ++i) // calculate displacements
  {
    displacements[i] = displ;
    displ += recieveCounts[i];
  }

  int startFacet = displacements[rank];
  int endFacet = startFacet + recieveCounts[rank];

  startTime = MPI_Wtime();
  qvalues = loop_over_facets(startFacet, endFacet); // perform main loop
  endTime = MPI_Wtime();

  std::array<int, 5> particleStats = integrator->particle_stats(); 
  std::array<int, 5> totalParticleStats;
  if (rank != 0)
  { 
    MPI_Gatherv(qvalues.data(), qvalues.size(), MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    // MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);
  }
  else
  {
    totalParticleStats = integrator->particle_stats();
    MPI_Gatherv(qvalues.data(), qvalues.size(), MPI_DOUBLE, rootQvalues.data(), recieveCounts.data(), displacements.data(), MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    // MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, rootQvalues.data(), qvalues.size(), MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    // for (auto i:rootQvalues) { temp << i << std::endl; }

    for(int i=1; i<nprocs; ++i)
    {
      MPI_Recv(particleStats.data(), particleStats.size(), MPI_INT, i, 11, MPI_COMM_WORLD, &mpiStatus);
      for (int j=0; j<particleStats.size(); ++j){
        totalParticleStats[j] += particleStats[j]; 
      }
    }
    attach_mesh_attribute("Heatflux", targetFacets, rootQvalues);
    write_out_mesh(meshWriteOptions::BOTH, targetFacets); 
  }


  mpi_particle_stats();
  print_particle_stats(totalParticleStats);

  endTime = MPI_Wtime();
  out_mpi_wtime("Execute_mpi()", endTime - startTime);


}

// get surfaces attributed to the aegis_target group set in cubit
void ParticleSimulation::implicit_complement_testing()
{

  std::vector<EntityHandle> children;
  DAG->moab_instance()->get_child_meshsets(implComplementVol, children, 1);
  std::vector<EntityHandle> vertices;

  for (const auto &i:children)
  {
    std::vector<EntityHandle> temp;
    DAG->moab_instance()->get_entities_by_type(i, MBVERTEX, vertices, false); 
    vertices.insert(vertices.begin(), temp.begin(), temp.end());
  }
  std::vector<double> vertexCoords(vertices.size()*3);

  DAG->moab_instance()->get_coords(&vertices[0], vertices.size(), vertexCoords.data());

  // DAG->moab_instance()->list_entities(children);

  std::ofstream implcitComplCoords("implcit_complement.txt");
  for (int i=0; i<vertexCoords.size(); i+=3)
  {
    implcitComplCoords << vertexCoords[i] << " " << vertexCoords[i+1] << " " << vertexCoords[i+2] << std::endl;
  }

  std::cout << "Number of Nodes in implicit complement = " << vertexCoords.size() << std::endl;


  return;

}

  
// create a moab tag and attach that data to each element in the moab::Range provided
void ParticleSimulation::attach_mesh_attribute(const std::string &tagName, moab::Range &entities, std::vector<double> &dataToAttach)
{
  moab::Tag tag;
  double tagValue;

  DAG->moab_instance()->tag_get_handle(tagName.c_str(), 1, MB_TYPE_DOUBLE, tag, MB_TAG_CREAT | MB_TAG_DENSE, &tagValue);
  DAG->moab_instance()->tag_set_data(tag, entities, &dataToAttach[0]);

}

// write out the mesh and target mesh with attributed data
// meshWriteOptions::FULL -- Write out the entire mesh 
// meshWriteOptions::TARGET -- Only write out target mesh with heatfluxes (useful if full mesh particularly large)
// meshWriteOptions::PARTIAL -- Write out specified entities
// default -- Full mesh is written
void ParticleSimulation::write_out_mesh(meshWriteOptions option, moab::Range rangeofEntities)
{
  DAG->remove_graveyard();
  EntityHandle targetMeshset;
  DAG->moab_instance()->create_meshset(MESHSET_SET, targetMeshset);
  DAG->moab_instance()->add_entities(targetMeshset, targetFacets);

  switch (option)
  {
  case meshWriteOptions::FULL:
    DAG->write_mesh("aegis_full.vtk", 1);
    break;
  
  case meshWriteOptions::TARGET:
    DAG->moab_instance()->write_mesh("aegis_target.vtk", &targetMeshset, 1);
    break;

  case meshWriteOptions::BOTH:
    DAG->moab_instance()->write_mesh("aegis_target.vtk", &targetMeshset, 1);
    DAG->write_mesh("aegis_full.vtk", 1);
    break;

  case meshWriteOptions::PARTIAL:
    if (!rangeofEntities.empty()){
      EntityHandle meshset;
      DAG->moab_instance()->create_meshset(MESHSET_SET, meshset);
      DAG->moab_instance()->add_entities(meshset, rangeofEntities);
      DAG->moab_instance()->write_mesh("aegis_partial.vtk", &meshset, 1);
    }
    else{
      log_string(LogLevel::ERROR, "No meshsets provided for partial mesh write out. Defaulting to full mesh");
      DAG->write_mesh("aegis_full.vtk", 1);
    }
    break;
  
  default: // default to full mesh and target
    DAG->moab_instance()->write_mesh("aegis_out_target.vtk", &targetMeshset, 1);
    DAG->write_mesh("aegis_full.vtk", 1);
    break;
  }

}
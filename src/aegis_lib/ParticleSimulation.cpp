#include "ParticleSimulation.h"
#include <mpi.h>

// setup AEGIS simulation
ParticleSimulation::ParticleSimulation(std::shared_ptr<JsonHandler> configFile,
                                       std::shared_ptr<EquilData> equilibirum,
                                       std::shared_ptr<SurfaceIntegrator> integrator)
{
  set_mpi_params();
  read_params(configFile);
  equilibrium = equilibirum;
  DAG = std::make_unique<moab::DagMC>();
  vtkInterface = std::make_unique<VtkInterface>(configFile);
  init_geometry();
  _integrator = integrator;
  _integrator->set_facets(targetFacets);
}

// Call different execute functions depending on which execute type selected
void
ParticleSimulation::Execute()
{
  string_to_lowercase(exeType);

  if (exeType == "serial" || nprocs == 1)
  {
    log_string(LogLevel::WARNING, "Executing in serial mode...");
    if (nprocs > 1)
    {
      log_string(LogLevel::ERROR, "Aegis was invoked with mpirun but config suggests serial run. "
                                  "Defaulting to MPI with no load balancing...");
      Execute_mpi();
    }
    else
    {
      Execute_serial();
    }
  }

  else if (exeType == "mpi")
  {
    log_string(LogLevel::WARNING, "Executing MPI with no load balancing...");
    Execute_mpi();
  }

  else if (exeType == "mpi_dynamic" || exeType == "dynamic")
  {
    log_string(LogLevel::WARNING, "Executing dynamic MPI load balancing...");
    Execute_dynamic_mpi();
  }

  else
  {
    log_string(LogLevel::FATAL, "No suitable execution type ('serial', 'mpi', "
                                "'mpi_dynamic') provided exiting...");
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::exit(1);
  }
}

// perform AEGIS simulation solve with dynamic mpi load balancing
void
ParticleSimulation::Execute_dynamic_mpi()
{

  MPI_Status mpiStatus;

  const int noMoreWork = -1;

  std::vector<double> particleHeatfluxes(num_particles_launched());
  std::vector<double> facetHeatfluxes(target_num_facets());

  double mainParticleLoopStart = MPI_Wtime();

  if (rank == 0)
  { // handler process...
    handler(particleHeatfluxes);
    // sum particles over triangles

    auto particleIter = particleHeatfluxes.begin();
    for (unsigned int i = 0; i < facetHeatfluxes.size(); ++i)
    {
      facetHeatfluxes[i] = std::reduce(particleIter, (particleIter + nParticlesPerFacet));
      particleIter += nParticlesPerFacet; // step iterator by n particles
    }

    // this could be a method in the worker-handler class
    if (debugDynamicBatching)
    {
      std::ofstream handlerOutputFile;
      std::stringstream fileName;
      fileName << "handler_rank_0.txt";
      handlerOutputFile.open(fileName.str());
      int depositCounter = 0;
      for (const auto i : particleHeatfluxes)
      {
        handlerOutputFile << i << "\n";
        if (i > 0)
        {
          depositCounter += 1;
        }
      }
      handlerOutputFile << std::endl;
      handlerOutputFile << "Number of depositing particles = " << depositCounter << std::endl;
    }
  }

  else
  { // worker processes...
    // abstract this out to its own worker-handler classes?
    worker();
  }

  std::array<int, 4> particleStats = _integrator->particle_stats();
  std::array<int, 4> totalParticleStats;

  if (rank != 0)
  {
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);
  }
  else
  {
    totalParticleStats = _integrator->particle_stats();
  }

  for (int i = 1; i < nprocs; ++i)
  {
    if (rank == 0)
    {
      MPI_Recv(particleStats.data(), particleStats.size(), MPI_INT, i, 11, MPI_COMM_WORLD,
               &mpiStatus);
      for (int j = 0; j < particleStats.size(); ++j)
      {
        totalParticleStats[j] += particleStats[j];
      }
    }
  }

  // write out data and print final

  if (rank == 0)
  {
    mainParticleLoopTime = MPI_Wtime() - mainParticleLoopStart;

    double aegisMeshWriteStart = MPI_Wtime();
    if (writeParticleLaunchPos == true)
    {
      write_particle_launch_positions(particleHeatfluxes);
    }
    attach_mesh_attribute("Heatflux", targetFacets, facetHeatfluxes);
    write_out_mesh(meshWriteOptions::BOTH, targetFacets);
    aegisMeshWriteTime = MPI_Wtime() - aegisMeshWriteStart;
  }

  print_particle_stats(totalParticleStats);
}

void
ParticleSimulation::handler(std::vector<double> & handlerQVals)
{
  MPI_Status status;
  MPI_Request request;
  int mpiDataTag = 100;
  int mpiIndexTag = 101;
  int workerStartIndex = 0;
  int noMoreWork = -1;

  std::vector<double> workerQVals(dynamicBatchSize);
  std::cout << "Dynamic task scheduling with " << (nprocs - 1) << " processes, each handling "
            << dynamicBatchSize << " particles" << std::endl;

  int particlesHandled = 0;
  for (int procID = 1; procID < nprocs; procID++)
  { // Send the initial indexes for each process statically
    MPI_Send(&particlesHandled, 1, MPI_INT, procID, procID, MPI_COMM_WORLD);
    particlesHandled += dynamicBatchSize;
  }

  int inactiveWorkers = 0;
  do
  {
    int avaialbleProcess = 0;
    MPI_Recv(&avaialbleProcess, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
             &status); // get avaialble process
    if (particlesHandled < handlerQVals.size())
    { // send that process the next index to start on
      // in the total list
      MPI_Isend(&particlesHandled, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD, &request);

      MPI_Irecv(workerQVals.data(), workerQVals.size(), MPI_DOUBLE, MPI_ANY_SOURCE, mpiDataTag,
                MPI_COMM_WORLD, &request);
      MPI_Irecv(&workerStartIndex, 1, MPI_INT, MPI_ANY_SOURCE, mpiIndexTag, MPI_COMM_WORLD,
                &request);
      // printf("Time taken to recieve next data = %f \n", MPI_Wtime());

      particlesHandled += dynamicBatchSize;

      MPI_Wait(&request, MPI_STATUS_IGNORE);
      double * ptr = handlerQVals.data() + workerStartIndex;
      memcpy(ptr, workerQVals.data(), sizeof(double) * workerQVals.size());
    }
    else
    { // send message to workers to tell them no more work available
      MPI_Send(&noMoreWork, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD);
      inactiveWorkers++;

      // recieve final work arrays
      int finalArraySize;
      MPI_Probe(MPI_ANY_SOURCE, mpiDataTag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &finalArraySize);
      std::vector<double> finalWorkArray(finalArraySize);

      MPI_Recv(finalWorkArray.data(), finalWorkArray.size(), MPI_DOUBLE, status.MPI_SOURCE,
               mpiDataTag, MPI_COMM_WORLD, &status);
      MPI_Recv(&workerStartIndex, 1, MPI_INT, status.MPI_SOURCE, mpiIndexTag, MPI_COMM_WORLD,
               &status);

      double * ptr = handlerQVals.data() + workerStartIndex;
      memcpy(ptr, finalWorkArray.data(), sizeof(double) * finalWorkArray.size());
    }

  } while (inactiveWorkers < nprocs - 1); // while some workers are still active
}

void
ParticleSimulation::worker()
{
  int noMoreWork = -1;
  MPI_Status status;
  MPI_Request request;
  int mpiDataTag = 100;
  int mpiIndexTag = 101;
  int counterDeposit = 0;
  int index = 0;

  std::ofstream workerOutputFile;
  if (debugDynamicBatching)
  {
    std::stringstream fileName;
    fileName << "rank" << rank << ".txt";
    workerOutputFile.open(fileName.str());
  }

  std::vector<double> workerQVals;
  int particleStartIndex = 0;
  MPI_Recv(&particleStartIndex, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);

  unsigned int start = particleStartIndex;
  unsigned int end = start + dynamicBatchSize;

  workerQVals = loop_over_particles(start, end); // process initial facets
  index = start;

  if (debugDynamicBatching)
  {
    workerOutputFile << "Loop over [" << start << ":" << end << "] particles: \n";
    index = start;
    for (auto i : workerQVals)
    {
      workerOutputFile << "[" << index + 1 << "] " << i << "\n";
      if (i > 0)
      {
        counterDeposit += 1;
      }
      index += 1;
    }
    workerOutputFile << std::endl;
  }

  MPI_Send(workerQVals.data(), workerQVals.size(), MPI_DOUBLE, 0, mpiDataTag, MPI_COMM_WORLD);
  MPI_Send(&start, 1, MPI_INT, 0, mpiIndexTag, MPI_COMM_WORLD);
  do
  {
    MPI_Send(&rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Recv(&particleStartIndex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

    if (particleStartIndex + dynamicBatchSize > num_particles_launched())
    {
      dynamicBatchSize = num_particles_launched() - particleStartIndex;
    }

    start = particleStartIndex;
    end = start + dynamicBatchSize;

    if (particleStartIndex != noMoreWork)
    { // processing the next set of available work
      // that has been dynamically allocated
      workerQVals = loop_over_particles(start, end);
      // send local (to worker) array back to handler process
      MPI_Send(workerQVals.data(), workerQVals.size(), MPI_DOUBLE, 0, mpiDataTag, MPI_COMM_WORLD);
      MPI_Send(&start, 1, MPI_INT, 0, mpiIndexTag, MPI_COMM_WORLD);
      index = start;

      if (debugDynamicBatching)
      {
        workerOutputFile << "Loop over [" << start << ":" << end << "] particles: \n";
        for (auto i : workerQVals)
        {
          workerOutputFile << "[" << index + 1 << "] " << i << "\n";
          if (i > 0)
          {
            counterDeposit += 1;
          }
          index += 1;
        }
        workerOutputFile << std::endl;
      }
    }
  } while (particleStartIndex != noMoreWork); // while there is work available

  workerOutputFile << "Particle Stats:" << std::endl;
  workerOutputFile << "Depositing particles = " << counterDeposit << std::endl;
}

// Read parameters from aegis_settings.json config file
void
ParticleSimulation::read_params(const std::shared_ptr<JsonHandler> & configFile)
{

  if (configFile->data().contains("aegis_params"))
  {
    nlohmann::json aegisParamsData = configFile->data()["aegis_params"];
    JsonHandler aegisParams(aegisParamsData);

    dagmcInputFile = aegisParams.get_required<std::string>("DAGMC");

    trackStepSize = aegisParams.get_optional<double>("step_size").value_or(trackStepSize);

    maxTrackSteps = aegisParams.get_optional<int>("max_steps").value_or(maxTrackSteps);

    particleLaunchPos =
        aegisParams.get_optional<std::string>("launch_pos").value_or(particleLaunchPos);

    if (aegisParamsData.contains("monte_carlo_params"))
    {
      nlohmann::json monteCarloParamsData =
          configFile->data()["aegis_params"]["monte_carlo_params"];
      JsonHandler monteCarloParams(monteCarloParamsData);

      nParticlesPerFacet = monteCarloParams.get_optional<int>("number_of_particles_per_facet")
                               .value_or(nParticlesPerFacet);

      writeParticleLaunchPos =
          monteCarloParams.get_optional<bool>("write_particle_launch_positions")
              .value_or(writeParticleLaunchPos);
    }

    if (particleLaunchPos == "fixed")
    {
      nParticlesPerFacet = 1;
    }

    noMidplaneTermination =
        aegisParams.get_optional<bool>("force_no_deposition").value_or(noMidplaneTermination);

    coordinateConfig =
        aegisParams.get_optional<std::string>("coordinate_system").value_or(coordinateConfig);

    exeType = aegisParams.get_optional<std::string>("execution_type").value_or(exeType);

    if (aegisParamsData.contains("dynamic_batch_params"))
    {
      nlohmann::json dynamicBatchParamsData =
          configFile->data()["aegis_params"]["dynamic_batch_params"];
      JsonHandler dynamicBatchParams(dynamicBatchParamsData);

      dynamicBatchSize =
          dynamicBatchParams.get_optional<int>("batch_size").value_or(dynamicBatchSize);

      workerProfiling =
          dynamicBatchParams.get_optional<bool>("worker_profiling").value_or(workerProfiling);

      debugDynamicBatching =
          dynamicBatchParams.get_optional<bool>("debug").value_or(debugDynamicBatching);
    }

    if (aegisParamsData.contains("target_surfs"))
    {
      for (auto i : aegisParamsData["target_surfs"])
      {
        vectorOfTargetSurfs.push_back(i);
      }
    }
    
    auto _str = aegisParams.get_optional<std::string>("mesh_units");
    if (_str) { // Check if the optional has a value and then set to enum
      if (*_str == "m") meshDimension = MeshUnits::M;
      else if (*_str == "cm") meshDimension = MeshUnits::CM;
      else if (*_str == "mm") meshDimension = MeshUnits::MM;
    }

    _str = aegisParams.get_optional<std::string>("write_mesh");
    if (_str) { // Check if the optional has a value and then set to enum
      if (*_str == "target") writeMesh = meshWriteOptions::TARGET;
      else if (*_str == "full") writeMesh = meshWriteOptions::FULL;
      else if (*_str == "partial") writeMesh = meshWriteOptions::PARTIAL;
    }

  }

  std::string launchType = "Launching particles from triangles at ";
  launchType += (particleLaunchPos == "fixed") ? "centroid position" : "random position";
  log_string(LogLevel::WARNING, launchType);
}

// select coordinate system from json input
void
ParticleSimulation::select_coordinate_system()
{
  string_to_lowercase(coordinateConfig);
  if (coordinateConfig == "cart" || coordinateConfig == "cartesian" || coordinateConfig == "xyz")
  {
    coordSys = coordinateSystem::CARTESIAN;
    log_string(LogLevel::WARNING, "Tracking in CARTESIAN coordinates...");
  }
  else if (coordinateConfig == "pol" || coordinateConfig == "polar" || coordinateConfig == "rz")
  {
    coordSys = coordinateSystem::POLAR;
    log_string(LogLevel::WARNING, "Tracking in POLAR coordinates...");
  }
  else if (coordinateConfig == "flux" || coordinateConfig == "psi")
  {
    coordSys = coordinateSystem::FLUX;
    log_string(LogLevel::WARNING, "Tracking in FLUX coordinates...");
  }
  else
  {
    coordSys = coordinateSystem::CARTESIAN;
    log_string(LogLevel::WARNING, "Invalid coordinate system, defaulting to cartesian...");
  }
}


void ParticleSimulation::scale_mesh(MeshUnits from)
{
  float scalingFactor = 1.0f;
  switch (from)
  {
  case MeshUnits::CM:
    scalingFactor = 0.01f;
    break;

  case MeshUnits::MM:
    scalingFactor = 0.001f;
    break;
  }

  for (auto &i:nodeCoords)
  {
    i = i*scalingFactor;
  }
  DAG->moab_instance()->set_coords(&nodesList[0], nodesList.size(), nodeCoords.data());
}

// initialise CAD geometry for AEGIS and magnetic field equilibrium for CAD
// Geometry
void
ParticleSimulation::init_geometry()
{
  double dagmcMeshReadStart = MPI_Wtime();

  // setup dagmc instance
  DAG->load_file(dagmcInputFile.c_str());
  DAG->init_OBBTree();
  DAG->setup_geometry(surfsList, volsList);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, facetsList);

  DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, nodesList, true);
  nodeCoords.resize(nodesList.size() * 3);
  DAG->moab_instance()->get_coords(&nodesList[0], nodesList.size(), nodeCoords.data());

  if (meshDimension != MeshUnits::UNSET) scale_mesh (meshDimension); // scale mesh to Metres 

  totalNumberOfFacets = facetsList.size();

  implComplementVol = DAG->entity_by_index(3, volsList.size());

  if (!DAG->is_implicit_complement(implComplementVol))
  {
    std::cout << "Particle not in implicit complement. Check volumes." << std::endl;
  }
  dagmcMeshReadTime = MPI_Wtime() - dagmcMeshReadStart;

  // setup B Field data
  double prepSurfacesStart = MPI_Wtime();

  std::vector<double> vertexCoordinates;
  DAG->moab_instance()->get_vertex_coordinates(vertexCoordinates);
  int numNodes = vertexCoordinates.size() / 3;
  if (rank == 0)
  {
    std::cout << "NUMBER OF NODES = " << numNodes << std::endl;
  }

  std::vector<std::vector<double>> vertexList(numNodes, std::vector<double>(3));

  for (int i = 0; i < numNodes; i++)
  {
    vertexList[i][0] = vertexCoordinates[i];
    vertexList[i][1] = vertexCoordinates[i + numNodes];
    vertexList[i][2] = vertexCoordinates[i + numNodes * 2];
  }

  equilibrium->psi_limiter(vertexList);
  select_coordinate_system();

  // if flux selected default to cartesian
  if (coordSys == coordinateSystem::FLUX)
  {
    coordSys = coordinateSystem::CARTESIAN;
    log_string(LogLevel::WARNING, "Flux coords tracking not currently "
                                  "implemented. Defaulting to CARTESIAN...");
  }

  // transform coordinate system of mesh
  if (coordSys != coordinateSystem::CARTESIAN)
  {
    mesh_coord_transform(coordSys);
  }

  targetFacets = select_target_surface();
  prepSurfacesTime = MPI_Wtime() - prepSurfacesStart;

  double setupArrayOfParticlesTimeStart = MPI_Wtime();
  setup_sources();
  setupArrayOfParticlesTime = MPI_Wtime() - setupArrayOfParticlesTimeStart;
}

// loop over facets in target surfaces
std::vector<double>
ParticleSimulation::loop_over_particles(int startIndex, int endIndex)
{
  double startTime = MPI_Wtime();
  std::vector<double> heatfluxVals;
  heatfluxVals.reserve(endIndex);

  // main loop over all facets
  for (int i = startIndex; i < endIndex; ++i)
  {
    auto particle = listOfParticles[i];
    _integrator->set_launch_position(particle.parent_entity_handle(), particle.get_pos());
    terminationState particleState = cartesian_particle_track(particle);
    double heatflux = (particleState == terminationState::DEPOSITING) ? particle.heatflux() : 0.0;
    heatfluxVals.emplace_back(heatflux);
  }

  double endTime = MPI_Wtime();
  if (workerProfiling)
  {
    printf("Loop over facets [%d:%d] from rank[%d] time taken = %fs \n", startIndex, endIndex, rank,
           endTime - startTime);
    fflush(stdout);
  }
  return heatfluxVals;
}

void
ParticleSimulation::setup_sources()
{
  std::vector<double> triangleCoords(9);
  std::vector<double> trianglePtsA(3), trianglePtsB(3), trianglePtsC(3);

  // Preallocate exactly num_particles_launched() elements
  listOfParticles.reserve(num_particles_launched());

  // Index-based modification
  size_t particleIndex = 0;
  for (const auto &facetEH : targetFacets) {
    std::vector<moab::EntityHandle> triangleNodes;
    DAG->moab_instance()->get_adjacencies(&facetEH, 1, 0, false, triangleNodes);
    DAG->moab_instance()->get_coords(&triangleNodes[0], triangleNodes.size(),
                                     triangleCoords.data());

    for (int j = 0; j < 3; j++) {
      trianglePtsA[j] = triangleCoords[j];
      trianglePtsB[j] = triangleCoords[j + 3];
      trianglePtsC[j] = triangleCoords[j + 6];
    }

    TriangleSource triangle(trianglePtsA, trianglePtsB, trianglePtsC, facetEH, particleLaunchPos);
    if (!equilibrium->check_if_in_bfield(triangle.launch_pos())) {
      log_string(LogLevel::INFO,
                 "Triangle start outside of magnetic field. Skipping to next triangle...");
      continue;
    }
    triangle.set_heatflux_params(equilibrium, "exp");
    double psid = triangle.get_psi() + equilibrium->psibdry;
    double heatflux = equilibrium->omp_power_dep(psid, triangle.BdotN(), "exp");
    triangle.update_heatflux(heatflux);
    double heatfluxPerParticle = heatflux / nParticlesPerFacet;

    for (int i = 0; i < nParticlesPerFacet; ++i) {
      auto launchPos = (particleLaunchPos == "fixed") ? triangle.centroid() : triangle.random_pt();
      ParticleBase particle(coordinateSystem::CARTESIAN, launchPos, heatfluxPerParticle, facetEH);
      particle.align_dir_to_surf(triangle.BdotN());
      listOfParticles.emplace_back(particle);
    }
  }
}

// loop over single particle track from launch to termination
terminationState
ParticleSimulation::cartesian_particle_track(ParticleBase & particle)
{
  double euclidDistToNextSurf = 0.0;
  nextSurf = 0;
  particle.set_dir(equilibrium);

  vtkInterface->init_new_vtkPoints();
  vtkInterface->insert_next_point_in_track(particle.get_xyz_pos());
  moab::DagMC::RayHistory history;

  for (int step = 0; step < maxTrackSteps; ++step)
  {
    trackLength += trackStepSize;
    iterationCounter++;

    DAG->ray_fire(implComplementVol, particle.pos.data(), particle.dir.data(), nextSurf,
                  nextSurfDist, &history, trackStepSize, rayOrientation);
    numberOfRayFireCalls++;

    if (nextSurf != 0)
    {
      particle.update_vectors(nextSurfDist); // update position to surface intersection point
      terminate_particle_shadow(particle);
      return terminationState::SHADOWED;
    }
    else
    {
      particle.update_vectors(trackStepSize); // update position by stepsize
    }

    vtkInterface->insert_next_point_in_track(particle.get_xyz_pos());

    if (!equilibrium->check_if_in_bfield(particle.pos))
    {
      terminate_particle_lost(particle);
      return terminationState::LOST;
    }

    particle.set_dir(equilibrium);
    particle.check_if_midplane_crossed(equilibrium->get_midplane_params());

    if (particle.atMidplane != 0 && !noMidplaneTermination)
    {
      terminate_particle_depositing(particle);
      return terminationState::DEPOSITING;
    }
  }

  particle.set_facet_history(history);
  // max length of particle track reached
  terminate_particle_maxlength(particle);
  return terminationState::MAXLENGTH;
}

// terminate particle after reaching midplane
void
ParticleSimulation::terminate_particle_depositing(ParticleBase & particle)
{
  std::stringstream terminationLogString;
  double heatflux = particle.heatflux();
  vtkInterface->write_particle_track(branchDepositingPart, heatflux);
  _integrator->count_particle(particle.parent_entity_handle(), terminationState::DEPOSITING,
                              heatflux);
  terminationLogString << "Midplane reached. Depositing power after travelling " << trackLength
                       << " units";
  log_string(LogLevel::INFO, terminationLogString.str());
}

// terminate particle after hitting shadowing geometry
void
ParticleSimulation::terminate_particle_shadow(ParticleBase & particle)
{
  std::stringstream terminationLogString;
  double heatflux = 0.0;
  vtkInterface->write_particle_track(branchShadowedPart, heatflux);
  _integrator->count_particle(particle.parent_entity_handle(), terminationState::SHADOWED,
                              heatflux);
  terminationLogString << "Surface " << nextSurf << " hit after travelling " << trackLength
                       << " units";
  log_string(LogLevel::INFO, terminationLogString.str());
}

// terminate particle after leaving magnetic field
void
ParticleSimulation::terminate_particle_lost(ParticleBase & particle)
{
  std::stringstream terminationLogString;
  double heatflux = 0.0;
  vtkInterface->write_particle_track(branchLostPart, heatflux);
  _integrator->count_particle(particle.parent_entity_handle(), terminationState::LOST, heatflux);
  terminationLogString << "Particle leaving magnetic field after travelling " << trackLength
                       << " units";
  log_string(LogLevel::INFO, terminationLogString.str());
}

// terminate particle after travelling user specified max distance
void
ParticleSimulation::terminate_particle_maxlength(ParticleBase & particle)
{
  std::stringstream terminationLogString;
  double heatflux = 0.0;
  vtkInterface->write_particle_track(branchMaxLengthPart, heatflux);
  _integrator->count_particle(particle.parent_entity_handle(), terminationState::MAXLENGTH,
                              heatflux);
  terminationLogString << "Fieldline trace reached maximum length before intersection";
}

// Get triangles from the surface(s) of interest
moab::Range
ParticleSimulation::select_target_surface()
{

  moab::Range targetFacets; // range containing all of the triangles in the
                            // surface of interest
  // can specify particular surfaces of interest
  EntityHandle targetSurf; // surface of interest

  std::unordered_map<EntityHandle, moab::Range> surfFacets;
  moab::Range surfFacetsKeys;
  int numTargetFacets = 0;

  if (!vectorOfTargetSurfs.empty())
  {
    std::stringstream surfaceIDsOut;
    surfaceIDsOut << "Surface IDs provided. Number of facets in each surface "
                     "with global IDs: ";
    for (auto & surfID : vectorOfTargetSurfs)
    {
      targetSurf = DAG->entity_by_id(2, surfID); // surface of interest
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
    surfaceIDsOut << "No surface IDs provided. Launching from all surfaces by "
                     "default. Facets per "
                     "surface ID: ";

    for (auto & surfEH : surfsList)
    {
      DAG->moab_instance()->get_entities_by_type(surfEH, MBTRI, surfFacets[surfEH]);
      surfFacetsKeys.insert(surfEH);
      if (rank == 0)
      {
        surfaceIDsOut << "[" << surfEH << "] "
                      << " - " << surfFacets[surfEH].size() << ",";
      }
      numTargetFacets += surfFacets[surfEH].size();
      targetFacets.merge(surfFacets[surfEH]);
    }
    log_string(LogLevel::WARNING, surfaceIDsOut.str());
    log_string(LogLevel::WARNING, "WARNING - Will take a significant amount of "
                                  "time for larger geometry sets");
  }

  targetNumFacets = targetFacets.size();
  int totalParticles =
      (particleLaunchPos == "mc") ? nParticlesPerFacet * targetNumFacets : targetNumFacets;
  if (rank == 0)
  {
    LOG_WARNING << "Total Number of Triangles particles launched from = " << numTargetFacets;
    LOG_WARNING << "Number of particles launched from each facet = " << nParticlesPerFacet
                << std::endl;
    LOG_WARNING << "Total Number of Particles launched = " << totalParticles << std::endl;
  }

  return targetFacets;
}

// return total number of facets in target surface(s) of interest
int
ParticleSimulation::target_num_facets()
{
  return targetNumFacets;
}

int
ParticleSimulation::num_particles_launched()
{
  return (targetNumFacets * nParticlesPerFacet);
}

// print stats for the entire run
void
ParticleSimulation::print_particle_stats(std::array<int, 4> particleStats)
{
  unsigned int particlesCounted = 0;

  for (const auto i : particleStats)
  {
    particlesCounted += i;
  }

  if (rank == 0)
  {
    LOG_WARNING << "Total number of triangles in geometry (Shadow + Target) = "
                << totalNumberOfFacets;
    LOG_WARNING << "Number of particles launched = " << particlesCounted;
    LOG_WARNING << "Number of particles depositing power from omp = " << particleStats[0];
    LOG_WARNING << "Number of shadowed particle intersections = " << particleStats[1];
    LOG_WARNING << "Number of particles lost from magnetic domain = " << particleStats[2];
    LOG_WARNING << "Number of particles terminated upon reaching max tracking length = "
                << particleStats[3];

    LOG_WARNING << "Number of particles not accounted for = "
                << (num_particles_launched() - particlesCounted);
  }
  // std::cout << "Number of Ray fire calls = " << numberOfRayFireCalls <<
  // std::endl;
}

// run AEGIS simulation on single core
void
ParticleSimulation::Execute_serial()
{
  double mainParticleLoopStart = MPI_Wtime();

  vtkInterface->init();
  moab::ErrorCode rval;

  int start = 0;
  int end = num_particles_launched();

  std::vector<double> particleHeatfluxes;
  double parallelProfilingStart = MPI_Wtime();
  particleHeatfluxes = loop_over_particles(start, end);
  double parallelProfilingEnd = MPI_Wtime();
  std::cout << "Parallel Compute Runtime = " << parallelProfilingEnd - parallelProfilingStart
            << std::endl;

  if (particleHeatfluxes.empty())
  {
    log_string(LogLevel::ERROR, "Error - loop over facets returned no "
                                "heatfluxes, please check logfile. Exiting...");
    std::exit(EXIT_FAILURE);
  }

  std::vector<double> facetHeatfluxes(target_num_facets());
  auto particleIter = particleHeatfluxes.begin();
  for (unsigned int i = 0; i < facetHeatfluxes.size(); ++i)
  {
    facetHeatfluxes[i] = std::reduce(particleIter, (particleIter + nParticlesPerFacet));
    particleIter += nParticlesPerFacet; // step iterator by n particles
  }

  // std::vector<double> psiStart(target_num_facets());
  // std::vector<double> bnList(target_num_facets());
  // std::vector<std::vector<double>> normalList(target_num_facets());

  // for (int i = 0; i < listOfTriangles.size(); ++i)
  // {
  //   psiStart[i] = listOfTriangles[i].get_psi();
  //   bnList[i] = listOfTriangles[i].BdotN();
  //   normalList[i] = listOfTriangles[i].get_normal();
  // }

  attach_mesh_attribute("Heatflux", targetFacets, facetHeatfluxes);
  // attach_mesh_attribute("Psi Start", targetFacets, psiStart);
  // attach_mesh_attribute("BdotN", targetFacets, bnList);
  // attach_mesh_attribute("Normal", targetFacets, normalList);
  write_out_mesh(writeMesh);
  mainParticleLoopTime = MPI_Wtime() - mainParticleLoopStart;

  // write out data and print final

  std::array<int, 4> particleStats = _integrator->particle_stats();

  vtkInterface->write_multiBlockData("particle_tracks.vtm");

  print_particle_stats(particleStats);
}

void
ParticleSimulation::Execute_mpi()
{

  double mainParticleLoopStart = MPI_Wtime();

  vtkInterface->init();
  MPI_Status mpiStatus;

  int totalFacets = target_num_facets();
  int nFacetsPerProc = totalFacets / nprocs;
  int remainder = totalFacets % nprocs;

  int rootRank = 0;

  std::vector<double> qvalues;                  // qvalues buffer local to each processor
  std::vector<double> rootQvalues(totalFacets); // total qvalues buffer on root process for IO

  std::vector<int> recieveCounts(nprocs, nFacetsPerProc);
  auto recvItr = recieveCounts.begin();

  int remaindersHandled = 0;
  while (remaindersHandled < remainder) // distribute remainders across processes
  {
    *recvItr += 1;
    std::advance(recvItr, 1);
    remaindersHandled++;
    if (recvItr == recieveCounts.end())
    {
      recvItr = recieveCounts.begin();
    }
  }

  std::vector<int> displacements(nprocs, 0);
  int displ = 0;

  for (int i = 0; i < nprocs; ++i) // calculate displacements
  {
    displacements[i] = displ;
    displ += recieveCounts[i];
  }

  int startIndex = displacements[rank];
  int endIndex = startIndex + recieveCounts[rank];

  double parallelProfilingStart = MPI_Wtime();
  qvalues = loop_over_particles(startIndex, endIndex); // perform main loop
  double endTime = MPI_Wtime();

  std::array<int, 4> particleStats = _integrator->particle_stats();
  std::array<int, 4> totalParticleStats;
  if (rank != 0)
  {
    MPI_Gatherv(qvalues.data(), qvalues.size(), MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE,
                rootRank, MPI_COMM_WORLD);
    // MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, NULL, 0,
    // MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);
  }
  else
  {
    totalParticleStats = _integrator->particle_stats();
    MPI_Gatherv(qvalues.data(), qvalues.size(), MPI_DOUBLE, rootQvalues.data(),
                recieveCounts.data(), displacements.data(), MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
    // MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE,
    // rootQvalues.data(), qvalues.size(), MPI_DOUBLE, root_rank,
    // MPI_COMM_WORLD); for (auto i:rootQvalues) { temp << i << std::endl; }

    for (int i = 1; i < nprocs; ++i)
    {
      MPI_Recv(particleStats.data(), particleStats.size(), MPI_INT, i, 11, MPI_COMM_WORLD,
               &mpiStatus);
      for (int j = 0; j < particleStats.size(); ++j)
      {
        totalParticleStats[j] += particleStats[j];
      }
    }

    attach_mesh_attribute("Heatflux", targetFacets, rootQvalues);
    write_out_mesh(writeMesh, targetFacets);
    mainParticleLoopTime = MPI_Wtime() - mainParticleLoopStart;
  }

  print_particle_stats(totalParticleStats);
}

void
ParticleSimulation::implicit_complement_testing()
{

  std::vector<EntityHandle> children;
  DAG->moab_instance()->get_child_meshsets(implComplementVol, children, 1);
  std::vector<EntityHandle> vertices;

  for (const auto & i : children)
  {
    std::vector<EntityHandle> temp;
    DAG->moab_instance()->get_entities_by_type(i, MBVERTEX, vertices, false);
    vertices.insert(vertices.begin(), temp.begin(), temp.end());
  }
  std::vector<double> vertexCoords(vertices.size() * 3);
  std::vector<double> vertexCoordsFlux(vertices.size() * 3);

  DAG->moab_instance()->get_coords(&vertices[0], vertices.size(), vertexCoords.data());

  // DAG->moab_instance()->list_entities(children);

  std::ofstream implicitComplCoordsxyz("implcit_complement_xyz.txt");
  std::ofstream implicitComplCoordsFlux("implcit_complement_flux.txt");
  std::vector<double> temp1, temp2;

  for (int i = 0; i < vertexCoords.size(); i += 3)
  {
    implicitComplCoordsxyz << vertexCoords[i] << " " << vertexCoords[i + 1] << " "
                           << vertexCoords[i + 2] << std::endl;
    temp1 =
        CoordTransform::cart_to_polar(vertexCoords[i], vertexCoords[i + 1], vertexCoords[i + 2]);
    temp2 = CoordTransform::polar_to_flux(temp1, equilibrium);
    vertexCoordsFlux[i] = temp2[0];
    vertexCoordsFlux[i + 1] = temp2[1];
    vertexCoordsFlux[i + 2] = temp2[2];
    implicitComplCoordsFlux << temp2[0] << " " << temp2[1] << " " << temp2[2] << std::endl;
  }

  // int ctr = 101;
  // for (const auto &i:vertexCoordsFlux)
  // {
  //   implicitComplCoordsFlux << i[0] << " " << i[1] << " " << i[2] <<
  //   std::endl; EntityHandle vertexHandle = ctr;
  //   fluxDAG->moab_instance()->create_vertex(i.data(), vertexHandle);
  //   EntityHandle tri_conn[] = {i[0], vertex1, vertex2, vertex3};
  //   EntityHandle quad_handle = 0;
  //   create_element( MeshQuad, quad_conn, 4, quad_handle );
  //   ctr +=1;
  // }

  std::cout << "Number of Nodes in implicit complement = " << vertexCoords.size() << std::endl;

  return;
}

// create a moab tag and attach that data to each element in the moab::Range
// provided
void
ParticleSimulation::attach_mesh_attribute(const std::string & tagName, moab::Range & entities,
                                          std::vector<double> & dataToAttach)
{
  moab::Tag tag;
  double tagValue;

  DAG->moab_instance()->tag_get_handle(tagName.c_str(), 1, MB_TYPE_DOUBLE, tag,
                                       MB_TAG_CREAT | MB_TAG_DENSE, &tagValue);
  DAG->moab_instance()->tag_set_data(tag, entities, &dataToAttach[0]);
}

// overload for attaching vectors to mesh entitities
void
ParticleSimulation::attach_mesh_attribute(const std::string & tagName, moab::Range & entities,
                                          std::vector<std::vector<double>> & dataToAttach)
{
  moab::Tag tag;
  double tagValue;

  DAG->moab_instance()->tag_get_handle(tagName.c_str(), 3, MB_TYPE_DOUBLE, tag,
                                       MB_TAG_CREAT | MB_TAG_DENSE, &tagValue);

  DAG->moab_instance()->tag_set_data(tag, entities, &dataToAttach[0][0]);
}
// write out the mesh and target mesh with attributed data
// meshWriteOptions::FULL -- Write out the entire mesh
// meshWriteOptions::TARGET -- Only write out target mesh with heatfluxes
// (useful if full mesh particularly large) meshWriteOptions::BOTH -- Both FULL
// and TARGET meshes are written out meshWriteOptions::PARTIAL -- Write out
// specified entities default -- Full mesh is written
void
ParticleSimulation::write_out_mesh(meshWriteOptions option, moab::Range rangeofEntities)
{
  // mesh_coord_transform(coordinateSystem::CARTESIAN);
  // get target meshset for aegis_target outputs
  DAG->remove_graveyard();
  EntityHandle targetMeshset;
  DAG->moab_instance()->create_meshset(MESHSET_SET, targetMeshset);
  DAG->moab_instance()->add_entities(targetMeshset, targetFacets);

  // remove file extension and path from output file string
  std::string aegisOut = dagmcInputFile;
  aegisOut = aegisOut.substr(0, aegisOut.find_last_of("."));
  aegisOut = aegisOut.substr(aegisOut.find_last_of("/") + 1, aegisOut.length());

  std::string aegisOutFull = "aegis_full_" + aegisOut + ".vtk";
  std::string aegisOutTarget = "aegis_target_" + aegisOut + ".vtk";
  std::string aegisOutPartial = "aegis_partial_" + aegisOut + ".vtk";

  switch (option)
  {
    case meshWriteOptions::FULL:
      DAG->write_mesh(aegisOutFull.c_str(), 1);
      break;

    case meshWriteOptions::TARGET:
      DAG->moab_instance()->write_mesh(aegisOutTarget.c_str(), &targetMeshset, 1);
      break;

    case meshWriteOptions::BOTH:
      DAG->moab_instance()->write_mesh(aegisOutTarget.c_str(), &targetMeshset, 1);
      DAG->write_mesh(aegisOutFull.c_str(), 1);
      break;

    case meshWriteOptions::PARTIAL:
      if (!rangeofEntities.empty())
      {
        EntityHandle meshset;
        DAG->moab_instance()->create_meshset(MESHSET_SET, meshset);
        DAG->moab_instance()->add_entities(meshset, rangeofEntities);
        DAG->moab_instance()->write_mesh(aegisOutPartial.c_str(), &meshset, 1);
      }
      else
      {
        log_string(LogLevel::ERROR, "No meshsets provided for partial mesh write "
                                    "out. Defaulting to full mesh");
        DAG->write_mesh("aegis_full.vtk", 1);
      }
      break;

    default: // default to full mesh and target
      DAG->moab_instance()->write_mesh("aegis_out_target.vtk", &targetMeshset, 1);
      DAG->write_mesh("aegis_full.vtk", 1);
      break;
  }
}

void
ParticleSimulation::mesh_coord_transform(coordinateSystem coordSys)
{
  std::vector<double> temp(3);
  int counter = 0;
  switch (coordSys)
  {
    case coordinateSystem::CARTESIAN:
      for (int i = 0; i < nodeCoords.size(); i += 3)
      {
        temp = {nodeCoords[i], nodeCoords[i + 1], nodeCoords[i + 2]};
        temp = CoordTransform::polar_to_cart(temp);
        DAG->moab_instance()->set_coords(&nodesList[counter], 1, temp.data());
        counter += 1;
      }
      break;

    case coordinateSystem::POLAR:
      for (int i = 0; i < nodeCoords.size(); i += 3)
      {
        temp = {nodeCoords[i], nodeCoords[i + 1], nodeCoords[i + 2]};
        temp = CoordTransform::cart_to_polar(temp);
        DAG->moab_instance()->set_coords(&nodesList[counter], 1, temp.data());
        counter += 1;
      }
      break;

    case coordinateSystem::FLUX:
      // for (int i=0; i<nodeCoords.size(); i+=3)
      // {
      //   temp = {nodeCoords[i], nodeCoords[i+1], nodeCoords[i+2]};
      //   temp = CoordTransform::cart_to_polar(temp);
      //   temp = CoordTransform::polar_to_flux(temp, equilibrium);
      //   DAG->moab_instance()->set_coords(&nodesList[counter], 1, temp.data());
      //   counter +=1;
      // }
      break;

    default:
      std::cout << "No coordinate system provided for mesh coordinate transform" << std::endl;
      break;
  }
}

void
ParticleSimulation::write_particle_launch_positions(std::vector<double> & particleHeatfluxes)
{
  std::ofstream particleLaunchPosOut("particle_launch_positions.txt");
  particleLaunchPosOut << "X,Y,Z,Heatflux_per_particle" << std::endl;
  for (int i = 0; i < listOfParticles.size(); ++i)
  {
    particleLaunchPosOut << listOfParticles[i].launchPos[0] << ","
                         << listOfParticles[i].launchPos[1] << ","
                         << listOfParticles[i].launchPos[2] << "," << particleHeatfluxes[i] << "\n";
  }
  particleLaunchPosOut << std::flush;
}

// return timings for main parallel loop, mesh setup and mesh write out
std::map<std::string, double>
ParticleSimulation::get_profiling_times()
{
  std::map<std::string, double> profilingTimes;
  profilingTimes.insert(std::make_pair("DAGMC Mesh read runtime = ", dagmcMeshReadTime));
  profilingTimes.insert(
      std::make_pair("Preparing surfaces for particles runtime = ", prepSurfacesTime));
  profilingTimes.insert(
      std::make_pair("Pool of particles generation runtime = ", setupArrayOfParticlesTime));
  profilingTimes.insert(
      std::make_pair("Main particle tracking loop runtime = ", mainParticleLoopTime));
  profilingTimes.insert(std::make_pair("Mesh Write out runtime = ", aegisMeshWriteTime));

  return profilingTimes;
}

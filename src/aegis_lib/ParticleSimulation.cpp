#include "ParticleSimulation.h"
#include <mpi.h>

// setup AEGIS simulation
ParticleSimulation::ParticleSimulation(std::shared_ptr<InputJSON> configFile,
                                       std::shared_ptr<EquilData> equil)
{
  set_mpi_params();

  JSONsettings = configFile;

  read_params(JSONsettings);
  equilibrium = equil;
  DAG = std::make_unique<moab::DagMC>();
  vtkInterface = std::make_unique<VtkInterface>(JSONsettings);

  init_geometry();
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

  else if (exeType == "mpi_padded" || exeType == "padded")
  {
    log_string(LogLevel::WARNING, "Executing padded MPI with no load balancing...");
    Execute_padded_mpi();
  }

  else if (exeType == "mpi_dynamic" || exeType == "dynamic")
  {
    log_string(LogLevel::WARNING, "Executing dynamic MPI load balancing...");
    Execute_dynamic_mpi();
  }

  else
  {
    log_string(LogLevel::FATAL, "No suitable execution type ('serial', 'mpi', 'mpi_padded', "
                                "'mpi_dynamic') provided exiting...");
    std::exit(1);
  }
}

// perform AEGIS simulation solve with dynamic mpi load balancing
void
ParticleSimulation::Execute_dynamic_mpi()
{

  MPI_Status mpiStatus;

  vtkInterface->init();
  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  const int noMoreWork = -1;

  std::vector<double> handlerQVals(num_facets());
  // int counter = 0;
  setup_sources();

  if (rank == 0)
  { // handler process...
    handler(handlerQVals);

    // for (const auto i:handlerQVals)
    // {
    //   file << i << "\n";
    //   if (i > 0 ) { depositCounter +=1; }
    // }
    // file << std::endl;
    // file << "Number of depositing particles = " << depositCounter <<
    // std::endl;
  }

  else
  { // worker processes...
    worker();
  }

  std::array<int, 5> particleStats = integrator->particle_stats();
  std::array<int, 5> totalParticleStats;

  if (rank != 0)
  {
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);
  }
  else
  {
    totalParticleStats = integrator->particle_stats();
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
    attach_mesh_attribute("Heatflux", targetFacets, handlerQVals);
    write_out_mesh(meshWriteOptions::BOTH, targetFacets);
  }

  mpi_particle_stats();

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
  int counter = 0;

  std::vector<double> workerQVals(dynamicTaskSize);
  std::cout << "Dynamic task scheduling with " << (nprocs - 1) << " processes, each handling "
            << dynamicTaskSize << " facets" << std::endl;
  int particlesHandled = 0;
  int batchesComplete = 0;
  int totalBatches = num_facets() / dynamicTaskSize;

  for (int procID = 1; procID < nprocs; procID++)
  { // Send the initial indexes for each process statically
    MPI_Send(&particlesHandled, 1, MPI_INT, procID, procID, MPI_COMM_WORLD);
    particlesHandled += dynamicTaskSize;
    batchesComplete += 1;
  }

  int recvCount = 0;

  int inactiveWorkers = 0;

  do
  {
    int avaialbleProcess = 0;
    MPI_Recv(&avaialbleProcess, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
             &status); // get avaialble process
    if (particlesHandled < handlerQVals.size())
    { // send that process the next index to start on
      // in the total list
      MPI_Send(&particlesHandled, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD);

      MPI_Recv(workerQVals.data(), workerQVals.size(), MPI_DOUBLE, MPI_ANY_SOURCE, mpiDataTag,
               MPI_COMM_WORLD, &status);
      MPI_Recv(&workerStartIndex, 1, MPI_INT, MPI_ANY_SOURCE, mpiIndexTag, MPI_COMM_WORLD, &status);
      // printf("Time taken to recieve next data = %f \n", MPI_Wtime());

      particlesHandled += dynamicTaskSize;
      batchesComplete += 1;
      // progress indicator

      if (batchesComplete == totalBatches / 4)
      {
        std::cout << "25% ... " << std::flush;
      }
      else if (batchesComplete == totalBatches / 2)
      {
        std::cout << "50% ... " << std::flush;
      }
      else if (batchesComplete == 3 * totalBatches / 4)
      {
        std::cout << "75% ... " << std::flush;
      }

      counter += workerQVals.size();
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

  // std::ofstream file;
  // std::stringstream fileName;
  // fileName << "rank" << rank << ".txt";
  // file.open(fileName.str());

  std::vector<double> workerQVals;
  int startIndex = 0;
  MPI_Recv(&startIndex, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);

  int start = startIndex;
  int end = start + dynamicTaskSize;

  workerQVals = loop_over_facets(start, end); // process initial facets
  // file << "Loop over [" << start << ":" << end << "] facets: \n";
  // index = start;
  // for (auto i:workerQVals)
  // {
  //   file << "[" << index+1 << "] " <<  i << "\n";
  //   if (i > 0) { counterDeposit +=1; }
  //   index +=1;
  // }
  // file << std::endl;

  MPI_Send(workerQVals.data(), workerQVals.size(), MPI_DOUBLE, 0, mpiDataTag, MPI_COMM_WORLD);
  MPI_Send(&start, 1, MPI_INT, 0, mpiIndexTag, MPI_COMM_WORLD);
  do
  {
    MPI_Send(&rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Recv(&startIndex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

    if (startIndex + dynamicTaskSize > num_facets())
    {
      dynamicTaskSize = num_facets() - startIndex;
    }

    start = startIndex;
    end = start + dynamicTaskSize;

    if (startIndex != noMoreWork)
    { // processing the next set of available work
      // that has been dynamically allocated
      workerQVals = loop_over_facets(start, end);
      // send local (to worker) array back to handler process
      MPI_Send(workerQVals.data(), workerQVals.size(), MPI_DOUBLE, 0, mpiDataTag, MPI_COMM_WORLD);
      MPI_Send(&start, 1, MPI_INT, 0, mpiIndexTag, MPI_COMM_WORLD);

      // file << "Loop over [" << start << ":" << end << "] facets: \n";
      // index = start;
      // for (auto i:workerQVals)
      // {
      //   file << "[" << index+1 << "] " <<  i << "\n";
      //   if (i > 0) { counterDeposit +=1; }
      //   index +=1;
      // }
      // file << std::endl;
    }
  } while (startIndex != noMoreWork); // while there is work available

  // file << "Particle Stats:" << std::endl;
  // file << "Depositing particles = " << counterDeposit << std::endl;
}

// Read parameters from aegis_settings.json config file
void
ParticleSimulation::read_params(const std::shared_ptr<InputJSON> & inputs)
{

  json aegisNamelist;
  if (inputs->data.contains("aegis_params"))
  {
    aegisNamelist = inputs->data["aegis_params"];
    dagmcInputFile = aegisNamelist["DAGMC"];
    trackStepSize = aegisNamelist["step_size"];
    maxTrackSteps = aegisNamelist["max_steps"];
    particleLaunchPos = aegisNamelist["launch_pos"];
    noMidplaneTermination = aegisNamelist["force_no_deposition"];
    coordInputStr = aegisNamelist["coordinate_system"];
    exeType = aegisNamelist["execution_type"];
    dynamicTaskSize = aegisNamelist["dynamic_task_size"];

    select_coordinate_system();

    if (aegisNamelist.contains("target_surfs"))
    {
      for (auto i : aegisNamelist["target_surfs"])
      {
        vectorOfTargetSurfs.push_back(i);
      }
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
  string_to_lowercase(coordInputStr);
  if (coordInputStr == "cart" || coordInputStr == "cartesian" || coordInputStr == "xyz")
  {
    coordSys = coordinateSystem::CARTESIAN;
    log_string(LogLevel::WARNING, "Tracking in CARTESIAN coordinates...");
  }
  else if (coordInputStr == "pol" || coordInputStr == "polar" || coordInputStr == "rz")
  {
    coordSys = coordinateSystem::POLAR;
    log_string(LogLevel::WARNING, "Tracking in POLAR coordinates...");
  }
  else if (coordInputStr == "flux" || coordInputStr == "psi")
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

// initialise CAD geometry for AEGIS and magnetic field equilibrium for CAD
// Geometry
void
ParticleSimulation::init_geometry()
{
  // setup dagmc instance
  DAG->load_file(dagmcInputFile.c_str());
  DAG->init_OBBTree();
  DAG->setup_geometry(surfsList, volsList);
  DAG->moab_instance()->get_entities_by_type(0, MBTRI, facetsList);

  DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, nodesList, true);
  nodeCoords.resize(nodesList.size() * 3);
  DAG->moab_instance()->get_coords(&nodesList[0], nodesList.size(), nodeCoords.data());

  if (rank == 0)
  {
    std::cout << "Number of Triangles in Geometry " << facetsList.size() << std::endl;
  }
  implComplementVol = DAG->entity_by_index(3, volsList.size());

  if (!DAG->is_implicit_complement(implComplementVol))
  {
    std::cout << "Particle not in implicit complement. Check volumes." << std::endl;
  }

  // setup B Field data

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

  // if flux selected default to cartesian
  if (coordSys == coordinateSystem::FLUX)
  {
    coordSys = coordinateSystem::CARTESIAN;
    log_string(LogLevel::WARNING, "Flux coords tracking not currently "
                                  "implemented. Defaulting to CARTESIAN...");
  }
}

// loop over facets in target surfaces
std::vector<double>
ParticleSimulation::loop_over_facets(int startFacet, int endFacet)
{

  // transform coordinate system of mesh
  if (coordSys != coordinateSystem::CARTESIAN)
  {
    mesh_coord_transform(coordSys);
  }

  double startTime = MPI_Wtime();
  std::vector<double> heatfluxVals;

  std::vector<double> triCoords(9);
  std::vector<double> triA(3), triB(3), triC(3);

  for (int i = startFacet; i < endFacet; ++i)
  { // loop over all facets
    auto tri = listOfTriangles[i];

    if (i >= targetFacets.size()) // handle padded values
    {
      tri.update_heatflux(-1);
      integrator->count_particle(0, terminationState::PADDED, 0.0);
      continue;
    }

    auto particle = std::make_unique<ParticleBase>(coordSys);
    particle->set_pos(tri.launch_pos());

    integrator->set_launch_position(tri.entity_handle(), tri.launch_pos());

    particle->check_if_in_bfield(equilibrium); // if out of bounds skip to next triangle
    if (particle->outOfBounds)
    {
      log_string(LogLevel::INFO, "Particle start is out of magnetic field bounds. Skipping to "
                                 "next triangle. Check "
                                 "correct eqdsk is being used for the given geometry");
      continue;
    }

    vtkInterface->init_new_vtkPoints();
    vtkInterface->insert_next_point_in_track(particle->launchPos);

    particle->set_dir(equilibrium);

    particle->align_dir_to_surf(tri.BdotN());

    // Start ray tracing
    DagMC::RayHistory history;
    trackLength = trackStepSize;

    particle->update_vectors(trackStepSize, equilibrium);

    vtkInterface->insert_next_point_in_track(particle->posXYZ);

    terminationState particleState = loop_over_particle_track(tri, particle, history);

    // if particle not depositing set heatflux to 0.0
    if (particleState != terminationState::DEPOSITING)
    {
      tri.update_heatflux(0.0);
    }

    heatfluxVals.push_back(tri.heatflux());
  }

  double endTime = MPI_Wtime();
  bool profileDiagnostics = true;
  if (profileDiagnostics)
  {
    printf("Loop over facets [%d:%d] from rank[%d] time taken = %fs \n", startFacet, endFacet, rank,
           endTime - startTime);
    fflush(stdout);
  }
  return heatfluxVals;
}

void
ParticleSimulation::setup_sources()
{

  std::vector<double> triCoords(9);
  std::vector<double> triA(3), triB(3), triC(3);
  for (const auto & facet : targetFacets)
  {
    std::vector<moab::EntityHandle> triNodes;
    DAG->moab_instance()->get_adjacencies(&facet, 1, 0, false, triNodes);
    DAG->moab_instance()->get_coords(&triNodes[0], triNodes.size(), triCoords.data());

    for (int j = 0; j < 3; j++)
    {
      triA[j] = triCoords[j];
      triB[j] = triCoords[j + 3];
      triC[j] = triCoords[j + 6];
    }

    TriangleSource tri(triA, triB, triC, facet, particleLaunchPos);
    tri.set_heatflux_params(equilibrium, "exp");
    listOfTriangles.push_back(tri);
  }
}

// loop over single particle track from launch to termination
terminationState
ParticleSimulation::loop_over_particle_track(TriangleSource & tri,
                                             std::unique_ptr<ParticleBase> & particle,
                                             DagMC::RayHistory & history)
{
  double euclidDistToNextSurf = 0.0;
  nextSurf = 0;

  for (int step = 0; step < maxTrackSteps; ++step)
  {
    trackLength += trackStepSize;
    iterationCounter++;

    DAG->ray_fire(implComplementVol, particle->posXYZ.data(), particle->dir.data(), nextSurf,
                  nextSurfDist, &history, trackStepSize, rayOrientation);
    numberOfRayFireCalls++;

    if (nextSurf != 0)
    {
      particle->update_vectors(nextSurfDist); // update position to surface intersection point
      terminate_particle(tri, history, terminationState::SHADOWED);
      return terminationState::SHADOWED;
    }
    else
    {
      particle->update_vectors(trackStepSize); // update position by stepsize
    }

    vtkInterface->insert_next_point_in_track(particle->posXYZ);
    particle->check_if_in_bfield(equilibrium);
    if (particle->outOfBounds)
    {
      terminate_particle(tri, history, terminationState::LOST);
      return terminationState::LOST;
    }

    else
    {
      particle->set_dir(equilibrium);
      particle->align_dir_to_surf(tri.BdotN());
    }
    particle->check_if_midplane_crossed(equilibrium->get_midplane_params());

    if (particle->atMidplane != 0 && !noMidplaneTermination)
    {
      terminate_particle(tri, history, terminationState::DEPOSITING);
      return terminationState::DEPOSITING;
    }

    if (step == (maxTrackSteps - 1))
    {
      terminate_particle(tri, history, terminationState::MAXLENGTH);
      return terminationState::MAXLENGTH;
    }
  }

  return terminationState::LOST; // default return lost particle
}

// terminate particle depending on 1 of 4 termination states:
//
// DEPOSITING - Particle reaches outer midplane and deposits power on facet.
// Heatflux = Q SHADOWED - Particle hits another piece of geometry. Heatflux = 0
// LOST - Particle leaves magnetic field so trace stops. Heatflux = 0
// MAX LENGTH - Particle reaches maximum user set length before anything else.
// Heatflux = 0
void
ParticleSimulation::terminate_particle(TriangleSource & tri, DagMC::RayHistory & history,
                                       terminationState termination)
{
  double heatflux;
  std::stringstream ss;

  switch (termination)
  {
    case terminationState::DEPOSITING:
      heatflux = tri.heatflux();
      vtkInterface->write_particle_track(branchDepositingPart, heatflux);
      integrator->count_particle(tri.entity_handle(), termination, heatflux);
      log_string(LogLevel::INFO, "Midplane reached. Depositing power");
      break;

    case terminationState::SHADOWED:
      heatflux = 0.0;
      history.get_last_intersection(intersectedFacet);
      vtkInterface->write_particle_track(branchShadowedPart, heatflux);
      integrator->count_particle(tri.entity_handle(), termination, heatflux);
      ss << "Surface " << nextSurf << " hit after travelling " << trackLength << " units";
      log_string(LogLevel::INFO, ss.str());
      break;

    case terminationState::LOST:
      heatflux = 0.0;
      vtkInterface->write_particle_track(branchLostPart, heatflux);
      integrator->count_particle(tri.entity_handle(), termination, heatflux);
      log_string(LogLevel::INFO, "TRACE STOPPED BECAUSE LEAVING MAGNETIC FIELD");
      break;

    case terminationState::MAXLENGTH:
      heatflux = 0.0;
      vtkInterface->write_particle_track(branchMaxLengthPart, heatflux);
      integrator->count_particle(tri.entity_handle(), termination, heatflux);
      log_string(LogLevel::INFO, "Fieldline trace reached maximum length before intersection");
  }
}

void
ParticleSimulation::ray_hit_on_launch(std::unique_ptr<ParticleBase> & particle,
                                      DagMC::RayHistory & history)
{
  DAG->ray_fire(implComplementVol, particle->launchPos.data(), particle->dir.data(), nextSurf,
                nextSurfDist, &history, trackStepSize, rayOrientation);
  if (nextSurf != 0)
  {
    history.get_last_intersection(intersectedFacet);
    EntityHandle nextVolEH;
    DAG->next_vol(nextSurf, implComplementVol, nextVolEH);
    if (rank == 0)
    {
      LOG_INFO << "---- RAY HIT ON LAUNCH [" << facetCounter << "] ----";
    }
  }
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
  if (rank == 0)
  {
    LOG_WARNING << "Total Number of Triangles rays launched from = " << numTargetFacets;
  }

  numFacets = targetFacets.size();
  return targetFacets;
}

// return total number of facets in target surface(s) of interest
int
ParticleSimulation::num_facets()
{
  return numFacets;
}

// print stats for the entire run
void
ParticleSimulation::print_particle_stats(std::array<int, 5> particleStats)
{
  unsigned int particlesCounted = 0;

  for (const auto i : particleStats)
  {
    particlesCounted += i;
  }

  particlesCounted -= particleStats[4]; // remove padded particles

  if (rank == 0)
  {
    LOG_WARNING << "Number of particles launched = " << particlesCounted << std::endl;
    LOG_WARNING << "Number of particles depositing power from omp = " << particleStats[0]
                << std::endl;
    LOG_WARNING << "Number of shadowed particle intersections = " << particleStats[1] << std::endl;
    LOG_WARNING << "Number of particles lost from magnetic domain = " << particleStats[2]
                << std::endl;
    LOG_WARNING << "Number of particles terminated upon reaching max tracking length = "
                << particleStats[3] << std::endl;
    LOG_WARNING << "Number of padded null particles = " << particleStats[4] << std::endl;
    LOG_WARNING << "Number of particles not accounted for = " << (num_facets() - particlesCounted)
                << std::endl;
  }
  // std::cout << "Number of Ray fire calls = " << numberOfRayFireCalls <<
  // std::endl;
}

// print individual particle stats for each MPI rank
void
ParticleSimulation::mpi_particle_stats()
{

  for (int i = 1; i < nprocs; ++i)
  {
    if (rank == i)
    {
      std::cout << std::endl
                << "process " << i << " has the following particle stats:" << std::endl;
      std::array localRankParticleStats = integrator->particle_stats();

      std::cout << "DEPOSITING - " << localRankParticleStats[0] << std::endl;
      std::cout << "SHADOWED - " << localRankParticleStats[1] << std::endl;
      std::cout << "LOST - " << localRankParticleStats[2] << std::endl;
      std::cout << "MAX LENGTH - " << localRankParticleStats[3] << std::endl;

      int totalParticlesHandled = localRankParticleStats[0] + localRankParticleStats[1] +
                                  localRankParticleStats[2] + localRankParticleStats[3];

      std::cout << "TOTAL - " << totalParticlesHandled << std::endl;
    }
  }
}

// run AEGIS simulation on single core
void
ParticleSimulation::Execute_serial()
{

  vtkInterface->init();

  moab::ErrorCode rval;

  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  // implicit_complement_testing();

  int start = 0;
  int end = targetFacets.size();

  std::vector<double> qvalues;
  setup_sources();
  qvalues = loop_over_facets(start, end);
  if (qvalues.empty())
  {
    log_string(LogLevel::ERROR, "Error - loop over facets returned no "
                                "heatfluxes, please check logfile. Exiting...");
    std::exit(EXIT_FAILURE);
  }

  std::vector<double> psiStart(num_facets());
  std::vector<double> bnList(num_facets());
  std::vector<std::vector<double>> normalList(num_facets());
  for (int i = 0; i < listOfTriangles.size(); ++i)
  {
    psiStart[i] = listOfTriangles[i].get_psi();
    bnList[i] = listOfTriangles[i].BdotN();
    normalList[i] = listOfTriangles[i].get_normal();
  }

  attach_mesh_attribute("Heatflux", targetFacets, qvalues);
  attach_mesh_attribute("Psi Start", targetFacets, psiStart);
  attach_mesh_attribute("BdotN", targetFacets, bnList);
  // attach_mesh_attribute("Normal", targetFacets, normalList);
  write_out_mesh(meshWriteOptions::BOTH);

  // write out data and print final

  std::array<int, 5> particleStats = integrator->particle_stats();

  vtkInterface->write_multiBlockData("particle_tracks.vtm");

  print_particle_stats(particleStats);
}

void
ParticleSimulation::Execute_padded_mpi()
{

  vtkInterface->init();
  MPI_Status mpiStatus;

  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  int totalFacets = num_facets();
  int remainder = totalFacets % nprocs;
  int nPadded = 0;
  if (remainder != 0)
  {
    int paddedTotalFacets = num_facets() + (nprocs - remainder);
    nPadded = paddedTotalFacets - totalFacets;
    totalFacets = paddedTotalFacets;
  }
  int nFacetsPerProc = totalFacets / nprocs;
  int startFacet = rank * nFacetsPerProc;
  int endFacet = startFacet + nFacetsPerProc;
  int rootRank = 0;

  std::vector<double> qvalues;                  // qvalues buffer local to each processor
  std::vector<double> rootQvalues(totalFacets); // total qvalues buffer on root process for IO
  setup_sources();

  qvalues = loop_over_facets(startFacet, endFacet); // perform main loop

  std::array<int, 5> particleStats = integrator->particle_stats();
  std::array<int, 5> totalParticleStats;
  std::vector<double> nonRootQvals(nFacetsPerProc);
  if (rank != 0)
  {
    MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, rootRank,
               MPI_COMM_WORLD);
    MPI_Send(particleStats.data(), particleStats.size(), MPI_INT, 0, 11, MPI_COMM_WORLD);
  }
  else
  {
    totalParticleStats = integrator->particle_stats();
    MPI_Gather(qvalues.data(), qvalues.size(), MPI_DOUBLE, rootQvalues.data(), qvalues.size(),
               MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
    if (rootQvalues.size() > num_facets())
    {
      rootQvalues.resize(rootQvalues.size() - nPadded);
    }

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
    write_out_mesh(meshWriteOptions::BOTH, targetFacets);
  }

  mpi_particle_stats();
  print_particle_stats(totalParticleStats);
}

void
ParticleSimulation::Execute_mpi()
{
  startTime = MPI_Wtime();
  vtkInterface->init();
  MPI_Status mpiStatus;

  targetFacets = select_target_surface();
  integrator = std::make_unique<SurfaceIntegrator>(targetFacets);

  int totalFacets = num_facets();
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

  int startFacet = displacements[rank];
  int endFacet = startFacet + recieveCounts[rank];

  startTime = MPI_Wtime();
  setup_sources();

  qvalues = loop_over_facets(startFacet, endFacet); // perform main loop
  endTime = MPI_Wtime();

  std::array<int, 5> particleStats = integrator->particle_stats();
  std::array<int, 5> totalParticleStats;
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
    totalParticleStats = integrator->particle_stats();
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
    write_out_mesh(meshWriteOptions::BOTH, targetFacets);
  }

  mpi_particle_stats();
  print_particle_stats(totalParticleStats);

  endTime = MPI_Wtime();
}

// get surfaces attributed to the aegis_target group set in cubit
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
    temp1 = CoordTransform::cart_to_polar(vertexCoords[i], vertexCoords[i + 1], vertexCoords[i + 2],
                                          "forwards");
    temp2 = CoordTransform::polar_to_flux(temp1, "forwards", equilibrium);
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

void
ParticleSimulation::attach_mesh_attribute(const std::string & tagName, moab::Range & entities,
                                          std::vector<std::vector<double>> & dataToAttach)
{
  moab::Tag tag;
  double tagValue;

  DAG->moab_instance()->tag_get_handle(tagName.c_str(), 1, MB_TYPE_OPAQUE, tag,
                                       MB_TAG_CREAT | MB_TAG_DENSE, &tagValue);
  DAG->moab_instance()->tag_set_data(tag, entities, &dataToAttach[0]);
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
      if (!rangeofEntities.empty())
      {
        EntityHandle meshset;
        DAG->moab_instance()->create_meshset(MESHSET_SET, meshset);
        DAG->moab_instance()->add_entities(meshset, rangeofEntities);
        DAG->moab_instance()->write_mesh("aegis_partial.vtk", &meshset, 1);
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
        temp = CoordTransform::cart_to_polar(temp, "backwards");
        DAG->moab_instance()->set_coords(&nodesList[counter], 1, temp.data());
        counter += 1;
      }
      break;

    case coordinateSystem::POLAR:
      for (int i = 0; i < nodeCoords.size(); i += 3)
      {
        temp = {nodeCoords[i], nodeCoords[i + 1], nodeCoords[i + 2]};
        temp = CoordTransform::cart_to_polar(temp, "forwards");
        DAG->moab_instance()->set_coords(&nodesList[counter], 1, temp.data());
        counter += 1;
      }
      break;

    case coordinateSystem::FLUX:
      // for (int i=0; i<nodeCoords.size(); i+=3)
      // {
      //   temp = {nodeCoords[i], nodeCoords[i+1], nodeCoords[i+2]};
      //   temp = CoordTransform::cart_to_polar(temp, "forwards");
      //   temp = CoordTransform::polar_to_flux(temp, "forwards", equilibrium);
      //   DAG->moab_instance()->set_coords(&nodesList[counter], 1, temp.data());
      //   counter +=1;
      // }
      break;

    default:
      std::cout << "No coordinate system provided for mesh coordinate transform" << std::endl;
      break;
  }
}
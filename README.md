# Aegis 
Aegis is a Monte-Carlo particle tracking tool for calculating power deposition loads due to particles on surface accurate representations of fusion plasma First Wall Shielding. Aegis is part of Aurora-Multiphysics simulation code group. It makes use of DAGMC/Double-Down as a ray tracing tool that will be used to calculate neutral particle trajectories. Charged particle trajectories are tracked as following magnetic field lines calculated from toroidal flux surface data.

The following dependencies are **required**:
- **Embree v3.6.1** (Intel Embree Ray Tracer)  
- **Double-Down v1.0.0** (A double precision interface to Embree) 
- **MOAB Version 5.2.0** (Mesh Oriented datABase) 
- **DAGMC** (Direct Accelerated Geometry Monte Carlo code)
- **VTK** (Visualisation ToolKit) - For producing more advaced visualisations within ParaView (may become optional in the future)

# Building the repo
Embree should be built with the following additional flags:

    cmake -DEMBREE_ISPC_SUPPORT=0
Double-Down should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DEMBREE_DIR=\${EMBREE_DIR}
MOAB should be configured with the following additional flags:

    ./src/configure --with-hdf5=\${HDF5_DIR} --enable-optimize --enable-shared --disable-debug
DAGMC should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DDOUBLE_DOWN=on -DDOUBLE_DOWN_DIR=\${DOUBLEDOWN_DIR} -DBUILD_TALLY=ON

This repo's submodules (`nlohmann_json`) should be installed with:

    git submodule update --init --recursive

And finally AEGIS can be configured and built with:

    cmake -DDAGMC_DIR=\${DAGMC_DIR} 
    make 
    make test

This produces an `aegis` executable located in `bin`. The `inres1` case can be run as an example by navigating to the directory and calling 

    `../bin/aegis aegis_settings.json`

AEGIS requires a json config file to set various runtime parameters before execution. An example of such a config file can be found in /inres1/aegis_settings.json. AEGIS will search for a config file named "aegis_settings.json" in the current working directory otherwise, a different json file can be passed as a positional argument to the aegis executable. An exhaustive list of options will be provided in documentation. However, the required parameters and some of the more useful optional parameters are listed below:

Required parameters: 
- The path to a DAGMC h5m file in `aegis_params{"DAGMC":}`
- The path to a g-eqdsk file in `equil_params{"eqdsk":}`
- The power at the Scrape-Off-Layer being mapped onto components in `equil_params{"power_sol":}`
- The Scrape-Off-Layer width in `equil_params{"lambda_q":}`
- The radial coordinate of the outer midplane in `equil_params{"r_outrbdry":}`

Optional parameters for the `aegis_params` block:
- `"step_size":` - The straight-line distance a particle moves along fieldline (length of ray-tracing query)
- `"max_steps":` - The max number of steps a particle travels before it is considered `LOST`
- `"launch_pos"` - The launch position of particles on surface elements. Can be either `fixed` (triangle barycentre) or `random` (randomly sampled from element surface) 
   - `"monte_carlo_params":{"number_of_particles_per_facet"}` - If `launch_pos=random` the number of particles sampled can be controlled with this parameter. Deafult = 1
- `"target_surfs":` - The surfaces from which to launch particles in the DAGMC geometry. Values should be provided as a comma seperated array like `[1,2,3,4,5,6]` corresponding to surface IDs in the geometry. If left empty, particles are launched from every surface. This is useful for defining target surfaces and leaving the rest as "shadowing" surfaces 
- `"coordinate_system":` - Specify the tracking coordinate system (currently only cartesian tracking is available)
- `"execution_type":` - Specify whether the simulation should run in parallel or serial. Currently there are 4 available modes `serial`, `mpi` (non-load balanced mpi), `mpi_dynamic` (load balanced mpi). Default is `serial`, either `serial` or `mpi_dynamic` should be preferred.
   - `"dynamic_batch_params":{"batch_size":}` - Specify the number of particles in each batch passed to worker processes. Default = 16.  

Optional parameters for the `equil_params` block:
- `"draw_equil_rz":` - Write out a file containing the 2D RZ grid of equilibrium data read in from the g-eqdsk
- `"draw_equil_xyz":` - Write out a file containing xyz coordinates of the equilibirum grid rotated around in 360 degrees

Optional parameters for the `vtk_params` block:
- `"draw_particle_tracks":` - Draw vtkPolylines for each particle history from birth to termination. A vtk multiblock file will be produced with the extension `.vtm` which can be opened with paraview to view the individual particle tracks grouped by their termination state. 

In order to produce visualisations intended to be used with ParaView, Aegis also makes use of the `VTK library` which should be installed locally. Currently the ability to produce a heatflux distribution across a CAD surface is dependent on the `VTK library`. In the future this may be abstracted out, however the ability to produce particle track plots will likely remain dependent on the `VTK libraries`. The CMakeLists.txt file included will pull the the necessary modules from VTK to produce visualisations as shown:

Some example outputs from an Aegis run with the same magnetic equilibrium and geometry are shown below: 

**Left: Heatflux deposited across target surface with shadowing geometry shown around target**

<p float="left">
  <img src="https://github.com/aurora-multiphysics/aegis/blob/main/gh_images/heatflux_deposited.png" alt="Power Deposited" width="400"/>
  <img src="https://github.com/aurora-multiphysics/aegis/blob/main/gh_images/particle_tracks.png" alt="Particle Tracks" width="400" /> 
</p>

**Right: Individual particle tracks (launched from cells in CAD mesh) coloured by their respective Heatflux with OMP and 2D magnetic equilibrium constructed from G-eqdsk read**


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


In order to produce visualisations intended to be used with ParaView, Aegis also makes use of the `VTK library` which should be installed locally. Currently the ability to produce a heatflux distribution across a CAD surface is dependent on the `VTK library`. In the future this may be abstracted out, however the ability to produce particle track plots will likely remain dependent on the `VTK libraries`. The CMakeLists.txt file included will pull the the necessary modules from VTK to produce visualisations as shown:

Some example outputs from an Aegis run with the same magnetic equilibrium and geometry are shown below: 

**Left: Heatflux deposited across target surface with shadowing geometry shown around target**

<p float="left">
  <img src="https://github.com/aurora-multiphysics/aegis/blob/main/gh_images/heatflux_deposited.png" alt="Power Deposited" width="400"/>
  <img src="https://github.com/aurora-multiphysics/aegis/blob/main/gh_images/particle_tracks.png" alt="Particle Tracks" width="400" /> 
</p>

**Right: Individual particle tracks (launched from cells in CAD mesh) coloured by their respective Heatflux with OMP and 2D magnetic equilibrium constructed from G-eqdsk read**

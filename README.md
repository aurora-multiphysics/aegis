# Aegis 
Aegis is a Monte-Carlo particle tracking tool for calculating power deposition loads (inspired by the SMARDDA FORTRAN code used for the same case) due to particles on surface accurate representations of fusion plasma First Wall Shielding. Makes use of DAGMC/Double-Down as a ray tracing tool that will be used to calculate neutral particle trajectories. Charged particle trajectories will eventually be calculated from field-line tracing and solving of the ODEs which govern these particles trajectories. 

The following dependancies are **required**:
- **Embree v3.6.1** (Intel Embree Ray Tracer)  
- **Double-Down v1.0.0** (A double precision interface to Embree) 
- **MOAB Version 5.2.0** (Mesh Oriented datABase) 
- **DAGMC** (Direct Accelerated Geometry Monte Carlo code)

The following dependancies are **optional**:
- **VTK** (Visualisation ToolKit) - For producing more advaced visualisations within ParaView

# Building the repo
Embree should be built with the following additional flags:

    cmake -DEMBREE_ISPC_SUPPORT=0
Double-Down should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DEMBREE_DIR=\${EMBREE_DIR}
MOAB should be configured with the following additional flags:

    ./src/configure --with-hdf5=\${HDF5_DIR} --enable-optimize --enable-shared --disable-debug
DAGMC should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DDOUBLE_DOWN=on -DDOUBLE_DOWN_DIR=\${DOUBLEDOWN_DIR} -DBUILD_TALLY=ON

And finally this repo should be built with the following flag:

    cmake -DDAGMC_DIR=\${DAGMC_DIR}

In order to produce visualisations to be used with ParaView, Aegis also makes use of the `VTK library` which should be installed locally. The CMakeLists.txt file included will pull the the necessary modules from VTK to produce visualisations as shown:

Some example outputs from an Aegis run with the same magnetic equilibrium and geometry are shown below: 

**Left: Heatflux deposited across the cells in a CAD mesh with OMP in the background**

<p float="left">
  <img src="https://github.com/Waqar-ukaea/aegis/blob/main/gh_images/heatflux_deposited.png" alt="Power Deposited" width="400"/>
  <img src="https://github.com/Waqar-ukaea/aegis/blob/main/gh_images/particle_tracks.png" alt="Particle Tracks" width="390" /> 
</p>

**Right: Individual particle tracks (launched from cells in CAD mesh) coloured by their respective Heatflux with OMP in the background**

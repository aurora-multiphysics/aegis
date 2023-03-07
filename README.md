# dagmc-SMARDDA (name undecided)
A fusion simulation tool to be used to caluclate power loading on first wall components. Makes use of DAGMC/Double-Down as a ray tracing tool that will be used to 
calculate neutral particle trajectories. Charged particle trajectories will mostly likely need to be calculated from field-line tracing and solving of the ODEs which
govern these particles trajectories. The simulations will be performed on 3D CAD surface geometries, with the intention of helping engineers with their design process. 
The following dependancies are required:
- **Embree v3.6.1** (Intel Embree Ray Tracer)  
- **Double-Down v1.0.0** (A double precision interface to Embree) 
- **MOAB Version 5.2.0** (Mesh Oriented datABase) 
- **DAGMC** (Direct Accelerated Geometry Monte Carlo code)
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

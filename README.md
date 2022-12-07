# dagmc_testing (name undecided)
A ray tracing tool that will eventually be used to calculate power loads on 3D CAD geometry, with the intention of use in Nuclear Fusion reactor component design aid. Makes use of Embree/Double-down/DAGMC as ray tracing engine. The following dependancies are required:
- **Embree v3.6.1** (Intel Embree Ray Tracer)  
- **Double-Down v1.0.0** (A double precision interface to Embree) 
- **MOAB Version 5.2.0** (Mesh Oriented datABase) 
- **DAGMC** (Direct Accelerated Geometry Monte Carlo code)
# Building the repo
Embree should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DEMBREE_DIR=\${EMBREE_DIR} -DEMBREE_ISPC_SUPPORT=0
Double-Down should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DEMBREE_DIR=\${EMBREE_DIR}
MOAB should be configured with the following additional flags:

    ./src/configure --with-hdf5=\${HDF5_DIR} --enable-optimize --enable-shared --disable-debug
DAGMC should be built with the following additional flags:

    cmake -DMOAB_DIR=\${MOAB_DIR} -DDOUBLE_DOWN=on -DDOUBLE_DOWN_DIR=\${DOUBLEDOWN_DIR} -DBUILD_TALLY=ON


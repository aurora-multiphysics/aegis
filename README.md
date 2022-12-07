# dagmc_testing (name undecided)
A ray tracing tool that will eventually be used to calculate power loads on 3D CAD geometry. Makes use of Embree/Double-down/DAGMC as ray tracing engine. 
Dependancies required:
- **Embree v3.6.1** (Intel Embree Ray Tracer)
- **Double-Down v1.0.0** (A double precision interface to Embree)
- **MOAB Version 5.2.0** (Mesh Oriented datABase)
- **DAGMC** (Direct Accelerated Geometry Monte Carlo code)
# Building the repo 
    cmake -DDAGMC_DIR=~/[path_to_DAGMC_dir]/lib/cmake/dagmc 
    make
    make test


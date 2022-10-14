# Get base image 
FROM ubuntu:20.04

ARG WORKDIR="opt"
ARG DAGMC_CMAKE_CONFIG="opt/dagmc_bld/DAGMC/lib/cmake/dagmc"

# Install apt dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y \
    apt-get install git \
    apt-get install python3 \
    apt-get install pip \
    apt-get install libeigen3-dev \
    apt-get install libhdf5-dev \
    apt-get install autoconf \
    apt-get install make \
    apt-get install cmake \
    apt-get install libtool \
    apt-get install g++ \ 
    apt-get install gfortran \
    apt-get install libblas-dev \ 
    apt-get install libatlas-base-dev \
    apt-get install liblapack-dev \

# Setting up dagmc_bld directory for DAGMC and MOAB
RUN cd/$WORKDIR && \
    mkdir dagmc_bld && \
    cd dagmc_bld && \
    mkdir -p MOAB\bld && 

# Install MOAB 
RUN cd /$WORKDIR/dagmc_bld/MOAB/ && \
    git clone https://bitbucket.org/fathomteam/moab && \
    cd moab/ && \
    git checkout Version5.1.0 && \
    autoreconf -fi && \
    cd ../ && \
    ln -s moab src && \ 
    cd bld/ && \
    ../src/configure --enable-optimize --enable-shared --disable-debug \
    		     --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
	             --prefix=$HOME/dagmc_bld/MOAB && \
    make && \ 
    make check && \
    make install && \
    echo 'export PATH=$PATH:$HOME/dagmc_bld/HDF5/bin' >> ~/.bashrc && \
    echo 'export PATH=$PATH:$HOME/dagmc_bld/MOAB/bin' >> ~/.bashrc && \
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/HDF5/lib' >> ~/.bashrc && \
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/MOAB/lib' >> ~/.bashrc 

#Install DAGMC
RUN cd /$WORKDIR/dagmc_bld/ && \
    mkdir DAGMC && \
    cd DAGMC && \
    git clone https://github.com/svalinn/DAGMC && \
    cd DAGMC && \
    git checkout develop && \
    git submodule update --init && \
    cd $WORKDIR/dagmc_bld/DAGMC/ && \
    ls -s DAGMC src && \
    mkdir bld && \
    cd bld/ && \ 
    INSTALL_PATH=$WORKDIR/dagmc_bld/DAGMC && \
    cmake ../src -DMOAB_DIR=$WORKDIR/dagmc_bld/MOAB \
                 -DBUILD_TALLY=ON \
                 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH && \
    make && \
    make install && \
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/dagmc_bld/dagmc/lib' >> ~/.bashrc 
   
    
    

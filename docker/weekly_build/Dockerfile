# Get base image 
FROM ubuntu:20.04

ARG WORKDIR="opt"


# Install apt dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y \
    docker \
    git \
    python3 \
    pip \
    libeigen3-dev \
    libhdf5-dev \
    autoconf \
    make \
    cmake \
    libtool \
    g++ \ 
    gcc \
    gfortran \
    libblas-dev \ 
    libatlas-base-dev \
    liblapack-dev \
    build-essential libboost-system-dev \
    build-essential libboost-thread-dev \
    build-essential libboost-program-options-dev \
    build-essential libboost-test-dev \
    libboost-log-dev \ 
    libblas-dev \
    liblapack-dev \
    libhdf5-dev \
    libeigen3-dev \
    libglfw3-dev  \
    libtbb-dev \
    pkg-config \
    rsync \
    bison \
    flex

# Setting up dagmc_bld directory for DAGMC and MOAB
RUN cd /$WORKDIR && \
    mkdir dagmc_bld && \
    cd dagmc_bld && \
    mkdir -p MOAB/bld

# Install MOAB 
RUN cd /$WORKDIR/dagmc_bld/MOAB/ && \
    git clone https://bitbucket.org/fathomteam/moab && \
    cd moab/ && \
    git checkout Version5.2.0 && \
    autoreconf -fi && \
    cd ../ && \
    ln -s moab src && \ 
    cd bld/ && \
    ../src/configure --enable-optimize --enable-shared --disable-debug \
    		     --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
	             --prefix=/$WORKDIR/dagmc_bld/MOAB && \
    make && \ 
    make check && \
    make install && \
    echo 'export PATH=$PATH:$/WORKDIR/dagmc_bld/HDF5/bin' >> ~/.bashrc && \
    echo 'export PATH=$PATH:$/WORKDIR/dagmc_bld/MOAB/bin' >> ~/.bashrc && \
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$/WORKDIR/dagmc_bld/HDF5/lib' >> ~/.bashrc && \
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$/WORKDIR/dagmc_bld/MOAB/lib' >> ~/.bashrc 

# Install Embree 
RUN cd /$WORKDIR/dagmc_bld/ && \
    mkdir EMBREE && \
    cd EMBREE && \
    git clone https://github.com/embree/embree.git && \
    cd embree && \
    git checkout v3.6.1 && \
    cmake -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_C_COMPILER=/usr/bin/gcc -DEMBREE_ISPC_SUPPORT=0 && \
    make && \
    make install
    
# Install Double-Down
RUN cd /$WORKDIR/dagmc_bld/ && \
    mkdir DOUBLE-DOWN && \
    cd DOUBLE-DOWN && \
    git clone https://github.com/pshriwise/double-down && \
    cd double-down && \
    cmake -DMOAB_DIR=/$WORKDIR/dagmc_bld/MOAB -DEMBBREE_DIR=/$WORKDIR/dagmc_bld/EMBREE/embree && \
    make && \
    make install && \
    make test
    

# Install DAGMC
RUN cd /$WORKDIR/dagmc_bld/ && \
    mkdir DAGMC && \
    cd DAGMC && \
    git clone https://github.com/svalinn/DAGMC && \
    cd DAGMC && \
    git checkout develop && \
    git submodule update --init && \
    cd /$WORKDIR/dagmc_bld/DAGMC/ && \
    ln -s DAGMC src && \
    mkdir bld && \
    cd bld/ && \ 
    INSTALL_PATH=/$WORKDIR/dagmc_bld/DAGMC && \
    cmake ../src -DMOAB_DIR=/$WORKDIR/dagmc_bld/MOAB \
    		 -DDOUBLE_DOWN=ON \
		 -DDOUBLE_DOWN_DIR=/$WORKDIR/dagmc_bld/DOUBLE-DOWN/double-down \
                 -DBUILD_TALLY=ON \
                 -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH && \
    make && \
    make install && \
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$/WORKDIR/dagmc_bld/dagmc/lib' >> ~/.bashrc && \
    cd ../lib/cmake/dagmc && \
    DAGMC_CMAKE_CONFIG=$(pwd)
   
# Install VTK 
RUN cd /$WORKDIR/ && \
    mkdir VTK && \
    git clone https://github.com/Kitware/VTK.git && \
    cd VTK && \
    mkdir VTK-build && \
    cd VTK-build && \
    cmake ../ && \ 
    make && \
    make install
    

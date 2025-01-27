#!/bin/bash
#set -ue

export WORKDIR=$PWD
mkdir aegis-deps
export compile_cores=24

#function load_modules() {
#  module purge
#  module load hdf5/openmpi/gcc/9.3/openmpi-4.0.4/1.12.0
#  module load eigen
#}

install_apt_deps() {

  # Update the package list
  apt-get update

  # Install required packages
  apt-get install -y \
    wget \
    python3 \
    pip \
    libeigen3-dev \
    libhdf5-serial-dev \
    autoconf \
    make \
    cmake \
    libtool \
    libopenmpi-dev \
    g++ \
    gcc \
    gfortran \
    libblas-dev \
    libatlas-base-dev \
    liblapack-dev \
    build-essential \
    libboost-system-dev \
    libboost-thread-dev \
    libboost-program-options-dev \
    libboost-test-dev \
    libboost-log-dev \
    libtbb-dev \
    libx11-dev
    pkg-config
}

function build_moab(){
  cd $WORKDIR/aegis-deps
  mkdir MOAB
  cd MOAB
  git clone https://bitbucket.org/fathomteam/moab
  cd moab
#  git checkout Version5.2.0
  autoreconf -fi
  cd ../
  ln -s moab src
  mkdir bld
  cd bld 
  ../src/configure --enable-optimize --enable-shared --disable-debug \
    		           --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial
  make -j"$compile_cores"
  make check -j"$compile_cores"
  sudo make install
}

function build_embree(){
  cd $WORKDIR/aegis-deps
  mkdir EMBREE 
  cd EMBREE
  wget https://github.com/embree/embree/releases/download/v4.3.3/embree-4.3.3.x86_64.linux.tar.gz
  tar xzf embree-4.3.3.x86_64.linux.tar.gz
  rm embree-4.3.3.x86_64.linux.tar.gz
  source embree-vars.sh

}

function build_double_down(){
  cd $WORKDIR/aegis-deps
  mkdir DOUBLE-DOWN
  cd DOUBLE-DOWN
  git clone https://github.com/pshriwise/double-down
  cd double-down
  cmake -DMOAB_DIR=$WORKDIR/aegis-deps/MOAB -DEMBREE_DIR=$WORKDIR/aegis-deps/EMBREE/ -DCMAKE_INSTALL_PREFIX=$WORKDIR/aegis-deps/DOUBLE-DOWN
  make -j"$compile_cores"
  make install
  make test
}

function build_dagmc(){
  cd $WORKDIR/aegis-deps
  mkdir DAGMC
  cd DAGMC
  git clone https://github.com/svalinn/DAGMC
  cd DAGMC
  git checkout develop
  git submodule update --init --recursive
  cd $WORKDIR/aegis-deps/DAGMC
  ln -s DAGMC src
  mkdir bld
  cd bld
  INSTALL_PATH=$WORKDIR/aegis-deps/DAGMC
  cmake ../src -DMOAB_DIR=$WORKDIR/aegis-deps/MOAB \
    		       -DDOUBLE_DOWN=ON \
               -DDOUBLE_DOWN_DIR=$WORKDIR/aegis-deps/DOUBLE-DOWN/double-down \
               -DBUILD_TALLY=ON \
               -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH
  make -j"$compile_cores"
  make install 
  # make test
}

function build_vtk(){
  cd $WORKDIR/aegis-deps
  mkdir VTK
  git clone https://github.com/Kitware/VTK.git
  cd VTK
  git checkout v9.3.1
  mkdir VTK-build
  cd VTK-build
  cmake ../ -DCMAKE_INSTALL_PREFIX=$WORKDIR/aegis-deps/VTK
  make -j"$compile_cores"
  make install
}

function build_aegis(){
  cd $WORKDIR
  git clone https://github.com/aurora-multiphysics/aegis.git
  cd aegis
  git submodule update --init --recursive 
  mkdir bld
  cd bld
  cmake ../ -DDAGMC_DIR=$WORKDIR/aegis-deps/DAGMC/ -DCMAKE_PREFIX_PATH=$WORKDIR/aegis-deps/VTK
  make 
}

#load_modules
install_apt_deps
build_moab
build_embree
build_double_down
build_dagmc
build_vtk
build_aegis

# Get base image 
FROM waqar1711/aegis:latest

ARG WORKDIR="opt"

# Install aegis  
RUN cd /$WORKDIR && \
    mkdir temp && \
    cd temp && \
    git clone https://github.com/aurora-multiphysics/aegis.git && \
    cd aegis && \
    mkdir bld && \
    cd bld && \
    git submodule update --init --recursive && \
    cmake ../ -DDAGMC_DIR=/$WORKDIR/dagmc_bld/DAGMC/lib/cmake/dagmc/ && \
    make && \
    make test 

    

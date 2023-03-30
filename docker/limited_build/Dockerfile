# Get base image 
FROM waqar1711/aegis:latest

ARG WORKDIR="opt"

# Install aegis  
RUN cd /$WORKDIR && \
    mkdir temp && \
    cd temp && \
    git clone https://github.com/Waqar-ukaea/aegis.git && \
    cd aegis && \
    mkdir bld && \
    cd bld && \
    cmake ../ -DDAGMC_DIR=/$WORKDIR/dagmc_bld/DAGMC/lib/cmake/dagmc/ && \
    make && \
    make test


    

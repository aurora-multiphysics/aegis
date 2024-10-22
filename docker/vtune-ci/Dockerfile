# Get base image 
FROM waqar1711/aegis:latest

ARG WORKDIR="opt"

# Install aegis  
RUN cd /$WORKDIR && \
    git clone --recurse https://github.com/aurora-multiphysics/aegis.git && \
    cd aegis && \
    mkdir bld && \
    cd bld && \
    cmake ../ -DDAGMC_DIR=/$WORKDIR/dagmc_bld/DAGMC/lib/cmake/dagmc/ && \
    make && \
    make test 

# Install Intel VTune
    RUN apt update && \
    apt install -y gpg-agent wget && \
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
    apt update && \
    apt install -y intel-oneapi-vtune

# Run Intel VTune hotspots analysis 
    RUN cd $WORKDIR && \
    mkdir artifacts && \
    cd artifacts && \
    source /opt/intel/oneapi/setvars.sh && \
    vtune -collect hotspots -strategy=:trace:trace --app-working-dir=$WORKDIR/aegis/inres1 $WORKDIR/aegis/bin/aegis && \
    vtune -report summary -format=html > hotspots_summary.html
    

# Get base image 
FROM waqar1711/aegis:latest

ARG WORKDIR="opt"

RUN apt-get update && apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates gpg-agent software-properties-common && \
    rm -rf /var/lib/apt/lists/*
# Repository to install Intel(R) oneAPI Libraries
RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB | gpg --dearmor | tee /usr/share/keyrings/intel-oneapi-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/intel-oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main " > /etc/apt/sources.list.d/oneAPI.list

RUN mkdir -p /opt/workdir
RUN apt-get update && apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates build-essential pkg-config gnupg libarchive13 openssh-server openssh-client wget net-tools git intel-oneapi-vtune  && \
    rm -rf /var/lib/apt/lists/*

# VTune environment variables
ENV LANG=C.UTF-8
ENV ONEAPI_ROOT='/opt/intel/oneapi'
ENV PATH='/opt/intel/oneapi/vtune/2024.2/bin64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
ENV PKG_CONFIG_PATH='/opt/intel/oneapi/vtune/2024.2/include/pkgconfig/lib64:'
ENV SETVARS_COMPLETED='0'
ENV VTUNE_PROFILER_2024_DIR='/opt/intel/oneapi/vtune/2024.2'
ENV VTUNE_PROFILER_DIR='/opt/intel/oneapi/vtune/2024.2'

# Install aegis  
RUN cd /$WORKDIR && \
    git clone --recurse https://github.com/aurora-multiphysics/aegis.git && \
    cd aegis && \
    mkdir bld && \
    cd bld && \
    cmake ../ -DDAGMC_DIR=/$WORKDIR/dagmc_bld/DAGMC/lib/cmake/dagmc/ && \
    make && \
    make test 

# Collect VTune results
RUN cd ${WORKDIR} && \
    . /opt/intel/oneapi/setvars.sh && \
    vtune -collect hotspots -strategy=:trace:trace --app-working-dir=/opt/aegis/inres1 /opt/aegis/bin/aegis && \
    vtune -report summary -format=html > hotspots_summary.html && \
    tar -cvf hotspots_results.tar r000hs/ 

# Get base image 
FROM waqar1711/smardda-dagmc:v1

ARG WORKDIR="opt"

# Install dagmc_testing  
RUN cd /$WORKDIR && \
    mkdir temp && \
    cd temp && \
    git clone https://github.com/Waqar-ukaea/aegis.git && \
    cd aegis && \
    cmake -DDAGMC_DIR=/$WORKDIR/dagmc_bld/DAGMC/lib/cmake/dagmc/ && \
    make && \
    make test


    

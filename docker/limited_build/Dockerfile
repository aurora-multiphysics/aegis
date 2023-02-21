# Get base image 
FROM waqar1711/smardda-dagmc:v1

ARG WORKDIR="opt"

# Install dagmc_testing  
RUN cd /$WORKDIR && \
    git clone https://github.com/Waqar-ukaea/dagmc_testing.git && \
    cd dagmc_testing && \
    cmake -DDAGMC_DIR=/$WORKDIR/dagmc_bld/DAGMC/lib/cmake/dagmc/ && \
    make && \
    make test


    

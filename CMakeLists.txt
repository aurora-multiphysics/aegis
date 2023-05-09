cmake_minimum_required(VERSION 3.8)

# no space in these names! used as shared object name and python module name
set(PROJECT_NAME aegis) 
set(PROJECT_LIB_NAME ${PROJECT_NAME}-lib)

project(${PROJECT_NAME} VERSION 0.0.0 LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 17)
#set (CMAKE_CXX_CLANG_TIDY 
#     clang-tidy;
#     -checks=*)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREAD ON)
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
set(Boost_PACKAGES system log log_setup)
set(BOOST_REQUESTED_VERSION 1.71)

find_package(Boost REQUIRED COMPONENTS ${Boost_PACKAGES})

find_package(OpenMP REQUIRED)


#include(${CMAKE_SOURCE_DIR}/cmake/SetBuildType.cmake)

# Set debugging symbols
set(CMAKE_BUILD_TYPE Debug)


# Import DAGMC
message(STATUS "Looking for DAGMC...")
set(DAGMC_DIR="~/dagmc_bld/DAGMC/" CACHE PATH "Path to the DAGMC build or install prefix.")

if(DAGMC_DIR) 
   find_package(DAGMC REQUIRED NAMES DAGMC HINTS "${DAGMC_DIR}" 
		"${DAGMC_DIR}/lib/cmake/dagmc" NO_DEFAULT PATH)
else()
find_package(DAGMC REQUIRED lib/cmake)
endif()
message(STATUS "Found DAGMC config in: ${DAGMC_DIR}")

#target_link_libraries(${PROJECT_NAME} dagmc-shared uwuw-shared )
#target_compile_definitions( PRIVATE DAGMC)
#target_link_libraries(libopenmc dagmc-shared uwuw-shared)


message(STATUS "DAGMC_INCLUDE_DIRS is set to: ${DAGMC_INCLUDE_DIRS}")
message(STATUS "DAGMC_LIBRARIES is set to: ${DAGMC_LIBRARIES}")

find_library(libopenmc.so ~/openmc/openmc/openmc/lib/)


enable_testing()

## put all targets in bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
## put all libraries in lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

############################################################
# Set the build type environment variable
############################################################
if(BUILD_TYPE STREQUAL "TEST")
  set(BUILD_TYPE_COMPILE_FLAGS "-g;-O0;--coverage")
  set(TEST_LIBRARIES "gcov")
elseif(BUILD_TYPE STREQUAL "DEBUG")
  set(BUILD_TYPE_COMPILE_FLAGS "-g -O0")
  set(TEST_LIBRARIES "")
elseif(BUILD_TYPE STREQUAL "RELEASE")
  set(BUILD_TYPE_COMPILE_FLAGS "-O2")
  set(TEST_LIBRARIES "")
endif()

############################################################
# Add directories with source files
############################################################
add_subdirectory(src)
add_subdirectory(test)

############################################################
# all configurable options
############################################################
#find_package(Doxygen)
#option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})
#option(BUILD_PYTHON "build python interface module (requires python-dev or python3-dev)" ON)

############################################################
# Create python wrapper by pybind11
############################################################
if(BUILD_PYTHON)
    add_subdirectory(python)
endif()

############################################################
# Create doxygen documentation `make doc`
############################################################
#if(BUILD_DOCUMENTATION)
#    if(NOT DOXYGEN_FOUND)
#        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
#    endif()
#
#    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
#    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
#    configure_file(${doxyfile_in} ${doxyfile} @ONLY)
#
#    add_custom_target(doc
#        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        COMMENT "Generating API documentation with Doxygen"
#        VERBATIM)
#endif()

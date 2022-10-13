cmake_minimum_required(VERSION 3.10)

# set the project name
project(
        dagmc_call 
        VERSION 1.0
        LANGUAGES CXX)

find_package(DAGMC REQUIRED lib/cmake)


add_executable(${PROJECT_NAME} dagmc_call.cpp)

target_link_libraries(dagmc_call dagmc-shared uwuw-shared )

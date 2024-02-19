file(GLOB_RECURSE lib_files "*.h" "*.hpp" "*.cpp")

message(STATUS "lib files:")
message(STATUS ${lib_files})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


# make a shared object library to link against for testing
add_library(${PROJECT_LIB_NAME} SHARED ${lib_files})
target_compile_options(${PROJECT_LIB_NAME} PUBLIC ${BUILD_TYPE_COMPILE_FLAGS})
set_target_properties(${PROJECT_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(${PROJECT_LIB_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_include_directories(${PROJECT_LIB_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC ${Boost_LIBRARIES} ${TEST_LIBRARIES})
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC dagmc-shared uwuw-shared )
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC MPI::MPI_CXX)
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC ${VTK_LIBRARIES})
target_link_libraries(${PROJECT_LIB_NAME} PUBLIC nlohmann_json::nlohmann_json)

# vtk_module_autoinit is needed
#vtk_module_autoinit(
#  TARGETS ${PROJECT_NAME_LIB}
#  MODULES ${VTK_LIBRARIES}
#)


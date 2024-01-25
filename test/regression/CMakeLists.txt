file(GLOB_RECURSE src_regression "*.cpp")


include_directories(${Boost_INCLUDE_DIRS})
add_executable(regression_tests ${src_regression})
add_compile_options(regression_tests ${BUILD_TYPE_COMPILER_FLAGS})
target_include_directories(regression_tests PUBLIC ${GTEST_INCLUDE_DIRS})
target_include_directories(regression_tests PUBLIC ${GMOCK_INCLUDE_DIRS})
target_include_directories(regression_tests PUBLIC ${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}_lib)
target_include_directories(regression_tests PUBLIC ${DAGMC_INCLUDE_DIRS})
target_include_directories(regression_tests PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(regression_tests PUBLIC ${VTK_LIBRARIES})

target_link_libraries(regression_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} pthread)
target_link_libraries(regression_tests ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/lib${PROJECT_NAME}.so)
target_link_libraries(regression_tests ${Boost_LIBRARIES})
target_link_libraries(regression_tests ${DAGMC_LIBRARY} dagmc-shared uwuw-shared )
target_link_libraries(regression_tests OpenMP::OpenMP_CXX)
target_link_libraries(regression_tests ${VTK_LIBRARIES})
#install(TARGETS regression_tests DESTINATION bin)

include(GoogleTest)
gtest_discover_tests(regression_tests WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test/data")

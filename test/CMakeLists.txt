Find_package(Threads REQUIRED)
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lpthread")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREAD ON)
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
set(Boost_PACKAGES system log log_setup)
set(BOOST_REQUESTED_VERSION 1.71)

find_package(Boost REQUIRED COMPONENTS ${Boost_PACKAGES})
find_package(MPI REQUIRED)

include(ExternalProject)

ExternalProject_Add(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  GIT_TAG main
  )

ExternalProject_Get_Property(googletest source_dir)
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)
set(GMOCK_INCLUDE_DIRS ${source_dir}/googlemock/include)

ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBRARY_PATH  ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
set(GMOCK_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a)
set(GTEST_MAIN_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
set(GMOCK_MAIN_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock_main.a)
set(GTEST_LIBRARY gtest)
set(GMOCK_LIBRARY gmock)
set(GTEST_MAIN_LIBRARY gtest_main)
set(GMOCK_MAIN_LIBRARY gmock_main)
add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
add_library(${GMOCK_LIBRARY} UNKNOWN IMPORTED)
add_library(${GTEST_MAIN_LIBRARY} UNKNOWN IMPORTED)
add_library(${GMOCK_MAIN_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
  IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(${GMOCK_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GMOCK_LIBRARY_PATH}
  IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(${GTEST_MAIN_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GTEST_MAIN_LIBRARY_PATH}
  IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(${GMOCK_MAIN_LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${GMOCK_MAIN_LIBRARY_PATH}
  IMPORTED_LINK_INTERFACE_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(${GTEST_LIBRARY} googletest)
add_dependencies(${GMOCK_LIBRARY} googletest)
add_dependencies(${GTEST_MAIN_LIBRARY} googletest)
add_dependencies(${GMOCK_MAIN_LIBRARY} googletest)


# Set up to generate coverage reports
find_package(PythonInterp REQUIRED)
find_program(GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/scripts/test)


IF(BUILD_TYPE STREQUAL "TEST")
  add_custom_target(gcovr
    # Create folder
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/gcovr

    # Running gcovr
    COMMAND ${Python_EXECUTABLE} ${GCOVR_PATH} --html --html-details -s
            -r ${PROJECT_SOURCE_DIR} --object-directory=${PROJECT_BINARY_DIR}
            -o gcovr/index.html
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Running gcovr to produce HTML code coverage report."
  )
  # Show info where to find the report
  add_custom_command(TARGET gcovr POST_BUILD
    COMMAND ;
    COMMENT "Open ./gcovr/index.html in your browser to view the full coverage report."
  )
ENDIF()

add_subdirectory(unit)
add_subdirectory(integration)
add_subdirectory(regression)


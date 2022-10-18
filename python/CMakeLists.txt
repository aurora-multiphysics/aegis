######################################
# pybind11 wrap
######################################

find_package(PythonInterp REQUIRED)
message(STATUS "Found Python: ${PYTHON_EXECUTABLE} (found version \"${PYTHON_VERSION_STRING}\")")

#first of all detect system installation, then check local folder,
find_package(pybind11 QUIET)
if(pybind11_FOUND)
    message("system wide pybind11 installation found\n")
    # systemwide installation by package manager with header in `/usr/include/`
else()
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include/")
        message("===============================================================\n"
                "pybind11 not found on system and ${CMAKE_CURRENT_SOURCE_DIR}, run git submodule add\n"
                "===============================================================\n")
        execute_process(COMMAND git submodule add --force https://github.com/pybind/pybind11.git
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    ERROR_VARIABLE error_var
                    RESULT_VARIABLE ret_var)
        message("rv='${ret_var}'\n error msg = '${error_var}'")
    else()
        message("pybind11 is found in the folder: ${CMAKE_CURRENT_SOURCE_DIR}\n")
    endif()

    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/pybind11")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include")
    SET(pybind11_FOUND TRUE)
endif()

if(pybind11_FOUND)
    # this will generate multiple versions of py module
    #set(PYBIND11_PYTHON_VERSION 2.7 3.6)
    # use the defined ${PROJECT_LIB_NAME} as python module name/target name
    set(PY_MODULE_NAME ${PROJECT_NAME}_wrap)
    pybind11_add_module(${PY_MODULE_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}_pybind11.cpp)

    target_include_directories(${PY_MODULE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)

    #target_link_directories()
    link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_libraries(${PY_MODULE_NAME} PRIVATE ${PROJECT_LIB_NAME})
    # the output wrap module are put into ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
endif()

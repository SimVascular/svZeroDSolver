# This file is part of svZeroDSolver licensed under the MIT License
# See the LICENSE.md file for license information

set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(svZeroDSolver)

include(FetchContent)

# Enable code coverage
# -----------------------------------------------------------------------------
set(ENABLE_COVERAGE OFF CACHE BOOL "Enable code coverage")
# coverage
if(ENABLE_COVERAGE)
  # set compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -coverage -g")

  # find required tools
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml REQUIRED)

  # add coverage target
  add_custom_target(coverage
    # gather data
    COMMAND ${LCOV} --directory . --capture --output-file coverage.info --ignore-errors gcov --exclude '/usr/include/*' --exclude '/usr/lib/*' --exclude '*/_deps/*'
    # generate report
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# -----------------------------------------------------------------------------
# Set the location to store the binaries and libraries created by this project.
# -----------------------------------------------------------------------------
#
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Need to add the -fPIC flag to build the interface shared library.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

# -----------------------------------------------------------------------------
# Fetch Eigen
# -----------------------------------------------------------------------------
# Eigen is a header-only C++ template library for linear algebra.
#
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)

set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
set( OFF)

FetchContent_MakeAvailable(Eigen)

# -----------------------------------------------------------------------------
# Fetch pybind11
# -----------------------------------------------------------------------------
# pybind11 is a library used to create Python bindings of existing C++ code.
#
set(PYBIND11_PYTHON_VERSION 3.8 CACHE STRING "")
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        master
)

FetchContent_GetProperties(pybind11)

if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

# -----------------------------------------------------------------------------
# Fetch nlohmann/json
# -----------------------------------------------------------------------------
# This is a C++ header-only library for reading JSON files.
#
FetchContent_Declare(
    json
    URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)

FetchContent_MakeAvailable(json)

# -----------------------------------------------------------------------------
# Fetch pybind11_json
# -----------------------------------------------------------------------------
# pybind11_json is an nlohmann::json to pybind11 bridge that converts 
# between nlohmann::json and py::object objects. 
#
FetchContent_Declare(
    pybind11_json
    GIT_REPOSITORY https://github.com/pybind/pybind11_json.git
    GIT_TAG        master
)

FetchContent_MakeAvailable(pybind11_json)

# -----------------------------------------------------------------------------
# Set executables for each application.
# -----------------------------------------------------------------------------
# There are two executables
#   1) svzerodsolver
#   2) svzerodcalibrator 
#
add_executable(svzerodsolver applications/svzerodsolver.cpp 
  $<TARGET_OBJECTS:svzero_algebra_library> 
  $<TARGET_OBJECTS:svzero_model_library> 
  $<TARGET_OBJECTS:svzero_solve_library> 
)

add_executable(svzerodcalibrator applications/svzerodcalibrator.cpp
  $<TARGET_OBJECTS:svzero_algebra_library> 
  $<TARGET_OBJECTS:svzero_optimize_library> 
  $<TARGET_OBJECTS:svzero_model_library> 
)

# -----------------------------------------------------------------------------
# Setup building of the pysvzerod Python extension module.
# -----------------------------------------------------------------------------
# Replace EXCLUDE_FROM_ALL with SHARED to test building the Python 
# shared library.
pybind11_add_module(pysvzerod EXCLUDE_FROM_ALL applications/pysvzerod.cpp)

# -----------------------------------------------------------------------------
# Add source sub-directories.
# -----------------------------------------------------------------------------
add_subdirectory("src/algebra")
add_subdirectory("src/interface")
add_subdirectory("src/model")
add_subdirectory("src/optimize")
add_subdirectory("src/solve")

# -----------------------------------------------------------------------------
# Set header file include directories.
# -----------------------------------------------------------------------------
#
target_include_directories(svzerodsolver PUBLIC
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model
  ${CMAKE_SOURCE_DIR}/src/solve
)

target_include_directories(svzerodcalibrator PUBLIC 
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model 
  ${CMAKE_SOURCE_DIR}/src/optimize
  ${CMAKE_SOURCE_DIR}/src/solve 
)

target_include_directories(pysvzerod PUBLIC 
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model 
  ${CMAKE_SOURCE_DIR}/src/optimize
  ${CMAKE_SOURCE_DIR}/src/solve 
)

# -----------------------------------------------------------------------------
# Set the libraries each application needs to link to.
# -----------------------------------------------------------------------------
# These are needed only for the header-only applications.
#
target_link_libraries(svzerodsolver PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodsolver PRIVATE nlohmann_json::nlohmann_json)

target_link_libraries(svzerodcalibrator PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodcalibrator PRIVATE nlohmann_json::nlohmann_json)

target_link_libraries(pysvzerod PRIVATE Eigen3::Eigen)
target_link_libraries(pysvzerod PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(pysvzerod PRIVATE pybind11_json)
target_link_libraries(pysvzerod PRIVATE svzero_algebra_library)
target_link_libraries(pysvzerod PRIVATE svzero_model_library)
target_link_libraries(pysvzerod PRIVATE svzero_optimize_library)
target_link_libraries(pysvzerod PRIVATE svzero_solve_library)


# Create distribution
set(ENABLE_DISTRIBUTION OFF CACHE BOOL "Enable installer build")
add_subdirectory("distribution")


# -----------------------------------------------------------------------------
# Enforce and check code format
# -----------------------------------------------------------------------------
# set paths of files to format
set(SDIR ${PROJECT_SOURCE_DIR}/src/**/)

# format the code
add_custom_target(codeformat
  COMMAND find ${SDIR}/*.h ${SDIR}/*.cpp | xargs clang-format -style=file:${PROJECT_SOURCE_DIR}/.clang-format -i)

# check code format
add_custom_target(codecheck
  COMMAND find ${SDIR}/*.h ${SDIR}/*.cpp | xargs clang-format -style=file:${PROJECT_SOURCE_DIR}/.clang-format --dry-run --Werror)

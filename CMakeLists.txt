cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(svZeroDSolver)
include(FetchContent)

# -----------------------------------------------------------------------------
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
    COMMAND ${LCOV} --directory . --capture --output-file coverage.info --ignore-errors inconsistent --ignore-errors gcov
    # exclude externals
    COMMAND ${LCOV} --remove coverage.info -o coverage.info --ignore-errors unused '/usr/*' '/opt/*' '_deps/*' '/Library' --ignore-errors inconsistent
    # generate report
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info --ignore-errors inconsistent --ignore-errors unmapped
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# -----------------------------------------------------------------------------
# Fetch Eigen
# -----------------------------------------------------------------------------
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
FetchContent_Declare(
    json
    URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)
FetchContent_MakeAvailable(json)

# -----------------------------------------------------------------------------
# Fetch pybind11_json
# -----------------------------------------------------------------------------
FetchContent_Declare(
    pybind11_json
    GIT_REPOSITORY https://github.com/pybind/pybind11_json.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(pybind11_json)

# -----------------------------------------------------------------------------
# Set executable and python library
# -----------------------------------------------------------------------------
include_directories(src)
add_executable(svzerodsolver applications/svzerodsolver.cpp)
add_executable(svzerodcalibrator applications/svzerodcalibrator.cpp)
pybind11_add_module(svzerodplus EXCLUDE_FROM_ALL applications/svzerodplus.cpp)
add_subdirectory("src/interface")

# -----------------------------------------------------------------------------
# Link libraries
# -----------------------------------------------------------------------------
target_link_libraries(svzerodsolver PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodsolver PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(svzerodcalibrator PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodcalibrator PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(svzerodplus PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodplus PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(svzerodplus PRIVATE pybind11_json)

cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(svZeroDSolver)
include(FetchContent)

# -----------------------------------------------------------------------------
# Fetch simdjson
# -----------------------------------------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
FetchContent_Declare(
  simdjson
  GIT_REPOSITORY https://github.com/simdjson/simdjson.git
  GIT_TAG  master
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(simdjson)

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
# Set executable and python library
# -----------------------------------------------------------------------------
include_directories(src)
add_executable(solver applications/solver.cpp)
add_executable(calibrator applications/calibrator.cpp)
pybind11_add_module(libsvzerodsolver src/python.cpp)
add_subdirectory("src/interface")

# -----------------------------------------------------------------------------
# Link libraries
# -----------------------------------------------------------------------------
target_link_libraries(solver PRIVATE Eigen3::Eigen)
target_link_libraries(solver PRIVATE simdjson)
target_link_libraries(calibrator PRIVATE Eigen3::Eigen)
target_link_libraries(calibrator PRIVATE simdjson)
target_link_libraries(calibrator PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(libsvzerodsolver PRIVATE Eigen3::Eigen)
target_link_libraries(libsvzerodsolver PRIVATE simdjson)

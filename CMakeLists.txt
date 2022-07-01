cmake_minimum_required(VERSION 3.23)
project(svZeroDSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(svzerodsolver src/main.cpp src/external/simdjson/singleheader/simdjson.cpp src/external/eigen/Eigen/Dense src/external/eigen/Eigen/Sparse src/external/eigen/Eigen/SparseLU)

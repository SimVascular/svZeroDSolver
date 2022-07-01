# Dependencies for macOS:
# brew install eigen

cmake_minimum_required(VERSION 3.23)
project(svZeroDSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
add_executable(svzerodsolver src/main.cpp src/io/simdjson.cpp)

target_link_libraries(svzerodsolver Eigen3::Eigen)

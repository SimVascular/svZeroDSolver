cmake_minimum_required(VERSION 3.23)
project(svZeroDSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(simdjson CONFIG REQUIRED)

add_executable(svzerodsolver src/main.cpp)

target_link_libraries(svzerodsolver Eigen3::Eigen)
target_link_libraries(svzerodsolver simdjson::simdjson)

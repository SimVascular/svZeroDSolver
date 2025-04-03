# This file is part of svZeroDSolver licensed under the MIT License
# See the LICENSE.md file for license information

# Build the optimize object files.

set(lib svzero_optimize_library)

set(CXXSRCS LevenbergMarquardtOptimizer.cpp calibrate.cpp)

set(HDRS LevenbergMarquardtOptimizer.h calibrate.h)

add_library(${lib} OBJECT ${CXXSRCS} )

target_include_directories( ${lib} PUBLIC 
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model 
  ${CMAKE_SOURCE_DIR}/src/solve 
)

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )


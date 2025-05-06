# This file is part of svZeroDSolver licensed under the MIT License
# See the LICENSE.md file for license information

# Build the algebra object files.

set(lib svzero_algebra_library)

set(CXXSRCS Integrator.cpp SparseSystem.cpp State.cpp )

set(HDRS Integrator.h SparseSystem.h State.h )

add_library(${lib} OBJECT ${CXXSRCS} )

target_include_directories(${lib} PUBLIC
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model
  ${CMAKE_SOURCE_DIR}/src/solve
)

target_link_libraries( ${lib} Eigen3::Eigen )


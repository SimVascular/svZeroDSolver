# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

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


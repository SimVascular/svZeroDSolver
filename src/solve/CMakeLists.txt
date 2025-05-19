# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

# Build the solver object files.

set(lib svzero_solve_library)

set(CXXSRCS 
  csv_writer.cpp 
  SimulationParameters.cpp 
  Solver.cpp
)

set(HDRS 
  csv_writer.h 
  debug.h 
  SimulationParameters.h 
  Solver.h 
)

add_library(${lib} OBJECT ${CXXSRCS} )

target_include_directories(${lib} PUBLIC
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model
  ${CMAKE_SOURCE_DIR}/src/solve
)

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )


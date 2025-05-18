# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

# Set the library name.
set(lib svzero_interface)

set(CXXSRCS interface.cpp )
set(HDRS interface.h )

add_library(${lib} SHARED ${CXXSRCS}) 

target_include_directories(${lib} PUBLIC 
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model 
  ${CMAKE_SOURCE_DIR}/src/solve 
)

# Set the libraries and object files that need to be linked 
# to when building the library.
#
target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json ) 
target_link_libraries( ${lib} svzero_algebra_library)
target_link_libraries( ${lib} svzero_model_library)
target_link_libraries( ${lib} svzero_solve_library)

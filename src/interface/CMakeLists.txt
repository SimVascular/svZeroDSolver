# This file is part of svZeroDSolver licensed under Stanford University, The Regents of the University of 
#                                                   California, and others.
# 
# See the LICENSE.md file for license information

# This CMake file is used to build the shared library 'svzero_interface'
# used to interface to the svZeroDSolver from other programs.

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

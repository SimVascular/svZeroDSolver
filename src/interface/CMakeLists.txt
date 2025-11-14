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

# Put the built library in:
set(_iface_outdir "${CMAKE_BINARY_DIR}/src/interface")
set_target_properties(${lib} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY  "${_iface_outdir}"  # .dll on Windows
  LIBRARY_OUTPUT_DIRECTORY  "${_iface_outdir}"  # .so/.dylib on Linux/macOS
  ARCHIVE_OUTPUT_DIRECTORY  "${_iface_outdir}"  # import/static .lib/.a
  OUTPUT_NAME               "svzero_interface"  # basename stays consistent
)

# Add the 'lib' prefix
if (WIN32)
  set_target_properties(${lib} PROPERTIES
    PREFIX "lib"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
  )
endif()
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

# Build the model object files.

set(lib svzero_model_library)

set(CXXSRCS 
  Block.cpp
  BloodVessel.cpp
  BloodVesselJunction.cpp
  ChamberElastanceInductor.cpp
  ClosedLoopCoronaryBC.cpp
  ClosedLoopCoronaryLeftBC.cpp
  ClosedLoopCoronaryRightBC.cpp
  ClosedLoopHeartPulmonary.cpp
  ClosedLoopRCRBC.cpp
  DOFHandler.cpp
  FlowReferenceBC.cpp
  Junction.cpp
  Model.cpp
  Node.cpp
  OpenLoopCoronaryBC.cpp
  Parameter.cpp
  PressureReferenceBC.cpp
  ResistanceBC.cpp
  ResistiveJunction.cpp
  ValveTanh.cpp
  WindkesselBC.cpp
)

set(HDRS 
  Block.h
  BlockType.h
  BloodVessel.h
  BloodVesselJunction.h
  ChamberElastanceInductor.h
  ClosedLoopCoronaryBC.h
  ClosedLoopCoronaryLeftBC.h
  ClosedLoopCoronaryRightBC.h
  ClosedLoopHeartPulmonary.h
  ClosedLoopRCRBC.h
  DOFHandler.h
  FlowReferenceBC.h
  Junction.h
  Model.h
  Node.h
  OpenLoopCoronaryBC.h
  Parameter.h
  PressureReferenceBC.h
  ResistanceBC.h
  ResistiveJunction.h
  ValveTanh.h
  WindkesselBC.h
)

add_library(${lib} OBJECT ${CXXSRCS} ${HDRS})

target_include_directories( ${lib} PUBLIC
  ${CMAKE_SOURCE_DIR}/src/algebra
  ${CMAKE_SOURCE_DIR}/src/model
  ${CMAKE_SOURCE_DIR}/src/solve
)

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )


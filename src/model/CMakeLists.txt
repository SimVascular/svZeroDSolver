
set(lib svzero_model_library)
message("##### svzero_model_library ######")

set(CXXSRCS 
  Block.cpp
  BloodVessel.cpp
  BloodVesselJunction.cpp
  ClosedLoopCoronaryBC.cpp
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
  WindkesselBC.cpp
)

set(HDRS 
Block.h
BlockType.h
BloodVessel.h
BloodVesselJunction.h
ClosedLoopCoronaryBC.h
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
WindkesselBC.h
)

add_library(${lib} STATIC ${CXXSRCS} )

include_directories( .. )

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )






set(lib svzero_algebra_library)

message("##### svzero_algebra_library ######")

set(CXXSRCS Integrator.cpp SparseSystem.cpp State.cpp )
set(HDRS Integrator.h SparseSystem.h State.h )

include_directories(..)

#add_library(${lib} STATIC ${CXXSRCS} )
add_library(${lib} SHARED ${CXXSRCS} )


target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )
target_link_libraries( ${lib} svzero_model_library)


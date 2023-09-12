
set(lib svzero_io_library)
message("##### svzero_io_library ######")

set(CXXSRCS SimulationParameters.cpp csvwriter.cpp io.cpp)

set(HDRS SimulationParameters.h csvwriter.h io.h)

add_library(${lib} SHARED ${CXXSRCS} )

include_directories( .. )

#target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )
#target_link_libraries( ${lib} svzero_algebra_library)
#target_link_libraries( ${lib} svzero_model_library)


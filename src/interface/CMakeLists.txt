
set(lib svzero_interface_library)
message("##### svzero_interface_library ######")

set(CXXSRCS interface.cpp )
set(HDRS interface.h )

add_library(${lib} SHARED ${CXXSRCS} )

include_directories( .. )

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json ) 

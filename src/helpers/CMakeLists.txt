
set(lib svzero_helpers_library)
message("##### svzero_helpers_library ######")

set(CXXSRCS helpers.cpp)

set(HDRS debug.h endswith.h helpers.h)

add_library(${lib} STATIC ${CXXSRCS} )

include_directories( .. )

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )





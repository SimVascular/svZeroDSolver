
set(lib svzero_interface_library)
message("##### svzero_interface_library ######")

set(CXXSRCS interface.cpp )
set(HDRS interface.h )

#find_package(Eigen3 REQUIRED NO_MODULE)
#find_package(simdjson REQUIRED)
#find_package(pybind11 REQUIRED)

#include_directories( ${EIGEN_INCLUDE_PATH} )

add_library(${lib} SHARED ${CXXSRCS} )

include_directories( .. )

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} simdjson::simdjson ) 
#target_link_libraries( ${lib} jsoncpp ) 



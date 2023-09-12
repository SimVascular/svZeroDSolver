
set(lib svzero_interface_library)
message("##### svzero_interface_library ######")

file(GLOB MODEL_SRC_FILES ../model/*.cpp)
file(GLOB IO_SRC_FILES ../io/*.cpp)
file(GLOB ALGEBRA_SRC_FILES ../algebra/*.cpp)

set(CXXSRCS interface.cpp )
set(HDRS interface.h )

add_library(${lib} SHARED ${CXXSRCS} ${MODEL_SRC_FILES} ${IO_SRC_FILES} ${ALGEBRA_SRC_FILES} )

include_directories( .. )

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json ) 
#target_link_libraries( ${lib} svzero_io_library)
#target_link_libraries( ${lib} svzero_model_library)


set(lib svzero_solve_library)
message("##### svzero_solve_library ######")

set(CXXSRCS Solver.cpp)

set(HDRS Solver.h)

add_library(${lib} SHARED ${CXXSRCS} )

include_directories( .. )

target_link_libraries( ${lib} Eigen3::Eigen )
target_link_libraries( ${lib} nlohmann_json::nlohmann_json )
target_link_libraries( ${lib} svzero_model_library) 
target_link_libraries( ${lib} svzero_algebra_library) 
target_link_libraries( ${lib} svzero_io_library) 
target_link_libraries( ${lib} svzero_optimize_library) 


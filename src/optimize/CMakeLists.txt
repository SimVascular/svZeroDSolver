
#set(lib svzero_optimize_library)
message("##### svzero_optimize_library ######")

set(CXXSRCS LevenbergMarquardtOptimizer.cpp calibrate.cpp)

set(HDRS LevenbergMarquardtOptimizer.h calibrate.h)

#add_library(${lib} STATIC ${CXXSRCS} )

include_directories( .. )

target_sources(svzerodsolver PRIVATE ${CXXSRCS})
target_sources(svzerodcalibrator PRIVATE ${CXXSRCS})

#target_link_libraries( ${lib} Eigen3::Eigen )
#target_link_libraries( ${lib} nlohmann_json::nlohmann_json )
#target_link_libraries( ${lib} svzero_model_library) 


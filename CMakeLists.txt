

project(svZeroDSolver)
cmake_minimum_required(VERSION 3.23)

add_compile_options(-std=c++17)

add_library(jsoncpp SHARED IMPORTED)
set_target_properties(jsoncpp PROPERTIES
    IMPORTED_LOCATION "/usr/local/lib/libjsoncpp.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include/json"
)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
add_executable(svzerodsolver src/main.cpp src/dofhandler.cpp src/node.cpp src/block.cpp src/junction.cpp src/bloodvessel.cpp src/flowreference.cpp src/rcrblockwithdistalpressure.cpp src/integrator.cpp src/system.cpp src/model.cpp src/parameter.cpp)

target_link_libraries(svzerodsolver Eigen3::Eigen)
target_link_libraries(svzerodsolver jsoncpp)

# brew install eigen
# brew install jsoncpp
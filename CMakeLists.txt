

project(svZeroDSolver)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
add_executable(svzerodsolver src/main.cpp src/dofhandler.cpp src/node.cpp src/block.cpp src/junction.cpp src/bloodvessel.cpp)

target_link_libraries(svzerodsolver Eigen3::Eigen)

# brew install eigen
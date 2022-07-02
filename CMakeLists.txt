cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(svZeroDSolver)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(simdjson REQUIRED)

add_library(jsoncpp SHARED IMPORTED)

if(APPLE)
    set_target_properties(jsoncpp PROPERTIES
        IMPORTED_LOCATION "/usr/local/lib/libjsoncpp.dylib"
        INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include/json"
    )
endif()

add_executable(svzerodsolver src/main.cpp)

target_link_libraries(svzerodsolver Eigen3::Eigen)
target_link_libraries(svzerodsolver simdjson::simdjson)
target_link_libraries(svzerodsolver jsoncpp)

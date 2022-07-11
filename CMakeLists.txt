cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(svZeroDSolver)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(simdjson REQUIRED)
find_package(pybind11 REQUIRED)

add_library(jsoncpp SHARED IMPORTED)

if(APPLE)
    set_target_properties(jsoncpp PROPERTIES
        IMPORTED_LOCATION "/usr/local/lib/libjsoncpp.dylib"
        INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include/json"
    )
elseif("${CMAKE_SYSTEM}" MATCHES "Linux")
    set_target_properties(jsoncpp PROPERTIES
        IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libjsoncpp.so"
        INTERFACE_INCLUDE_DIRECTORIES "/usr/include/jsoncpp/json"
    )
endif()

add_executable(svzerodsolver src/main.cpp)
pybind11_add_module(svzerodsolvercpp src/svzerodsolver.cpp)

target_link_libraries(svzerodsolver PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodsolver PRIVATE simdjson::simdjson)
target_link_libraries(svzerodsolver PRIVATE jsoncpp)
target_link_libraries(svzerodsolvercpp PRIVATE Eigen3::Eigen)
target_link_libraries(svzerodsolvercpp PRIVATE simdjson::simdjson)
target_link_libraries(svzerodsolvercpp PRIVATE jsoncpp)

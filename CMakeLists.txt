cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(svZeroDSolver)

set(CLUSTER OFF CACHE BOOL "Turn on cluster build.")

if (CLUSTER)
    set(EXTERNAL_LIB_HEADERS externals/eigen externals/jsoncpp/include/json externals/jsoncpp/include/ externals/simdjson/singleheader externals/pybind11)
    set(JSONCPP_INCLUDE_DIRS "$externals/jsoncpp/include/")
    file(GLOB JSONCPP_SOURCES
		"externals/jsoncpp/include/json/*.h"
		"externals/jsoncpp/src/lib_json/*.inl"
		"externals/jsoncpp/src/lib_json/*.h"
		"externals/jsoncpp/src/lib_json/*.cpp")
    add_executable(svzerodsolver src/main.cpp externals/simdjson/singleheader/simdjson.cpp ${JSONCPP_SOURCES})
    add_subdirectory(externals/pybind11)
    pybind11_add_module(svzerodsolvercpp src/svzerodsolver.cpp)
    target_include_directories(svzerodsolver PRIVATE ${EXTERNAL_LIB_HEADERS})
    target_include_directories(svzerodsolvercpp PRIVATE ${EXTERNAL_LIB_HEADERS})

else()

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

endif()
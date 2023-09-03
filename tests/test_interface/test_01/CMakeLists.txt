add_executable(svZeroD_interface_test01 ../LPNSolverInterface/LPNSolverInterface.cpp  main.cpp)
target_link_libraries(svZeroD_interface_test01 ${CMAKE_DL_LIBS})

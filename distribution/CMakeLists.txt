# Create distribution
if(ENABLE_DISTRIBUTION)
  include(CPackComponent)

  # set version (= today)
  string(TIMESTAMP TODAY "%Y-%m-%d")

  # add description (shown in installer)
  set(CPACK_PACKAGE_VERSION "${TODAY}")
  set(CPACK_PACKAGE_NAME "svZeroDSolver")
  set(CPACK_PACKAGE_VENDOR "SimVascular.org")
  set(CPACK_PACKAGE_CONTACT "SimVascular <sv@simvascular.org>")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "svZeroDSolver is a fast lumped-parameter solver for blood flow and pressure in hemodynamic networks")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/readme.txt")

  # define executale file name
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}")

  # define install location
  set(CPACK_PACKAGING_INSTALL_PREFIX "/usr/local/sv/svZeroDSolver/${CPACK_PACKAGE_VERSION}")

  # define install target within build folder
  set(CPACK_OUTPUT_FILE_PREFIX ${CMAKE_BINARY_DIR}/distribution)

  if(APPLE)
    set(CPACK_GENERATOR "productbuild")
  endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CPACK_GENERATOR DEB)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "")
  endif()

  # create the installer for both solver and calibrator
  install(TARGETS svzerodsolver svzerodcalibrator
          RUNTIME DESTINATION bin
          COMPONENT applications)

  # add the interface header
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../src/interface/interface.h DESTINATION include COMPONENT headers)

  # add interface library
  if(APPLE)
    install(FILES ${CMAKE_BINARY_DIR}/src/interface/libsvzero_interface.dylib DESTINATION lib COMPONENT libraries)
  elseif(LINUX)
    install(FILES ${CMAKE_BINARY_DIR}/src/interface/libsvzero_interface.so DESTINATION lib COMPONENT libraries)
  endif(APPLE)

  set(CPACK_COMPONENTS_ALL applications libraries headers)

  include(CPack)
endif()

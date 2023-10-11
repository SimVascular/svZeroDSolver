# Create distribution
if(ENABLE_DISTRIBUTION)
  include(CPackComponent)

  # add description
  set(CPACK_PACKAGE_VERSION "SimVascular-2023-10")
  set(CPACK_PACKAGE_NAME "svZeroDPlus")
  set(CPACK_PACKAGE_VENDOR "SimVascular.org")
  set(CPACK_PACKAGE_CONTACT "SimVascular<sv@simvascular.org>")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "svZeroDPlus is a fast lumped-parameter solver for blood flow and pressure in hemodynamic networks")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/readme.txt")
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}")

  if(APPLE)
    set(CPACK_GENERATOR "productbuild")
  elseif(LINUX)
    set(CPACK_GENERATOR DEB)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "")
  endif(APPLE)

  # define install targets
  set(CPACK_OUTPUT_FILE_PREFIX ${CMAKE_BINARY_DIR}/distribution)

  # create the installer (only solver)
  install(TARGETS svzerodsolver RUNTIME DESTINATION bin COMPONENT executable)
  include(CPack)
endif()
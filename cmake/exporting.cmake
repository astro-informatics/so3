# Exports so3 so other packages can access it
export(
  TARGETS astro-informatics-so3
  FILE "${PROJECT_BINARY_DIR}/astro-informatics-so3Targets.cmake"
  NAMESPACE astro-informatics-so3::)

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE astro-informatics-so3)
endif()

set(INCLUDE_INSTALL_DIR include/)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/so3Config.in.cmake
  "${PROJECT_BINARY_DIR}/astro-informatics-so3Config.cmake"
  INSTALL_DESTINATION lib/cmake/so3
  PATH_VARS INCLUDE_INSTALL_DIR)
write_basic_package_version_file(
  so3ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

if(NOT CONAN_EXPORTED)
  install(FILES "${PROJECT_BINARY_DIR}/astro-informatics-so3Config.cmake"
                "${PROJECT_BINARY_DIR}/astro-informatics-so3ConfigVersion.cmake"
          DESTINATION lib/cmake/astro-informatics-so3)
endif()

install(
  EXPORT so3Targets
  DESTINATION lib/cmake/astro-informatics-so3
  NAMESPACE astro-informatics-so3::)

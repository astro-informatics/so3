set(SO3_VERSION "@PROJECT_VERSION@")

@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(FFTW3 REQUIRED)
find_dependency(Ssht REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/astro-informatics-so3Targets.cmake")
set(SO3_LIBRARIES astro-informatics-so3::astro-informatics-so3)

check_required_components(SO3)

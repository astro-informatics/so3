cmake_minimum_required(VERSION 3.7)
option(openmp         "Enable OpenMP"                                   on)
option(fPIC           "Enable fPIC (only for gcc)"                      on)
option(tests          "Enable testing"                                  on)

project(so3 C)


set(SO3_VERSION "1.1b1")
set(SO3_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
find_package(OpenMP QUIET)
find_package(Ssht REQUIRED)

add_library(so3 STATIC src/c/so3_core.c src/c/so3_sampling.c src/c/so3_adjoint.c)
target_include_directories(so3 PUBLIC ${FFTW3_INCLUDE_DIRS})
target_include_directories(so3 PUBLIC ${Ssht_INCLUDE_DIR})
target_include_directories(so3 PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>
  $<INSTALL_INTERFACE:include/>)
target_include_directories(so3 PRIVATE ${PROJECT_SOURCE_DIR}/include/)

if(OpenMP_C_FOUND AND openmp)
  target_link_libraries(so3 OpenMP::OpenMP_C)
endif()
target_link_libraries(so3 ${Ssht_LIBRARY})
target_link_libraries(so3 ${FFTW3_DOUBLE_SERIAL_LIBRARY})
target_link_libraries(so3 m)
if(fPIC)
  set_target_properties(so3 PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS} -fPIC")
else()
  set_target_properties(so3 PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
endif()

install(TARGETS so3 EXPORT So3Targets ARCHIVE DESTINATION lib PUBLIC_HEADER)
install(FILES
  include/so3_error.h include/so3_sampling.h
  include/so3_core.h include/so3_types.h
  include/so3_adjoint.h
  include/so3.h
  DESTINATION include)

include("exporting")
#building about and test executables
add_executable(so3_about src/c/so3_about.c)
target_compile_definitions(so3_about PRIVATE SO3_VERSION="${SO3_VERSION}" SO3_BUILD="${SO3_BUILD}")
set_target_properties(so3_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if(tests)
  add_executable(so3_test src/c/so3_test.c)
  target_include_directories(so3_test PRIVATE ${PROJECT_SOURCE_DIR}/include/)
  target_link_libraries(so3_test so3)
  set_target_properties(so3_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  add_executable(so3_unittest src/c/unittest/so3_unittest.c)
  target_link_libraries(so3_unittest so3)
  target_include_directories(so3_unittest PRIVATE ${PROJECT_SOURCE_DIR}/include/)
  set_target_properties(so3_unittest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

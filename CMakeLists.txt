cmake_minimum_required(VERSION 3.7)
option(openmp         "Enable OpenMP"                  ON)
option(fPIC           "Enable fPIC (only for gcc)"     ON)
option(tests          "Enable testing"                 ON)
option(conan_fftw     "Download fftw using conan"      OFF)
option(conan_ssht     "Download ssht using conan"      OFF)
option(python         "Creates python package"         OFF)

project(so3 C)


set(SO3_VERSION "1.1b1")
set(SO3_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

if(conan_fftw)
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake")
    endif()
    include(${CMAKE_BINARY_DIR}/conan.cmake)

    if(NOT CONAN_OPTIONS)
        set(CONAN_OPTIONS
            "fftw:fPIC=True" "fftw:shared=False" "fftw:precision=double"
            "fftw:openmp=False")
    endif()
    if(NOT CONAN_BUILD)
        set(CONAN_BUILD "missing")
    endif()
    set(_GLIBCXX_USE_CXX11_ABI TRUE)

    conan_check(REQUIRED)
    conan_add_remote(
        URL https://api.bintray.com/conan/bincrafters/public-conan
        NAME bincrafters)
    conan_cmake_run(REQUIRES fftw/3.3.8@bincrafters/stable
        BASIC_SETUP
        OPTIONS "${CONAN_OPTIONS}"
        KEEP_RPATHS
        CMAKE_TARGETS
        NO_OUTPUT_DIRS
        BUILD ${CONAN_BUILD}
        # import headers to help vim find them
        IMPORTS "include, * -> conan/include")
endif()

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
find_package(Ssht PATHS ${CMAKE_CURRENT_BINARY_DIR}/ssht)
if(NOT Ssht_FOUND)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -S ${PROJECT_SOURCE_DIR}/cmake/ssht
            -B ${CMAKE_CURRENT_BINARY_DIR}/ssht
            -DFFTW3_LIBRARIES=${FFTW3_LIBRARIES}
            -DFFTW3_INCLUDE_DIR=${FFTW3_INCLUDE_DIR}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -G "${CMAKE_GENERATOR}"
            -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/ssht
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        COMMAND_ECHO STDOUT
    )
    execute_process(
        COMMAND ${CMAKE_COMMAND}
            --build ${CMAKE_CURRENT_BINARY_DIR}/ssht
        COMMAND_ECHO STDOUT
    )
    execute_process(
        COMMAND ${CMAKE_COMMAND}
            --install ${CMAKE_CURRENT_BINARY_DIR}/ssht
        COMMAND_ECHO STDOUT
    )
    find_package(Ssht PATHS ${CMAKE_CURRENT_BINARY_DIR}/ssht)
    set(SSHT_DOWNLOADED TRUE CACHE BOOL "ssht was downloaded")
endif()

add_library(
    so3
    STATIC
    src/c/so3_core.c src/c/so3_sampling.c src/c/so3_adjoint.c src/c/so3_conv.c
    src/c/so3_test_utils.c
)
target_include_directories(so3 PUBLIC ${FFTW3_INCLUDE_DIR})
if(SSHT_DOWNLOADED)
    target_include_directories(so3 PRIVATE ${Ssht_INCLUDE_DIRS})
else()
    target_include_directories(so3 PUBLIC ${Ssht_INCLUDE_DIRS})
endif()
target_link_libraries(so3 PUBLIC ssht)
if(fPIC)
  set_target_properties(so3 PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS} -fPIC")
else()
  set_target_properties(so3 PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
endif()

install(TARGETS so3 EXPORT So3Targets ARCHIVE DESTINATION lib PUBLIC_HEADER)
install(FILES
  src/c/so3_error.h src/c/so3_sampling.h
  src/c/so3_core.h src/c/so3_types.h
  src/c/so3_adjoint.h src/c/so3_conv.h
  src/c/so3_test.h src/c/so3_test_utils.h
  src/c/so3.h
  DESTINATION include)

include("exporting")
#building about and test executables
add_executable(so3_about src/c/so3_about.c)
target_compile_definitions(so3_about PRIVATE SO3_VERSION="${SO3_VERSION}" SO3_BUILD="${SO3_BUILD}")
set_target_properties(so3_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if(tests)
  add_executable(so3_test src/c/so3_test.c)
  if(SSHT_DOWNLOADED)
      target_include_directories(so3_test PRIVATE ${Ssht_INCLUDE_DIRS})
  endif()
  target_include_directories(so3_test PRIVATE ${FFTW3_INCLUDE_DIR})
  target_link_libraries(so3_test PRIVATE so3 ${FFTW3_LIBRARIES})
  set_target_properties(so3_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  add_executable(so3_unittest src/c/unittest/so3_unittest.c)
  if(SSHT_DOWNLOADED)
      target_include_directories(so3_unittest PRIVATE ${Ssht_INCLUDE_DIRS})
  endif()
  target_include_directories(so3_unittest PRIVATE ${FFTW3_INCLUDE_DIR})
  target_link_libraries(so3_unittest PRIVATE so3 ${FFTW3_LIBRARIES})
  set_target_properties(so3_unittest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

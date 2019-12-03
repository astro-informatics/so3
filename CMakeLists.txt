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

add_subdirectory(src/c)
# if(python)
    # add_subdirectory(src/pyso3)
# endif()

include("exporting")

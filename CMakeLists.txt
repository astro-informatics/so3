cmake_minimum_required(VERSION 3.12)
project(
  so3
  VERSION "1.3.7"
  DESCRIPTION "Fast and exact Wigner transforms"
  HOMEPAGE_URL "http://astro-informatics.github.io/so3/"
  LANGUAGES C)

option(conan_deps "Download ssht using conan" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

include(CTest)
if(conan_deps OR CONAN_EDITABLE_MODE)
  include("${PROJECT_SOURCE_DIR}/cmake/conan_dependencies.cmake")
endif()
if(EXISTS "${PROJECT_BINARY_DIR}/conan_paths.cmake")
  include("${PROJECT_BINARY_DIR}/conan_paths.cmake")
elseif(EXISTS "${PROJECT_BINARY_DIR}/FindFFTW3.cmake")
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}")
else()
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
find_package(ssht REQUIRED)
if(ssht_FOUND AND NOT TARGET ssht::ssht)
  message(FATAL_ERROR "NOT FOUND ${SSHT_LIBRARIES}")
endif()
find_package(FFTW3 REQUIRED)
find_library(MATH_LIBRARY m)

add_subdirectory(src/c)
if(BUILD_TESTING)
  include("${PROJECT_SOURCE_DIR}/cmake/fetch_cmocka.cmake")
  add_subdirectory(tests)
endif()

if(SKBUILD)
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/so3)
elseif(NOT CONAN_EXPORTED)
  include("${PROJECT_SOURCE_DIR}/cmake/exporting.cmake")
endif()

# only run documentation if this is not a sub-project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
    set(DOXYGEN_EXTRACT_ALL "YES")
    set(DOXYGEN_FILE_PATTERNS *.h *.c)
    doxygen_add_docs(docs src/c tests/ include)
  endif()
endif()

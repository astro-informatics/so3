cmake_minimum_required(VERSION 3.12)
project(
  so3
  VERSION "1.3.0"
  DESCRIPTION "Fast and exact Wigner transforms"
  HOMEPAGE_URL "http://astro-informatics.github.io/so3/"
  LANGUAGES C)

option(tests "Enable testing" ON)
option(python "Creates python package only" OFF)
option(conan_deps "Download ssht using conan" ON)
if(NOT WIN32)
  option(fPIC "Enable position independent code" ON)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(conan_deps)
  include("conan_dependencies")
  set(SSHT_TARGET CONAN_PKG::ssht)
else()
  find_package(Ssht REQUIRED)
  set(SSHT_TARGET ssht)
endif()

add_subdirectory(src/c)
if(tests)
  enable_testing()
  include("fetch_cmocka")
  add_subdirectory(tests)
endif()

if(NOT python)
  include("exporting")
else()
  include("python-setup")
  find_python()
  setup_skbuild()
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/so3)
endif()

find_package(Doxygen)
if(DOXYGEN_FOUND)
  doxygen_add_docs(docs WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/docs)
endif()

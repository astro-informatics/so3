cmake_minimum_required(VERSION 3.10)
option(fPIC "Enable position independent code" ON)
option(tests "Enable testing" ON)
option(conan_deps "Download ssht using conan" ON)

project(
  so3
  VERSION "1.2.1"
  DESCRIPTION "Fast and exact Wigner transforms"
  HOMEPAGE_URL "http://astro-informatics.github.io/so3/"
  LANGUAGES C)

if(tests)
    enable_testing()
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

set(CMAKE_C_STANDARD 99)
set(SO3_BUILD "CMake")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(conan_deps)
  include("conan_dependencies")
  set(SSHT_TARGET CONAN_PKG::ssht)
else()
  find_package(Ssht REQUIRED)
  set(SSHT_TARGET ssht)
endif()

add_subdirectory(src/c)

include("exporting")

find_package(Doxygen)
if(DOXYGEN_FOUND)
    doxygen_add_docs(docs WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/docs)
endif()

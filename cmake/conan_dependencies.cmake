set(_conan_provider_file "${CMAKE_BINARY_DIR}/conan_provider.cmake")
set(_conan_provider_url
    "https://raw.githubusercontent.com/conan-io/cmake-conan/0.19.0/conan_provider.cmake")
set(_conan_provider_sha256
    "d574ac4f1ad0784e743304fde4af6e6692191783afae4898974903b622928546")
set(_conan_provider_hash "")

if(EXISTS "${_conan_provider_file}")
  file(SHA256 "${_conan_provider_file}" _conan_provider_hash)
endif()

if(NOT "${_conan_provider_hash}" STREQUAL "${_conan_provider_sha256}")
  file(DOWNLOAD
       "${_conan_provider_url}"
       "${_conan_provider_file}"
       EXPECTED_HASH "SHA256=${_conan_provider_sha256}"
       TLS_VERIFY ON
       STATUS _conan_provider_status)
  list(GET _conan_provider_status 0 _conan_provider_status_code)
  list(GET _conan_provider_status 1 _conan_provider_status_message)
  if(NOT _conan_provider_status_code EQUAL 0)
    file(REMOVE "${_conan_provider_file}")
    message(FATAL_ERROR
      "Failed to download Conan dependency provider: ${_conan_provider_status_message}")
  endif()
endif()

list(PREPEND CMAKE_PROJECT_TOP_LEVEL_INCLUDES "${_conan_provider_file}")

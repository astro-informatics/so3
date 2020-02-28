add_library(
    so3
    STATIC
    so3_core.c so3_sampling.c so3_adjoint.c so3_conv.c so3_test_utils.c
)
target_include_directories(so3 PUBLIC ${FFTW3_INCLUDE_DIR})
if(SSHT_DOWNLOADED)
    target_include_directories(so3 PRIVATE ${Ssht_INCLUDE_DIRS})
else()
    target_include_directories(so3 PUBLIC ${Ssht_INCLUDE_DIRS})
endif()
target_link_libraries(so3 PUBLIC ssht ${FFTW3_LIBRARIES})
if(fPIC)
  set_target_properties(so3 PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS} -fPIC")
else()
  set_target_properties(so3 PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
endif()

install(TARGETS so3 EXPORT So3Targets ARCHIVE DESTINATION lib PUBLIC_HEADER)
install(FILES
  so3_error.h so3_sampling.h so3_core.h so3_types.h so3_adjoint.h so3_conv.h
  so3.h
  DESTINATION include)

#building about and test executables
add_executable(so3_about so3_about.c)
target_compile_definitions(so3_about PRIVATE SO3_VERSION="${SO3_VERSION}" SO3_BUILD="${SO3_BUILD}")
set_target_properties(so3_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if(tests)
  add_executable(so3_test so3_test.c)
  if(SSHT_DOWNLOADED)
      target_include_directories(so3_test PRIVATE ${Ssht_INCLUDE_DIRS})
  endif()
  target_include_directories(so3_test PRIVATE ${FFTW3_INCLUDE_DIR})
  target_link_libraries(so3_test PRIVATE so3 ${FFTW3_LIBRARIES})
  set_target_properties(so3_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  add_executable(so3_unittest unittest/so3_unittest.c)
  if(SSHT_DOWNLOADED)
      target_include_directories(so3_unittest PRIVATE ${Ssht_INCLUDE_DIRS})
  endif()
  target_include_directories(so3_unittest PRIVATE ${FFTW3_INCLUDE_DIR})
  target_link_libraries(so3_unittest PRIVATE so3 ${FFTW3_LIBRARIES})
  set_target_properties(so3_unittest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

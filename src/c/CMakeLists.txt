add_library(so3 STATIC so3_core.c so3_sampling.c so3_adjoint.c)
target_link_libraries(so3 PUBLIC ${SSHT_TARGET})
target_include_directories(
  so3 PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
             $<INSTALL_INTERFACE:include>)

if(fPIC)
  set_target_properties(so3 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

install(
  TARGETS so3
  EXPORT So3Targets
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER)

install(
  FILES ${PROJECT_SOURCE_DIR}/include/so3.h
        ${PROJECT_SOURCE_DIR}/include/so3_adjoint.h
        ${PROJECT_SOURCE_DIR}/include/so3_core.h
        ${PROJECT_SOURCE_DIR}/include/so3_error.h
        ${PROJECT_SOURCE_DIR}/include/so3_sampling.h
        ${PROJECT_SOURCE_DIR}/include/so3_types.h
  DESTINATION include/so3)

# building about and test executables
add_executable(so3_about so3_about.c)
target_compile_definitions(so3_about PRIVATE SO3_VERSION="${SO3_VERSION}"
                                             SO3_BUILD="${SO3_BUILD}")
set_target_properties(so3_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                           ${PROJECT_BINARY_DIR})
if(tests)
  add_executable(so3_test so3_test.c)
  target_link_libraries(so3_test so3)
  set_target_properties(so3_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                            ${PROJECT_BINARY_DIR})
  add_test(NAME so3_test COMMAND so3_test)

  add_executable(so3_unittest unittest/so3_unittest.c)
  target_link_libraries(so3_unittest so3)
  set_target_properties(so3_unittest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                ${PROJECT_BINARY_DIR})
  add_test(NAME so3_unittest COMMAND so3_unittest)
endif()

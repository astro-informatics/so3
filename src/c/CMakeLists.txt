add_library(so3 STATIC so3_core.c so3_sampling.c so3_adjoint.c so3_conv.c)
target_link_libraries(so3 PUBLIC ${SSHT_TARGET})
target_include_directories(
  so3
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
         $<INSTALL_INTERFACE:include>)
set_target_properties(so3 PROPERTIES C_STANDARD 99)
configure_file(${PROJECT_SOURCE_DIR}/include/so3/so3_version.in.h
               ${PROJECT_BINARY_DIR}/include/so3/so3_version.h)

if(fPIC)
  set_target_properties(so3 PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

install(
  TARGETS so3
  EXPORT So3Targets
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER)

install(
  FILES ${PROJECT_SOURCE_DIR}/include/so3/so3.h
        ${PROJECT_SOURCE_DIR}/include/so3/so3_adjoint.h
        ${PROJECT_SOURCE_DIR}/include/so3/so3_core.h
        ${PROJECT_SOURCE_DIR}/include/so3/so3_conv.h
        ${PROJECT_SOURCE_DIR}/include/so3/so3_error.h
        ${PROJECT_SOURCE_DIR}/include/so3/so3_sampling.h
        ${PROJECT_SOURCE_DIR}/include/so3/so3_types.h
        ${PROJECT_BINARY_DIR}/include/so3/so3_version.h
  DESTINATION include/so3)

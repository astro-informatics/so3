add_library(astro-informatics-so3 STATIC so3_core.c so3_sampling.c
                                         so3_adjoint.c so3_conv.c)
target_link_libraries(astro-informatics-so3 PUBLIC ssht::ssht FFTW3::FFTW3
                                                   ${MATH_LIBRARY})
target_include_directories(
  astro-informatics-so3
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
         $<INSTALL_INTERFACE:include>)
set_target_properties(
  astro-informatics-so3 PROPERTIES C_STANDARD 99 ARCHIVE_OUTPUT_DIRECTORY
                                                 ${PROJECT_BINARY_DIR}/lib)
target_compile_features(astro-informatics-so3 PUBLIC c_std_99)
configure_file(${PROJECT_SOURCE_DIR}/include/so3/so3_version.in.h
               ${PROJECT_BINARY_DIR}/include/so3/so3_version.h)

if(NOT SKBUILD)
  install(
    TARGETS astro-informatics-so3
    EXPORT so3Targets
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
endif()

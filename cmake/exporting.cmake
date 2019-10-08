# Exports so3 so other packages can access it
export(TARGETS so3 FILE "${PROJECT_BINARY_DIR}/So3Targets.cmake")

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE so3)
endif()

# First in binary dir
set(ALL_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")
configure_File(cmake/So3Config.in.cmake
  "${PROJECT_BINARY_DIR}/So3Config.cmake" @ONLY
)
configure_File(cmake/So3ConfigVersion.in.cmake
  "${PROJECT_BINARY_DIR}/So3ConfigVersion.cmake" @ONLY
)

# Then for installation tree
file(RELATIVE_PATH REL_INCLUDE_DIR
    "${CMAKE_INSTALL_PREFIX}/share/cmake/so3"
    "${CMAKE_INSTALL_PREFIX}/include/so3"
)
set(ALL_INCLUDE_DIRS "\${So3_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(cmake/So3Config.in.cmake
  "${PROJECT_BINARY_DIR}/CMakeFiles/So3Config.cmake" @ONLY
)

# Finally install all files
install(FILES
  "${PROJECT_BINARY_DIR}/CMakeFiles/So3Config.cmake"
  "${PROJECT_BINARY_DIR}/So3ConfigVersion.cmake"
    DESTINATION share/cmake/so3
    COMPONENT dev
)

install(EXPORT So3Targets DESTINATION share/cmake/so3 COMPONENT dev)

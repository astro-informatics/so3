# add command to cythonize results
add_cython_target(cython_sources pyso3.pyx C PY3)

Python3_add_library(pyso3 MODULE ${cython_sources})
target_include_directories(pyso3
    PUBLIC "${PROJECT_SOURCE_DIR}/include/c" ${NumPy_INCLUDE_DIRS})
if(SSHT_DOWNLOADED)
    target_include_directories(pyso3 PRIVATE ${Ssht_INCLUDE_DIRS})
else()
    target_include_directories(pyso3 PUBLIC ${Ssht_INCLUDE_DIRS})
endif()
target_link_libraries(pyso3 PUBLIC so3)
install(TARGETS pyso3 LIBRARY DESTINATION src/pyso3)

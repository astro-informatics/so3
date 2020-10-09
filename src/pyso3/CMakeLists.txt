# add command to cythonize results
add_cython_target(cython_sources pyso3.pyx C PY3)

python3_add_library(pyso3 MODULE WITH_SOABI ${cython_sources})
target_link_libraries(pyso3 PUBLIC so3)
target_include_directories(pyso3 PUBLIC ${NumPy_INCLUDE_DIRS})
if(WIN32)
    target_compile_definitions(pyso3 PUBLIC MS_WIN64)
endif()
install(TARGETS pyso3 LIBRARY DESTINATION src/pyso3)

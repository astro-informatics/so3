# add command to cythonize results
add_cython_target(cython_sources bindings.pyx C PY3)

python3_add_library(bindings MODULE WITH_SOABI ${cython_sources})
target_link_libraries(bindings PUBLIC so3)
target_include_directories(bindings PUBLIC ${NumPy_INCLUDE_DIRS})
if(WIN32)
    target_compile_definitions(bindings PUBLIC MS_WIN64)
endif()
install(TARGETS bindings LIBRARY DESTINATION src/so3)

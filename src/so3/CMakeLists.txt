# add command to cythonize results
add_cython_target(cython_sources bindings.pyx C PY3)

add_library(bindings MODULE ${cython_sources})
python_extension_module(bindings)
target_link_libraries(bindings astro-informatics-so3)
target_include_directories(bindings PUBLIC ${NumPy_INCLUDE_DIRS})
if(WIN32)
  target_compile_definitions(bindings PUBLIC MS_WIN64)
endif()
install(TARGETS bindings LIBRARY DESTINATION src/so3)

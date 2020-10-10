add_library(utilities OBJECT utilities.c)
target_link_libraries(utilities PUBLIC so3)
foreach(testname sampling so3 convolution)
  add_executable(test_${testname} test_${testname}.c)
  target_link_libraries(test_${testname} PRIVATE so3 cmocka)
  set_target_properties(
    test_${testname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
                                C_STANDARD 11)
  add_test(NAME test_${testname} COMMAND test_${testname})
endforeach()
foreach(testname so3 convolution)
    target_link_libraries(test_${testname} PRIVATE utilities)
endforeach()

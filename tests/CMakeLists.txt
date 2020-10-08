add_executable(test_sampling test_sampling.c)
target_link_libraries(test_sampling so3 cmocka)
set_target_properties(
  test_sampling PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
                           C_STANDARD 11)
add_test(NAME test_sampling COMMAND test_sampling)

add_executable(test_so3 test_so3.c utilities.c)
target_link_libraries(test_so3 so3 cmocka)
set_target_properties(test_so3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                          ${PROJECT_BINARY_DIR} C_STANDARD 11)
add_test(NAME test_so3 COMMAND test_so3)

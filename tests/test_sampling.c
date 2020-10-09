#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>

#include "so3/so3_sampling.h"
#include "so3/so3_types.h"

#include <cmocka.h>

typedef struct {
  int el, m, n, index;
} LMN2Index;

void sampling_assertions(
    LMN2Index const *tests, int n, so3_parameters_t const *parameters) {
  for (int i = 0; i < n; i += 1) {
    int ind;
    so3_sampling_elmn2ind(&ind, tests[i].el, tests[i].m, tests[i].n, parameters);
    assert_int_equal(ind, tests[i].index);
  }

  for (int i = 0; i < n; i += 1) {
    int el, m, n;
    so3_sampling_ind2elmn(&el, &m, &n, tests[i].index, parameters);
    assert_int_equal(el, tests[i].el);
    assert_int_equal(m, tests[i].m);
    assert_int_equal(n, tests[i].n);
  }
}

void test_indexing_padded_storage_with_zero_first(void **state) {
  // note that implementation does not depend on N
  so3_parameters_t parameters = {
      .L = 3,
      .N = -1,
      .n_order = SO3_N_ORDER_ZERO_FIRST,
      .storage = SO3_STORAGE_PADDED};
  LMN2Index tests[] = {
      {.el = 0, .m = 0, .n = 0, .index = 0},
      {.el = 2, .m = 1, .n = 0, .index = 7},
      {.el = 2, .m = 1, .n = -1, .index = 16},
      {.el = 2, .m = 1, .n = 1, .index = 25},
  };

  sampling_assertions(tests, sizeof(tests) / sizeof(tests[0]), &parameters);
}

void test_indexing_padded_storage_with_negative_first(void **state) {
  so3_parameters_t parameters = {
      .L = 3,
      .N = 3,
      .n_order = SO3_N_ORDER_NEGATIVE_FIRST,
      .storage = SO3_STORAGE_PADDED};
  LMN2Index tests[] = {
      {.el = 0, .m = 0, .n = 0, .index = 18},
      {.el = 2, .m = -2, .n = -2, .index = 4},
      {.el = 2, .m = 1, .n = 0, .index = 25},
      {.el = 2, .m = 1, .n = -1, .index = 16},
      {.el = 2, .m = 1, .n = 1, .index = 34},
  };

  sampling_assertions(tests, sizeof(tests) / sizeof(tests[0]), &parameters);
}

void test_indexing_compact_storage_with_zero_first(void **state) {
  // note that implementation does not depend on N
  so3_parameters_t parameters = {
      .L = 3,
      .N = -1,
      .n_order = SO3_N_ORDER_ZERO_FIRST,
      .storage = SO3_STORAGE_COMPACT};

  LMN2Index tests[] = {
      {.el = 0, .m = 0, .n = 0, .index = 0},
      {.el = 2, .m = 1, .n = 0, .index = 7},
      {.el = 2, .m = 1, .n = -1, .index = 15},
      {.el = 2, .m = 1, .n = 1, .index = 23},
      {.el = 2, .m = 1, .n = -2, .index = 28},
      {.el = 2, .m = 1, .n = 2, .index = 33},
  };

  sampling_assertions(tests, sizeof(tests) / sizeof(tests[0]), &parameters);
}

void test_indexing_compact_storage_with_negative_first(void **state) {
  so3_parameters_t parameters = {
      .L = 3,
      .N = 3,
      .n_order = SO3_N_ORDER_NEGATIVE_FIRST,
      .storage = SO3_STORAGE_COMPACT};
  LMN2Index tests[] = {
      {.el = 0, .m = 0, .n = 0, .index = 13},
      {.el = 2, .m = -2, .n = -2, .index = 0},
      {.el = 2, .m = 1, .n = 0, .index = 20},
      {.el = 2, .m = 1, .n = -1, .index = 11},
      {.el = 2, .m = 1, .n = 1, .index = 28},
      {.el = 2, .m = 1, .n = -2, .index = 3},
      {.el = 2, .m = 1, .n = 2, .index = 33},
  };

  sampling_assertions(tests, sizeof(tests) / sizeof(tests[0]), &parameters);
}

void sampling_real_assertions(
    LMN2Index const *tests, int n, so3_parameters_t const *parameters) {
  for (int i = 0; i < n; i += 1) {
    int ind;
    so3_sampling_elmn2ind_real(&ind, tests[i].el, tests[i].m, tests[i].n, parameters);
    assert_int_equal(ind, tests[i].index);
  }

  for (int i = 0; i < n; i += 1) {
    int el, m, n;
    so3_sampling_ind2elmn_real(&el, &m, &n, tests[i].index, parameters);
    assert_int_equal(el, tests[i].el);
    assert_int_equal(m, tests[i].m);
    assert_int_equal(n, tests[i].n);
  }
}

void test_indexing_real_padded_storage(void **state) {
  // note that implementation does not depend on N
  so3_parameters_t parameters = {.L = 3, .N = 3, .storage = SO3_STORAGE_PADDED};
  LMN2Index tests[] = {
      {.el = 0, .m = 0, .n = 0, .index = 0},
      {.el = 2, .m = 1, .n = 0, .index = 7},
      {.el = 2, .m = 1, .n = 1, .index = 16},
  };

  sampling_real_assertions(tests, sizeof(tests) / sizeof(tests[0]), &parameters);
}

void test_indexing_real_compact_storage(void **state) {
  so3_parameters_t parameters = {.L = 3, .N = 3, .storage = SO3_STORAGE_COMPACT};
  LMN2Index tests[] = {
      {.el = 0, .m = 0, .n = 0, .index = 0},
      {.el = 2, .m = 1, .n = 0, .index = 7},
      {.el = 2, .m = 1, .n = 1, .index = 15},
      {.el = 2, .m = 1, .n = 2, .index = 20},
  };

  sampling_real_assertions(tests, sizeof(tests) / sizeof(tests[0]), &parameters);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_indexing_padded_storage_with_zero_first),
      cmocka_unit_test(test_indexing_padded_storage_with_negative_first),
      cmocka_unit_test(test_indexing_compact_storage_with_zero_first),
      cmocka_unit_test(test_indexing_compact_storage_with_negative_first),
      cmocka_unit_test(test_indexing_real_padded_storage),
      cmocka_unit_test(test_indexing_real_compact_storage),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}

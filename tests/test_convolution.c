#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>

#include "so3/so3_conv.h"
#include "so3/so3_error.h"
#include "so3/so3_sampling.h"
#include "so3/so3_types.h"
#include "utilities.h"

#include <cmocka.h>

static void
gen_flmn(complex double *flmn, const so3_parameters_t *parameters, int seed) {
  if (parameters->reality)
    gen_flmn_real(flmn, parameters, seed);
  else
    gen_flmn_complex(flmn, parameters, seed);
}
static void
elmn2ind(int *ind, int el, int m, int n, const so3_parameters_t *parameters) {
  if (parameters->reality)
    so3_sampling_elmn2ind_real(ind, el, m, n, parameters);
  else
    so3_sampling_elmn2ind(ind, el, m, n, parameters);
}
static void
ind2elmn(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters) {
  if (parameters->reality)
    so3_sampling_ind2elmn_real(el, m, n, ind, parameters);
  else
    so3_sampling_ind2elmn(el, m, n, ind, parameters);
}

static void _n_loop_values_nmode_all(
    const so3_parameters_t *const parameters, int start, int stop, int inc) {
  int n_start_old, n_stop_old, n_inc_old;
  so3_sampling_n_loop_values(&n_start_old, &n_stop_old, &n_inc_old, parameters);
  assert_int_equal(start, n_start_old);
  assert_int_equal(stop, n_stop_old);
  assert_int_equal(inc, n_inc_old);
}

static void test_n_loop_values_real_nmode_all(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 1;
  parameters.n_mode = SO3_N_MODE_ALL;
  _n_loop_values_nmode_all(&parameters, 0, parameters.N - 1, 1);
}

static void test_n_loop_values_real_nmode_l(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 1;
  parameters.n_mode = SO3_N_MODE_L;
  _n_loop_values_nmode_all(&parameters, 0, parameters.N - 1, 1);
}

static void test_n_loop_values_real_nmode_even(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 1;
  parameters.n_mode = SO3_N_MODE_EVEN;
  _n_loop_values_nmode_all(
      &parameters,
      0,
      ((parameters.N - 1) % 2 == 0) ? parameters.N - 1 : parameters.N - 2,
      2);
}

static void test_n_loop_values_real_nmode_odd(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 1;
  parameters.n_mode = SO3_N_MODE_ODD;
  _n_loop_values_nmode_all(
      &parameters,
      1,
      ((parameters.N - 1) % 2 != 0) ? parameters.N - 1 : parameters.N - 2,
      2);
}

static void test_n_loop_values_real_nmode_maximum(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 1;
  parameters.n_mode = SO3_N_MODE_MAXIMUM;
  _n_loop_values_nmode_all(&parameters, parameters.N - 1, parameters.N - 1, 1);
}

static void test_n_loop_values_complex_nmode_all(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 0;
  parameters.n_mode = SO3_N_MODE_ALL;
  _n_loop_values_nmode_all(&parameters, -parameters.N + 1, parameters.N - 1, 1);
}

static void test_n_loop_values_complex_nmode_l(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 0;
  parameters.n_mode = SO3_N_MODE_L;
  _n_loop_values_nmode_all(&parameters, -parameters.N + 1, parameters.N - 1, 1);
}

static void test_n_loop_values_complex_nmode_even(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 0;
  parameters.n_mode = SO3_N_MODE_EVEN;
  _n_loop_values_nmode_all(
      &parameters,
      ((parameters.N - 1) % 2 == 0) ? -parameters.N + 1 : -parameters.N + 2,
      ((parameters.N - 1) % 2 == 0) ? parameters.N - 1 : parameters.N - 2,
      2);
}

static void test_n_loop_values_complex_nmode_odd(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 0;
  parameters.n_mode = SO3_N_MODE_ODD;
  _n_loop_values_nmode_all(
      &parameters,
      ((parameters.N - 1) % 2 != 0) ? -parameters.N + 1 : -parameters.N + 2,
      ((parameters.N - 1) % 2 != 0) ? parameters.N - 1 : parameters.N - 2,
      2);
}

static void test_n_loop_values_complex_nmode_maximum(void **state) {
  so3_parameters_t parameters = *(so3_parameters_t *)state;
  parameters.reality = 0;
  parameters.n_mode = SO3_N_MODE_MAXIMUM;
  _n_loop_values_nmode_all(
      &parameters,
      -parameters.N + 1,
      parameters.N - 1,
      parameters.N > 1 ? 2 * parameters.N - 2 : 1);
}

static void _harmonic_convolution(const so3_parameters_t *const parameters) {
  SO3_COMPLEX(double) * hlmn, *flmn, *glmn, h_value;
  int el, m, n;
  int seed = 1;
  int n_start, n_stop, n_inc;
  int ind_f, ind_g;

  const so3_parameters_t h_parameters =
      so3_conv_get_parameters_of_convolved_lmn(parameters, parameters);

  int flmn_length = so3_sampling_flmn_size(parameters);
  flmn = malloc(flmn_length * sizeof *flmn);
  SO3_ERROR_MEM_ALLOC_CHECK(flmn);

  int glmn_length = so3_sampling_flmn_size(parameters);
  glmn = malloc(glmn_length * sizeof *glmn);
  SO3_ERROR_MEM_ALLOC_CHECK(glmn);

  int hlmn_length = so3_sampling_flmn_size(&h_parameters);
  hlmn = malloc(hlmn_length * sizeof *hlmn);
  SO3_ERROR_MEM_ALLOC_CHECK(hlmn);

  gen_flmn(flmn, parameters, seed);
  gen_flmn(glmn, parameters, seed + 1);

  so3_conv_harmonic_convolution(
      hlmn, &h_parameters, flmn, parameters, glmn, parameters);

  for (int i = 0; i < flmn_length; i++) {
    h_value = 0;
    ind2elmn(&el, &m, &n, i, &h_parameters);
    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, parameters);
    for (int k = n_start; k <= n_stop; k += n_inc) {
      if (so3_sampling_is_elmn_non_zero(el, n, k, parameters)) {

        elmn2ind(&ind_f, el, m, k, parameters);
        elmn2ind(&ind_g, el, n, k, parameters);

        h_value += flmn[ind_f] * conj(glmn[ind_g]);
      }
    }
    assert_float_equal(creal(h_value), creal(hlmn[i]), 1e-12);
    assert_float_equal(cimag(h_value), cimag(hlmn[i]), 1e-12);
  }
}

static void test_harmonic_convolution_real(void **state) {
  so3_parameters_t *const parameters = (so3_parameters_t *)state;
  parameters->reality = 1;
  _harmonic_convolution(parameters);
}
static void test_harmonic_convolution_complex(void **state) {
  so3_parameters_t *const parameters = (so3_parameters_t *)state;
  parameters->reality = 0;
  _harmonic_convolution(parameters);
}

void test_conv_get_parameters_of_convolved_lmn(void **_) {
  const so3_parameters_t f_parameters = {
      .L0 = 1,
      .L = 3,
      .N = 2,
      .n_mode = 0,
      .reality = 0,
      .sampling_scheme = 0,
      .n_order = 0,
      .storage = 0,
      .steerable = 0,
  };
  so3_parameters_t g_parameters = f_parameters;
  g_parameters.L = 4;
  g_parameters.L0 = 0;
  g_parameters.N = 4;

  const so3_parameters_t h_parameters =
      so3_conv_get_parameters_of_convolved_lmn(&f_parameters, &g_parameters);

  assert_int_equal(h_parameters.L, 3);
  assert_int_equal(h_parameters.L0, 1);
  assert_int_equal(h_parameters.N, 3);
}

int main(void) {
  const so3_parameters_t states[] = {
      {.L0 = 0,
       .L = 8,
       .N = 8,
       .verbosity = 1,
       .n_mode = 0,
       .reality = 0,
       .sampling_scheme = 0,
       .n_order = 0,
       .storage = 0,
       .steerable = 0},
      {.L0 = 0,
       .L = 8,
       .N = 7,
       .verbosity = 1,
       .n_mode = 0,
       .reality = 0,
       .sampling_scheme = 0,
       .n_order = 0,
       .storage = 0,
       .steerable = 0},
      {.L0 = 0,
       .L = 8,
       .N = 0,
       .verbosity = 1,
       .n_mode = 0,
       .reality = 0,
       .sampling_scheme = 0,
       .n_order = 0,
       .storage = 0,
       .steerable = 0},
  };
  const struct CMUnitTest tests[] = {
      cmocka_unit_test_prestate(test_harmonic_convolution_real, (void **)states),
      cmocka_unit_test_prestate(test_harmonic_convolution_complex, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_real_nmode_all, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_real_nmode_l, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_real_nmode_even, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_real_nmode_odd, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_real_nmode_maximum, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_complex_nmode_all, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_complex_nmode_l, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_complex_nmode_even, (void **)states),
      cmocka_unit_test_prestate(test_n_loop_values_complex_nmode_odd, (void **)states),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_maximum, (void **)states),
      cmocka_unit_test_prestate(
          test_n_loop_values_real_nmode_all, (void **)(states + 1)),
      cmocka_unit_test_prestate(test_n_loop_values_real_nmode_l, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_real_nmode_even, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_real_nmode_odd, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_real_nmode_maximum, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_all, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_l, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_even, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_odd, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_maximum, (void **)(states + 1)),
      cmocka_unit_test_prestate(
          test_n_loop_values_real_nmode_maximum, (void **)(states + 2)),
      cmocka_unit_test_prestate(
          test_n_loop_values_complex_nmode_maximum, (void **)(states + 2)),
      cmocka_unit_test(test_conv_get_parameters_of_convolved_lmn),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}

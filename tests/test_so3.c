#include <assert.h>
#include <complex.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "so3/so3.h"
#include "utilities.h"

#include <cmocka.h>

typedef void (*inverse_real_t)(
    double *f, const complex double *flmn, const so3_parameters_t *parameters);
typedef void (*inverse_complex_t)(
    complex double *f, const complex double *flmn, const so3_parameters_t *parameters);
typedef void (*forward_real_t)(
    complex double *flmn, const double *f, const so3_parameters_t *parameters);
typedef void (*forward_complex_t)(
    complex double *flmn, const complex double *f, const so3_parameters_t *parameters);

typedef struct {
  so3_parameters_t params;
  int seed;
  double tolerance;
  inverse_real_t inverse_real;
  inverse_complex_t inverse_complex;
  forward_real_t forward_real;
  forward_complex_t forward_complex;
} SO3TestState;

void test_real_back_and_forth(void **_state) {
  SO3TestState const *state = *(const SO3TestState **)_state;

  so3_parameters_t const *parameters = &state->params;

  complex double *flmn_orig = calloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L, sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(flmn_orig);
  gen_flmn_real(flmn_orig, parameters, state->seed);

  double *f_image = calloc(
      (2 * parameters->L) * (parameters->L + 1) * (2 * parameters->N - 1),
      sizeof(double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_image);
  state->inverse_real(f_image, flmn_orig, parameters);

  complex double *f_back = calloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L, sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_back);
  state->forward_real(f_back, f_image, parameters);

  int const flmn_size = so3_sampling_flmn_size(parameters);

  for (int i = 0; i < flmn_size; i += 1) {
    assert_float_equal(creal(flmn_orig[i]), creal(f_back[i]), state->tolerance);
    assert_float_equal(cimag(flmn_orig[i]), cimag(f_back[i]), state->tolerance);
  }

  free(flmn_orig);
  free(f_back);
  free(f_image);
}

void test_real_direct_vs_ssht(void **_state) {
  SO3TestState const *state = *(const SO3TestState **)_state;

  so3_parameters_t const *parameters = &state->params;

  complex double *flmn_orig = calloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L, sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(flmn_orig);
  gen_flmn_real(flmn_orig, parameters, state->seed);

  double *f_direct = calloc(
      (2 * parameters->L) * (parameters->L + 1) * (2 * parameters->N - 1),
      sizeof(double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_direct);
  so3_core_inverse_direct_real(f_direct, flmn_orig, parameters);

  double *f_ssht = calloc(
      (2 * parameters->L) * (parameters->L + 1) * (2 * parameters->N - 1),
      sizeof(double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_ssht);
  so3_core_inverse_via_ssht_real(f_ssht, flmn_orig, parameters);

  int const flmn_size = so3_sampling_flmn_size(parameters);

  for (int i = 0; i < flmn_size; i += 1)
    assert_float_equal(f_direct[i], f_ssht[i], state->tolerance);

  free(flmn_orig);
  free(f_direct);
  free(f_ssht);
}

void test_back_and_forth(void **_state) {
  SO3TestState const *state = *(const SO3TestState **)_state;
  so3_parameters_t const *parameters = &state->params;
  if (state->forward_complex == &so3_core_forward_via_ssht &&
      parameters->steerable != 0 &&
      (parameters->n_mode == SO3_N_MODE_ALL || parameters->n_mode == SO3_N_MODE_EVEN ||
       parameters->n_mode == SO3_N_MODE_L))
    skip();

  complex double *flmn_orig = calloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L, sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(flmn_orig);
  gen_flmn_complex(flmn_orig, parameters, state->seed);

  complex double *f_image = calloc(
      (2 * parameters->L) * (parameters->L + 1) * (2 * parameters->N - 1),
      sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_image);
  state->inverse_complex(f_image, flmn_orig, parameters);

  complex double *f_back = calloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L, sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_back);
  state->forward_complex(f_back, f_image, parameters);

  int const flmn_size = so3_sampling_flmn_size(parameters);

  for (int i = 0; i < flmn_size; i += 1) {
    assert_float_equal(creal(flmn_orig[i]), creal(f_back[i]), state->tolerance);
    assert_float_equal(cimag(flmn_orig[i]), cimag(f_back[i]), state->tolerance);
  }

  free(flmn_orig);
  free(f_back);
  free(f_image);
}

void test_direct_vs_ssht(void **_state) {
  SO3TestState const *state = *(const SO3TestState **)_state;
  so3_parameters_t const *parameters = &state->params;
  if (parameters->n_mode == SO3_N_MODE_L)
    skip();
  else if (parameters->steerable != 0 && parameters->n_mode == SO3_N_MODE_ALL)
    skip();
  else if (
      parameters->steerable != 0 && parameters->n_mode == SO3_N_MODE_EVEN &&
      parameters->n_order == SO3_N_ORDER_ZERO_FIRST)
    skip();
  else if (
      parameters->steerable != 0 && parameters->n_mode == SO3_N_MODE_ODD &&
      parameters->n_order == SO3_N_ORDER_ZERO_FIRST)
    skip();
  else if (
      parameters->steerable != 0 && parameters->n_mode == SO3_N_MODE_ODD &&
      parameters->n_order == SO3_N_ORDER_NEGATIVE_FIRST &&
      parameters->storage == SO3_STORAGE_COMPACT)
    skip();
  else if (
      parameters->steerable != 0 && parameters->n_mode == SO3_N_MODE_EVEN &&
      parameters->n_order == SO3_N_ORDER_NEGATIVE_FIRST &&
      parameters->storage == SO3_STORAGE_COMPACT)
    skip();

  complex double *flmn_orig = calloc(
      (2 * parameters->N - 1) * parameters->L * parameters->L, sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(flmn_orig);
  gen_flmn_real(flmn_orig, parameters, state->seed);

  complex double *f_direct = calloc(
      (2 * parameters->L) * (parameters->L + 1) * (2 * parameters->N - 1),
      sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_direct);
  so3_core_inverse_direct(f_direct, flmn_orig, parameters);

  complex double *f_ssht = calloc(
      (2 * parameters->L) * (parameters->L + 1) * (2 * parameters->N - 1),
      sizeof(complex double));
  SO3_ERROR_MEM_ALLOC_CHECK(f_ssht);
  so3_core_inverse_via_ssht(f_ssht, flmn_orig, parameters);

  int const flmn_size = so3_sampling_flmn_size(parameters);

  for (int i = 0; i < flmn_size; i += 1) {
    assert_float_equal(creal(f_direct[i]), creal(f_ssht[i]), state->tolerance);
    assert_float_equal(cimag(f_direct[i]), cimag(f_ssht[i]), state->tolerance);
  }

  free(flmn_orig);
  free(f_direct);
  free(f_ssht);
}

static SO3TestState *parametrization(
    char const *name,
    so3_sampling_t sampling,
    so3_n_order_t order,
    so3_n_mode_t mode,
    so3_storage_t storage,
    _Bool steerable,
    _Bool real) {
  SO3TestState *state = (SO3TestState *)calloc(1, sizeof(SO3TestState));
  SO3_ERROR_MEM_ALLOC_CHECK(state);
  state->seed = 10;
  state->tolerance = 1e-8;
  state->params.L0 = 0;
  state->params.L = 8;
  state->params.N = 8;
  state->params.verbosity = 0;
  state->params.dl_method = SSHT_DL_TRAPANI;
  state->params.sampling_scheme = sampling;
  state->params.n_order = order;
  state->params.n_mode = mode;
  state->params.storage = storage;
  state->params.steerable = steerable;
  state->params.reality = real;
  if (strcmp("ssht", name) == 0) {
    state->inverse_real = &so3_core_inverse_via_ssht_real;
    state->forward_real = &so3_core_forward_via_ssht_real;
    state->inverse_complex = &so3_core_inverse_via_ssht;
    state->forward_complex = &so3_core_forward_via_ssht;
  } else {
    state->inverse_real = &so3_core_inverse_direct_real;
    state->forward_real = &so3_core_forward_direct_real;
    state->inverse_complex = &so3_core_inverse_direct;
    state->forward_complex = &so3_core_forward_direct;
  }
  return state;
}

char const *name_of_test(
    char const *prefix,
    so3_sampling_t sampling,
    so3_n_order_t order,
    so3_n_mode_t mode,
    so3_storage_t storage,
    _Bool steerable,
    _Bool real) {
  char const *samplings[] = {"MW", "MWSS"};
  char const *orders[] = {"zero-first", "negative-first"};
  char const *modes[] = {"all", "even", "odd", "maximum", "l"};
  char const *storages[] = {"padded", "compact"};
  char const *steerables[] = {"non-steerable", "steerable"};
  char const *reals[] = {"complex", "real"};
  char *result = calloc(256, sizeof(char));
  sprintf(
      result,
      "%s: signal=%s, sampling=%s, order=%s, mode=%s, storage=%s, %s",
      prefix,
      reals[real],
      samplings[sampling],
      orders[order],
      modes[mode],
      storages[storage],
      steerables[steerable]);
  return result;
}

int main(void) {
  struct CMUnitTest tests[180];
  memset(tests, 0, sizeof(tests));

  int i = 0;
  const so3_sampling_t sampling = SO3_SAMPLING_MW;
  const so3_n_order_t order = SO3_N_ORDER_NEGATIVE_FIRST;
  for (so3_n_mode_t mode = 0; mode < SO3_N_MODE_SIZE; mode += 1)
    for (so3_storage_t storage = 0; storage < SO3_STORAGE_SIZE; storage += 1)
      for (int steerable = 0; steerable < 2; steerable += 1, i += 1) {
        assert(i + 2 < sizeof(tests) / sizeof(tests[0]));
        tests[i].name = name_of_test(
            "back_and_forth: ssht", sampling, order, mode, storage, steerable, 1);
        tests[i].initial_state =
            parametrization("ssht", sampling, order, mode, storage, steerable, 1);
        tests[i].test_func = &test_real_back_and_forth;
      }

  for (so3_n_mode_t mode = 0; mode < SO3_N_MODE_SIZE; mode += 1)
    for (so3_storage_t storage = 0; storage < SO3_STORAGE_SIZE; storage += 1)
      for (int steerable = 0; steerable < 2; steerable += 1, i += 1) {
        tests[i].name = name_of_test(
            "back_and_forth: direct", sampling, order, mode, storage, steerable, 1);
        tests[i].initial_state =
            parametrization("direct", sampling, order, mode, storage, steerable, 1);
        tests[i].test_func = &test_real_back_and_forth;
      }

  for (so3_n_mode_t mode = 0; mode < SO3_N_MODE_SIZE; mode += 1)
    for (so3_storage_t storage = 0; storage < SO3_STORAGE_SIZE; storage += 1)
      for (int steerable = 0; steerable < 2; steerable += 1, i += 1) {
        tests[i].name = name_of_test(
            "direct vs ssht", sampling, order, mode, storage, steerable, 1);
        tests[i].initial_state =
            parametrization("", sampling, order, mode, storage, steerable, 1);
        tests[i].test_func = &test_real_direct_vs_ssht;
      }

  for (so3_n_mode_t mode = 0; mode < SO3_N_MODE_SIZE; mode += 1)
    for (so3_storage_t storage = 0; storage < SO3_STORAGE_SIZE; storage += 1)
      for (so3_n_order_t order = 0; order < SO3_N_ORDER_SIZE; order += 1)
        for (int steerable = 0; steerable < 2; steerable += 1, i += 1) {
          assert(i + 2 < sizeof(tests) / sizeof(tests[0]));
          tests[i].name = name_of_test(
              "back_and_forth: ssht", sampling, order, mode, storage, steerable, 0);
          tests[i].initial_state =
              parametrization("ssht", sampling, order, mode, storage, steerable, 0);
          tests[i].test_func = &test_back_and_forth;
        }

  for (so3_n_mode_t mode = 0; mode < SO3_N_MODE_SIZE; mode += 1)
    for (so3_storage_t storage = 0; storage < SO3_STORAGE_SIZE; storage += 1)
      for (so3_n_order_t order = 0; order < SO3_N_ORDER_SIZE; order += 1)
        for (int steerable = 0; steerable < 2; steerable += 1, i += 1) {
          tests[i].name = name_of_test(
              "back_and_forth: direct", sampling, order, mode, storage, steerable, 0);
          tests[i].initial_state =
              parametrization("direct", sampling, order, mode, storage, steerable, 0);
          tests[i].test_func = &test_back_and_forth;
        }

  for (so3_n_mode_t mode = 0; mode < SO3_N_MODE_SIZE; mode += 1)
    for (so3_storage_t storage = 0; storage < SO3_STORAGE_SIZE; storage += 1)
      for (so3_n_order_t order = 0; order < SO3_N_ORDER_SIZE; order += 1)
        for (int steerable = 0; steerable < 2; steerable += 1, i += 1) {
          tests[i].name = name_of_test(
              "direct vs ssht", sampling, order, mode, storage, steerable, 0);
          tests[i].initial_state =
              parametrization("", sampling, order, mode, storage, steerable, 0);
          tests[i].test_func = &test_direct_vs_ssht;
        }

  int result = cmocka_run_group_tests(tests, NULL, NULL);

  struct CMUnitTest *deletee = tests;
  for (int i = 0; i < sizeof(tests) / sizeof(tests[0]); i += 1) {
    if (tests[i].name == NULL)
      break;
    free((void *)tests[i].name);
    free(tests[i].initial_state);
    deletee += 1;
  }
  return result;
}

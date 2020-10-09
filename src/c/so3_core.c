// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_core.c
 * Core algorithms to perform Wigner transform on the rotation group SO(§).
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <complex.h> // Must be before fftw3.h
#include <fftw3.h>
#include <math.h>
#include <ssht/ssht.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "so3/so3_error.h"
#include "so3/so3_sampling.h"
#include "so3/so3_types.h"

#define MIN(a, b) ((a < b) ? (a) : (b))
#define MAX(a, b) ((a > b) ? (a) : (b))

typedef void (*inverse_complex_ssht)(
    complex double *, const complex double *, int, int, int, ssht_dl_method_t, int);
typedef void (*inverse_real_ssht)(
    double *, const complex double *, int, int, ssht_dl_method_t, int);
typedef void (*forward_complex_ssht)(
    complex double *, const complex double *, int, int, int, ssht_dl_method_t, int);
typedef void (*forward_real_ssht)(
    complex double *, const double *, int, int, ssht_dl_method_t, int);

/*!
 * Compute inverse Wigner transform for a complex signal via SSHT.
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  flmn Harmonic coefficients.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameter_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_inverse_via_ssht_real
 *                        \endlink instead for real signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_inverse_via_ssht(
    complex double *f, const complex double *flmn, const so3_parameters_t *parameters) {

  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  // Iterator
  int n;
  // Intermediate results
  complex double *fn, *ftemp;
  // Stride for several arrays
  int fn_n_stride;
  // FFTW-related variables
  int fftw_rank, fftw_howmany;
  int fftw_idist, fftw_odist;
  int fftw_istride, fftw_ostride;
  int fftw_n;
  complex double *fftw_target;
  fftw_plan plan;

  inverse_complex_ssht ssht;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Support N_MODE_L
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_inverse_via_ssht with storage method %d...\n",
          SO3_PROMPT,
          storage);
  }

  switch (sampling) {
  case SO3_SAMPLING_MW:
    fn_n_stride = L * (2 * L - 1);
    ssht = ssht_core_mw_lb_inverse_sov_sym;
    break;
  case SO3_SAMPLING_MW_SS:
    fn_n_stride = (L + 1) * 2 * L;
    ssht = ssht_core_mw_lb_inverse_sov_sym_ss;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid sampling scheme.");
  }

  // Compute fn(a,b)

  if (steerable) {
    // For steerable signals, we need to supersample in n/gamma,
    // in order to create a symmetric sampling.
    fftw_n = 2 * N; // Each transform is over 2*N

    // We need to perform the FFT into a temporary buffer, because
    // the result will be twice as large as the output we need.
    ftemp = malloc(2 * N * fn_n_stride * sizeof *ftemp);
    SO3_ERROR_MEM_ALLOC_CHECK(ftemp);

    fftw_target = ftemp;
  } else {
    fftw_n = 2 * N - 1; // Each transform is over 2*N-1

    fftw_target = f;
  }

  fn = calloc(fftw_n * fn_n_stride, sizeof *fn);
  SO3_ERROR_MEM_ALLOC_CHECK(fn);

  // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
  // necessary but still good practice.
  fftw_rank = 1;              // We compute 1d transforms
  fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

  // We want to transform columns
  fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
  fftw_istride = fftw_ostride =
      fn_n_stride; // Distance between two elements of the same column

  plan = fftw_plan_many_dft(
      fftw_rank,
      &fftw_n,
      fftw_howmany,
      fn,
      NULL,
      fftw_istride,
      fftw_idist,
      fftw_target,
      NULL,
      fftw_ostride,
      fftw_odist,
      FFTW_BACKWARD,
      FFTW_ESTIMATE);

  for (n = -N + 1; n <= N - 1; ++n) {
    int ind, offset, i, el;
    int L0e = MAX(L0, abs(n)); // 'e' for 'effective'
    double factor;
    complex double *flm;

    if ((n_mode == SO3_N_MODE_EVEN && n % 2) ||
        (n_mode == SO3_N_MODE_ODD && !(n % 2)) ||
        (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N - 1)) {
      continue;
    }

    flm = malloc(L * L * sizeof *flm);
    SO3_ERROR_MEM_ALLOC_CHECK(flm);

    switch (storage) {
    case SO3_STORAGE_PADDED:
      so3_sampling_elmn2ind(&ind, 0, 0, n, parameters);
      memcpy(flm, flmn + ind, L * L * sizeof(complex double));
      break;
    case SO3_STORAGE_COMPACT:
      so3_sampling_elmn2ind(&ind, abs(n), -abs(n), n, parameters);
      memcpy(flm + n * n, flmn + ind, (L * L - n * n) * sizeof(complex double));
      for (i = 0; i < n * n; ++i)
        flm[i] = 0.0;
      break;
    default:
      SO3_ERROR_GENERIC("Invalid storage method.");
    }

    el = L0e;
    i = offset = el * el;
    for (; el < L; ++el) {
      factor = sqrt((double)(2 * el + 1) / (16. * pow(SO3_PI, 3.)));
      for (; i < offset + 2 * el + 1; ++i)
        flm[i] *= factor;

      offset = i;
    }

    // The conditional applies the spatial transform, so that we store
    // the results in n-order 0, 1, 2, -2, -1
    offset = (n < 0 ? n + fftw_n : n);

    (*ssht)(fn + offset * fn_n_stride, flm, L0e, L, -n, dl_method, verbosity);

    if (n % 2)
      for (i = 0; i < fn_n_stride; ++i)
        fn[offset * fn_n_stride + i] = -fn[offset * fn_n_stride + i];

    if (verbosity > 0)
      printf("\n");

    free(flm);
  }

  fftw_execute(plan);
  fftw_destroy_plan(plan);

  if (steerable) {
    memcpy(f, ftemp, N * fn_n_stride * sizeof(complex double));
    free(ftemp);
  }

  free(fn);

  if (verbosity > 0)
    printf("%sInverse transform computed!\n", SO3_PROMPT);
}

/*!
 * Compute forward Wigner transform for a complex signal via SSHT.
 *
 * \param[out] flmn Harmonic coefficients. If \link so3_parameters_t::n_mode n_mode
 *                  \endlink is different from \link SO3_N_MODE_ALL \endlink,
 *                  this array has to be nulled before being past to the function.
 * \param[in] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_forward_via_ssht_real
 *                        \endlink instead for real signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_forward_via_ssht(
    complex double *flmn, const complex double *f, const so3_parameters_t *parameters) {
  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  // Iterator
  int i, n;
  // Intermediate results
  complex double *ftemp, *fn;
  // Stride for several arrays
  int fn_n_stride;
  // FFTW-related variables
  int fftw_rank, fftw_howmany;
  int fftw_idist, fftw_odist;
  int fftw_istride, fftw_ostride;
  int fftw_n;
  fftw_plan plan;

  forward_complex_ssht ssht;

  // for precomputation
  double factor;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Support N_MODE_L
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing forward transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_forward_via_ssht with storage method %d...\n",
          SO3_PROMPT,
          storage);
  }

  switch (sampling) {
  case SO3_SAMPLING_MW:
    fn_n_stride = L * (2 * L - 1);
    ssht = ssht_core_mw_lb_forward_sov_conv_sym;
    break;
  case SO3_SAMPLING_MW_SS:
    fn_n_stride = (L + 1) * 2 * L;
    ssht = ssht_core_mw_lb_forward_sov_conv_sym_ss;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid sampling scheme.");
  }

  if (steerable) {
    int g, offset;

    fn = calloc((2 * N - 1) * fn_n_stride, sizeof *fn);
    SO3_ERROR_MEM_ALLOC_CHECK(fn);

    for (n = -N + 1; n < N; n += 2) {
      // The conditional applies the spatial transform, because the fn
      // are to be stored in n-order 0, 1, 2, -2, -1
      offset = (n < 0 ? n + 2 * N - 1 : n);

      for (g = 0; g < N; ++g) {
        double gamma = g * SO3_PI / N;
        for (i = 0; i < fn_n_stride; ++i) {
          double weight = 2 * SO3_PI / N;
          fn[offset * fn_n_stride + i] +=
              weight * f[g * fn_n_stride + i] * cexp(-I * n * gamma);
        }
      }
    }
  } else {
    // Make a copy of the input, because input is const
    // This could potentially be avoided by copying the input into fn and using an
    // in-place FFTW. The performance impact has to be profiled, though.
    ftemp = malloc((2 * N - 1) * fn_n_stride * sizeof *ftemp);
    SO3_ERROR_MEM_ALLOC_CHECK(ftemp);
    memcpy(ftemp, f, (2 * N - 1) * fn_n_stride * sizeof(complex double));

    fn = malloc((2 * N - 1) * fn_n_stride * sizeof *fn);
    SO3_ERROR_MEM_ALLOC_CHECK(fn);

    // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
    // necessary but still good practice.
    fftw_rank = 1;
    fftw_n = 2 * N - 1;
    fftw_howmany = fn_n_stride;
    fftw_idist = fftw_odist = 1;
    fftw_istride = fftw_ostride = fn_n_stride;

    plan = fftw_plan_many_dft(
        fftw_rank,
        &fftw_n,
        fftw_howmany,
        ftemp,
        NULL,
        fftw_istride,
        fftw_idist,
        fn,
        NULL,
        fftw_ostride,
        fftw_odist,
        FFTW_FORWARD,
        FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    free(ftemp);

    factor = 2 * SO3_PI / (double)(2 * N - 1);
    for (i = 0; i < (2 * N - 1) * fn_n_stride; ++i)
      fn[i] *= factor;
  }

  for (n = -N + 1; n <= N - 1; ++n) {
    int ind, offset, el, sign;
    int L0e = MAX(L0, abs(n)); // 'e' for 'effective'

    complex double *flm = NULL;

    if ((n_mode == SO3_N_MODE_EVEN && n % 2) ||
        (n_mode == SO3_N_MODE_ODD && !(n % 2)) ||
        (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N - 1)) {
      continue;
    }

    if (storage == SO3_STORAGE_COMPACT)
      flm = malloc(L * L * sizeof *flm);

    // The conditional applies the spatial transform, because the fn
    // are stored in n-order 0, 1, 2, -2, -1
    offset = (n < 0 ? n + 2 * N - 1 : n);

    complex double *flm_block;
    complex double *fn_block = fn + offset * fn_n_stride;

    el = L0e;
    switch (storage) {
    case SO3_STORAGE_PADDED:
      so3_sampling_elmn2ind(&ind, 0, 0, n, parameters);
      flm_block = flmn + ind;
      i = offset = el * el;
      break;
    case SO3_STORAGE_COMPACT:
      flm_block = flm;
      i = offset = el * el - n * n;
      break;
    default:
      SO3_ERROR_GENERIC("Invalid storage method.");
    }

    (*ssht)(flm_block, fn_block, L0e, L, -n, dl_method, verbosity);

    if (storage == SO3_STORAGE_COMPACT) {
      so3_sampling_elmn2ind(&ind, abs(n), -abs(n), n, parameters);
      memcpy(flmn + ind, flm + n * n, (L * L - n * n) * sizeof(complex double));
    }

    if (n % 2)
      sign = -1;
    else
      sign = 1;

    for (; el < L; ++el) {
      factor = sign * sqrt(4.0 * SO3_PI / (double)(2 * el + 1));
      for (; i < offset + 2 * el + 1; ++i)
        flmn[ind + i] *= factor;

      offset = i;
    }

    if (storage == SO3_STORAGE_COMPACT)
      free(flm);

    if (verbosity > 0)
      printf("\n");
  }

  free(fn);

  if (verbosity > 0)
    printf("%sForward transform computed!\n", SO3_PROMPT);
}

/*!
 * Compute inverse Wigner transform for a real signal via SSHT.
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in] flmn Harmonic coefficients for n >= 0. Note that for n = 0, these have to
 *                 respect the symmetry flm0* = (-1)^(m+n)*fl-m0, and hence fl00 has to
 * be real. \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_inverse_via_ssht
 *                        \endlink instead for complex signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_inverse_via_ssht_real(
    double *f, const complex double *flmn, const so3_parameters_t *parameters) {
  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  // Iterator
  int n;
  // Intermediate results
  complex double *fn, *flm;
  double *ftemp;
  // Stride for several arrays
  int fn_n_stride;
  // FFTW-related variables
  int fftw_rank, fftw_howmany;
  int fftw_idist, fftw_odist;
  int fftw_istride, fftw_ostride;
  int fftw_n;
  double *fftw_target;
  fftw_plan plan;

  inverse_complex_ssht complex_ssht;
  inverse_real_ssht real_ssht;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Support N_MODE_L
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_inverse_via_ssht_real with storage method "
          "%d...\n",
          SO3_PROMPT,
          storage);
  }

  switch (sampling) {
  case SO3_SAMPLING_MW:
    fn_n_stride = L * (2 * L - 1);
    complex_ssht = ssht_core_mw_lb_inverse_sov_sym;
    real_ssht = ssht_core_mw_lb_inverse_sov_sym_real;
    break;
  case SO3_SAMPLING_MW_SS:
    fn_n_stride = (L + 1) * 2 * L;
    complex_ssht = ssht_core_mw_lb_inverse_sov_sym_ss;
    real_ssht = ssht_core_mw_lb_inverse_sov_sym_ss_real;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid sampling scheme.");
  }

  // Compute fn(a,b)

  // For steerable signals, we need to supersample in n/gamma,
  // in order to create a symmetric sampling.
  // Each transform is over fftw_n samples (logically; physically, fn for negative n
  // will be omitted)
  // if (steerable)
  //{
  // For steerable signals, we need to supersample in n/gamma,
  // in order to create a symmetric sampling.
  //    fftw_n = 2*N;

  // We need to perform the FFT into a temporary buffer, because
  // the result will be twice as large as the output we need.
  //    ftemp = malloc(2*N*fn_n_stride * sizeof *ftemp);
  //    SO3_ERROR_MEM_ALLOC_CHECK(ftemp);

  //    fftw_target = ftemp;
  //}
  // else
  //{
  fftw_n = 2 * N - 1;

  fftw_target = f;
  //}

  // Only need to store for non-negative n
  fn = calloc((fftw_n / 2 + 1) * fn_n_stride, sizeof *fn);
  SO3_ERROR_MEM_ALLOC_CHECK(fn);

  // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
  // necessary but still good practice.
  fftw_rank = 1;              // We compute 1d transforms
  fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

  // We want to transform columns
  fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
  fftw_istride = fftw_ostride =
      fn_n_stride; // Distance between two elements of the same column

  plan = fftw_plan_many_dft_c2r(
      fftw_rank,
      &fftw_n,
      fftw_howmany,
      fn,
      NULL,
      fftw_istride,
      fftw_idist,
      fftw_target,
      NULL,
      fftw_ostride,
      fftw_odist,
      FFTW_ESTIMATE);

  flm = malloc(L * L * sizeof *flm);
  SO3_ERROR_MEM_ALLOC_CHECK(flm);

  for (n = 0; n <= N - 1; ++n) {
    int ind, offset, i, el;
    int L0e = MAX(L0, abs(n)); // 'e' for 'effective'
    double factor;

    if ((n_mode == SO3_N_MODE_EVEN && n % 2) ||
        (n_mode == SO3_N_MODE_ODD && !(n % 2)) ||
        (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N - 1)) {
      continue;
    }

    switch (storage) {
    case SO3_STORAGE_PADDED:
      so3_sampling_elmn2ind_real(
          &ind, 0, 0, n, parameters); // L, N, SO3_STORE_NEG_FIRST_PAD);
      memcpy(flm, flmn + ind, L * L * sizeof(complex double));
      break;
    case SO3_STORAGE_COMPACT:
      so3_sampling_elmn2ind_real(
          &ind, n, -n, n, parameters); // L, N, SO3_STORE_NEG_FIRST_COMPACT);
      memcpy(flm + n * n, flmn + ind, (L * L - n * n) * sizeof(complex double));
      for (i = 0; i < n * n; ++i)
        flm[i] = 0.0;
      break;
    default:
      SO3_ERROR_GENERIC("Invalid storage method.");
    }

    el = L0e;
    i = offset = el * el;
    for (; el < L; ++el) {
      factor = sqrt((double)(2 * el + 1) / (16. * pow(SO3_PI, 3.)));
      for (; i < offset + 2 * el + 1; ++i)
        flm[i] *= factor;

      offset = i;
    }

    if (N > 1 || n) {
      (*complex_ssht)(fn + n * fn_n_stride, flm, L0e, L, -n, dl_method, verbosity);
    } else {
      double *fn_r;

      // Create an array of real doubles for n = 0
      fn_r = calloc(fn_n_stride, sizeof *fn_r);
      SO3_ERROR_MEM_ALLOC_CHECK(fn_r);

      (*real_ssht)(fn_r, flm, L0e, L, dl_method, verbosity);

      for (i = 0; i < fn_n_stride; ++i)
        fn[i] = fn_r[i];
    }

    if (n % 2)
      for (i = 0; i < fn_n_stride; ++i)
        fn[n * fn_n_stride + i] = -fn[n * fn_n_stride + i];

    if (verbosity > 0)
      printf("\n");
  }

  free(flm);

  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // if (steerable)
  //{
  //    memcpy(f, ftemp, N*fn_n_stride * sizeof *f);
  //    free(ftemp);
  //}

  free(fn);

  if (verbosity > 0)
    printf("%sInverse transform computed!\n", SO3_PROMPT);
}

/*!
 * Compute forward Wigner transform for a real signal via SSHT.
 *
 * \param[out] flmn Harmonic coefficients. If \link so3_parameters_t::n_mode n_mode
 *                  \endlink is different from \link SO3_N_MODE_ALL \endlink,
 *                  this array has to be nulled before being past to the function.
 * \param[in] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality \endlink flag
 *                        is ignored. Use \link so3_core_forward_via_ssht
 *                        \endlink instead for complex signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_forward_via_ssht_real(
    complex double *flmn, const double *f, const so3_parameters_t *parameters) {
  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  // Iterator
  int i, n;
  // Intermediate results
  double *ftemp;
  complex double *flm = NULL, *fn;
  // Stride for several arrays
  int fn_n_stride;
  // FFTW-related variables
  int fftw_rank, fftw_howmany;
  int fftw_idist, fftw_odist;
  int fftw_istride, fftw_ostride;
  int fftw_n;
  fftw_plan plan;

  forward_complex_ssht complex_ssht;
  forward_real_ssht real_ssht;

  // for precomputation
  double factor;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Support N_MODE_L
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  steerable = parameters->steerable;
  verbosity = parameters->verbosity;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing forward transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_forward_via_ssht_real with storage method "
          "%d...\n",
          SO3_PROMPT,
          storage);
  }

  switch (sampling) {
  case SO3_SAMPLING_MW:
    fn_n_stride = L * (2 * L - 1);
    complex_ssht = ssht_core_mw_lb_forward_sov_conv_sym;
    real_ssht = ssht_core_mw_lb_forward_sov_conv_sym_real;
    break;
  case SO3_SAMPLING_MW_SS:
    fn_n_stride = (L + 1) * 2 * L;
    complex_ssht = ssht_core_mw_lb_forward_sov_conv_sym_ss;
    real_ssht = ssht_core_mw_lb_forward_sov_conv_sym_ss_real;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid sampling scheme.");
  }

  // if (steerable)
  //{
  //    int g, offset;

  //    fn = calloc((2*N-1)*fn_n_stride, sizeof *fn);
  //    SO3_ERROR_MEM_ALLOC_CHECK(fn);

  //    for (n = -N+1; n < N; n+=2)
  //    {
  // The conditional applies the spatial transform, because the fn
  // are to be stored in n-order 0, 1, 2, -2, -1
  //        offset = (n < 0 ? n + 2*N-1 : n);

  //        for (g = 0; g < N; ++g)
  //        {
  //            double gamma = g * SO3_PI / N;
  //            for (i = 0; i < fn_n_stride; ++i)
  //            {
  //                double weight = 2*SO3_PI/N;
  //                fn[offset * fn_n_stride + i] += weight*f[g * fn_n_stride +
  //                i]*cexp(-I*n*gamma);
  //            }
  //        }
  //    }
  //}
  // else
  //{
  // Make a copy of the input, because input is const
  // This could potentially be avoided by copying the input into fn and using an
  // in-place FFTW. The performance impact has to be profiled, though.
  ftemp = malloc((2 * N - 1) * fn_n_stride * sizeof *ftemp);
  SO3_ERROR_MEM_ALLOC_CHECK(ftemp);
  memcpy(ftemp, f, (2 * N - 1) * fn_n_stride * sizeof(double));

  fn = malloc(N * fn_n_stride * sizeof *fn);
  SO3_ERROR_MEM_ALLOC_CHECK(fn);
  // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
  // necessary but still good practice.
  fftw_rank = 1;      // We compute 1d transforms
  fftw_n = 2 * N - 1; // Each transform is over 2*N-1 (logically; physically, fn for
                      // negative n will be omitted)
  fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

  // We want to transform columns
  fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
  fftw_istride = fftw_ostride =
      fn_n_stride; // Distance between two elements of the same column

  plan = fftw_plan_many_dft_r2c(
      fftw_rank,
      &fftw_n,
      fftw_howmany,
      ftemp,
      NULL,
      fftw_istride,
      fftw_idist,
      fn,
      NULL,
      fftw_ostride,
      fftw_odist,
      FFTW_ESTIMATE);

  fftw_execute(plan);
  fftw_destroy_plan(plan);

  free(ftemp);

  factor = 2 * SO3_PI / (double)(2 * N - 1);
  for (i = 0; i < N * fn_n_stride; ++i)
    fn[i] *= factor;
  //}

  if (storage == SO3_STORAGE_COMPACT)
    flm = malloc(L * L * sizeof *flm);

  for (n = 0; n <= N - 1; ++n) {
    int ind, offset, el, sign;
    int L0e = MAX(L0, abs(n)); // 'e' for 'effective'

    complex double *flm_block;

    if ((n_mode == SO3_N_MODE_EVEN && n % 2) ||
        (n_mode == SO3_N_MODE_ODD && !(n % 2)) ||
        (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N - 1)) {
      continue;
    }

    el = L0e;
    i = offset = el * el;
    switch (storage) {
    case SO3_STORAGE_PADDED:
      so3_sampling_elmn2ind_real(&ind, 0, 0, n, parameters);
      flm_block = flmn + ind;
      break;
    case SO3_STORAGE_COMPACT:
      flm_block = flm;

      offset -= n * n;
      i = offset;
      break;
    default:
      SO3_ERROR_GENERIC("Invalid storage method.");
    }

    if (N > 1) {
      (*complex_ssht)(
          flm_block, fn + n * fn_n_stride, L0e, L, -n, dl_method, verbosity);
    } else {
      // Now we know n = 0 in which case the reality conditions
      // for SO3 and SSHT coincide.
      int j;
      double *fn_r;

      // Create an array of real doubles for n = 0
      fn_r = malloc(fn_n_stride * sizeof *fn_r);
      SO3_ERROR_MEM_ALLOC_CHECK(fn_r);
      for (j = 0; j < fn_n_stride; ++j)
        fn_r[j] = creal(fn[j]);

      // Now use real SSHT transforms
      (*real_ssht)(flm_block, fn_r, L0e, L, dl_method, verbosity);

      free(fn_r);
    }

    if (storage == SO3_STORAGE_COMPACT) {
      so3_sampling_elmn2ind_real(&ind, n, -n, n, parameters);
      memcpy(flmn + ind, flm + n * n, (L * L - n * n) * sizeof(complex double));
    }

    if (n % 2)
      sign = -1;
    else
      sign = 1;

    for (; el < L; ++el) {
      factor = sign * sqrt(4.0 * SO3_PI / (double)(2 * el + 1));
      for (; i < offset + 2 * el + 1; ++i)
        flmn[ind + i] *= factor;

      offset = i;
    }

    if (verbosity > 0)
      printf("\n");
  }

  if (storage == SO3_STORAGE_COMPACT)
    free(flm);

  free(fn);

  if (verbosity > 0)
    printf("%sForward transform computed!\n", SO3_PROMPT);
}

/*!
 * Compute inverse Wigner transform for a complex signal directly (without using
 * SSHT).
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  flmn Harmonic coefficients.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameter_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_inverse_via_ssht_real
 *                        \endlink instead for real signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_inverse_direct(
    complex double *f, const complex double *flmn, const so3_parameters_t *parameters) {

  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Add optimisations for all n-modes.
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_inverse_direct with storage method %d...\n",
          SO3_PROMPT,
          storage);
  }

  // Iterators
  int el, m, n, mm; // mm for m'

  // Allocate memory.
  double *sqrt_tbl = calloc(2 * (L - 1) + 2, sizeof(*sqrt_tbl));
  SO3_ERROR_MEM_ALLOC_CHECK(sqrt_tbl);
  double *signs = calloc(L + 1, sizeof(*signs));
  SO3_ERROR_MEM_ALLOC_CHECK(signs);
  complex double *exps = calloc(4, sizeof(*exps));
  SO3_ERROR_MEM_ALLOC_CHECK(exps);

  // Perform precomputations.
  for (el = 0; el <= 2 * L - 1; ++el)
    sqrt_tbl[el] = sqrt((double)el);
  for (m = 0; m <= L - 1; m += 2) {
    signs[m] = 1.0;
    signs[m + 1] = -1.0;
  }
  int i;
  for (i = 0; i < 4; ++i)
    exps[i] = cexp(I * SO3_PION2 * i);

  // Compute Fmnm'
  // TODO: Currently m is fastest-varying, then n, then m'.
  // Should this order be changed to m-m'-n?
  complex double *Fmnm = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*Fmnm));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm);
  int m_offset = L - 1;
  int m_stride = 2 * L - 1;
  int n_offset = N - 1;
  int n_stride = 2 * N - 1;
  int mm_offset = L - 1;
  int mm_stride = 2 * L - 1;

  int n_start, n_stop, n_inc;

  double *dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SO3_ERROR_MEM_ALLOC_CHECK(dl);
  double *dl8 = NULL;
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SO3_ERROR_MEM_ALLOC_CHECK(dl8);
  }
  int dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  int dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);

  complex double *mn_factors = calloc((2 * L - 1) * (2 * N - 1), sizeof *mn_factors);
  SO3_ERROR_MEM_ALLOC_CHECK(mn_factors);

  // TODO: SSHT starts this loop from MAX(L0, abs(spin)).
  // Can we use a similar optimisation? el can probably
  // be limited by n, but then we'd need to switch the
  // loop order, which means we'd have to recompute the
  // Wigner plane for each n. That seems wrong?
  for (el = L0; el <= L - 1; ++el) {
    int eltmp;
    // Compute Wigner plane.
    switch (dl_method) {
    case SSHT_DL_RISBO:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_beta_risbo_eighth_table(
              dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, eltmp, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      } else {
        ssht_dl_beta_risbo_eighth_table(
            dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, el, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      }
      break;

    case SSHT_DL_TRAPANI:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, eltmp, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      } else {
        ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, el, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      }
      break;

    default:
      SO3_ERROR_GENERIC("Invalid dl method");
    }

    // Compute Fmnm' contribution for current el.

    // Factor which depends only on el.
    double elfactor = (2.0 * el + 1.0) / (8.0 * SO3_PI * SO3_PI);

    switch (n_mode) {
    case SO3_N_MODE_ALL:
      n_start = MAX(-N + 1, -el);
      n_stop = MIN(N - 1, el);
      n_inc = 1;
      break;
    case SO3_N_MODE_EVEN:
      n_start = MAX(-N + 1, -el);
      n_start += (-n_start) % 2;
      n_stop = MIN(N - 1, el);
      n_stop -= n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_ODD:
      n_start = MAX(-N + 1, -el);
      n_start += 1 + n_start % 2;
      n_stop = MIN(N - 1, el);
      n_stop -= 1 - n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_MAXIMUM:
      if (el < N - 1)
        continue;
      n_start = -N + 1;
      n_stop = N - 1;
      n_inc = MAX(1, 2 * N - 2);
      break;
    case SO3_N_MODE_L:
      if (el >= N)
        continue;
      n_start = -el;
      n_stop = el;
      n_inc = MAX(1, 2 * el);
      break;
    default:
      SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    // Factors which do not depend on m'.
    for (n = n_start; n <= n_stop; n += n_inc)
      for (m = -el; m <= el; ++m) {
        int ind;
        so3_sampling_elmn2ind(&ind, el, m, n, parameters);
        int mod = ((n - m) % 4 + 4) % 4;
        mn_factors[m + m_offset + m_stride * (n + n_offset)] = flmn[ind] * exps[mod];
      }

    for (mm = 0; mm <= el; ++mm) {
      // These signs are needed for the symmetry relations of
      // Wigner symbols.
      double elmmsign = signs[el] * signs[mm];

      // TODO: If the conditional for elnsign is a bottleneck
      // this loop can be split up just like the inner loop.
      for (n = n_start; n <= n_stop; n += n_inc) {
        double elnsign = n >= 0 ? 1.0 : elmmsign;
        // Factor which does not depend on m.
        double elnmm_factor =
            elfactor * elnsign * dl[abs(n) + dl_offset + mm * dl_stride];
        for (m = -el; m < 0; ++m)
          Fmnm
              [m + m_offset +
               m_stride * (n + n_offset + n_stride * (mm + mm_offset))] +=
              elnmm_factor * mn_factors[m + m_offset + m_stride * (n + n_offset)] *
              elmmsign * dl[-m + dl_offset + mm * dl_stride];
        for (m = 0; m <= el; ++m)
          Fmnm
              [m + m_offset +
               m_stride * (n + n_offset + n_stride * (mm + mm_offset))] +=
              elnmm_factor * mn_factors[m + m_offset + m_stride * (n + n_offset)] *
              dl[m + dl_offset + mm * dl_stride];
      }
    }
  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  switch (n_mode) {
  case SO3_N_MODE_ALL:
  case SO3_N_MODE_L:
    n_start = -N + 1;
    n_stop = N - 1;
    n_inc = 1;
    break;
  case SO3_N_MODE_EVEN:
    n_start = ((N - 1) % 2 == 0) ? -N + 1 : -N + 2;
    n_stop = ((N - 1) % 2 == 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_ODD:
    n_start = ((N - 1) % 2 != 0) ? -N + 1 : -N + 2;
    n_stop = ((N - 1) % 2 != 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_MAXIMUM:
    n_start = -N + 1;
    n_stop = N - 1;
    n_inc = MAX(1, 2 * N - 2);
    break;
  default:
    SO3_ERROR_GENERIC("Invalid n-mode.");
  }

  // Use symmetry to compute Fmnm' for negative m'.
  for (mm = -L + 1; mm < 0; ++mm)
    for (n = n_start; n <= n_stop; n += n_inc)
      for (m = -L + 1; m <= L - 1; ++m)
        Fmnm[m + m_offset + m_stride * (n + n_offset + n_stride * (mm + mm_offset))] =
            signs[abs(m + n) % 2] *
            Fmnm
                [m + m_offset +
                 m_stride * (n + n_offset + n_stride * (-mm + mm_offset))];

  // Apply phase modulation to account for sampling offset.
  for (mm = -L + 1; mm <= L - 1; ++mm) {
    complex double mmfactor = cexp(I * mm * SO3_PI / (2.0 * L - 1.0));
    for (n = n_start; n <= n_stop; n += n_inc)
      for (m = -L + 1; m <= L - 1; ++m)
        Fmnm[m + m_offset + m_stride * (n + n_offset + n_stride * (mm + mm_offset))] *=
            mmfactor;
  }

  // Allocate space for function values.
  complex double *fext = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*fext));
  SO3_ERROR_MEM_ALLOC_CHECK(fext);

  // Set up plan before initialising array.
  fftw_plan plan = fftw_plan_dft_3d(
      2 * N - 1, 2 * L - 1, 2 * L - 1, fext, fext, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Apply spatial shift.
  for (mm = -L + 1; mm <= L - 1; ++mm) {
    int mm_shift = mm < 0 ? 2 * L - 1 : 0;
    for (n = n_start; n <= n_stop; n += n_inc) {
      int n_shift = n < 0 ? 2 * N - 1 : 0;
      for (m = -L + 1; m <= L - 1; ++m) {
        int m_shift = m < 0 ? 2 * L - 1 : 0;
        fext[m + m_shift + m_stride * (mm + mm_shift + mm_stride * (n + n_shift))] =
            Fmnm
                [m + m_offset +
                 m_stride * (n + n_offset + n_stride * (mm + mm_offset))];
      }
    }
  }

  // Free Fmnm' memory.
  free(Fmnm);

  // Perform 3D FFT.
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // Extract f from the extended torus.
  int a, b, g;
  int a_stride = 2 * L - 1;
  int b_ext_stride = 2 * L - 1;
  int b_stride = L;
  for (g = 0; g < 2 * N - 1; ++g)
    for (b = 0; b < L; ++b)
      for (a = 0; a < 2 * L - 1; ++a)
        f[a + a_stride * (b + b_stride * (g))] =
            fext[a + a_stride * (b + b_ext_stride * (g))];

  // Free fext memory.
  free(fext);

  if (verbosity > 0)
    printf("%sInverse transform computed!\n", SO3_PROMPT);

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(mn_factors);
}

/*!
 * Compute forward Wigner transform for a complex signal directly (without using
 * SSHT).
 *
 * \param[out] flmn Harmonic coefficients. If \link so3_parameters_t::n_mode n_mode
 *                  \endlink is different from \link SO3_N_MODE_ALL \endlink,
 *                  this array has to be nulled before being passed to the function.
 * \param[in] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_forward_via_ssht_real
 *                        \endlink instead for real signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_forward_direct(
    complex double *flmn, const complex double *f, const so3_parameters_t *parameters) {
  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Add optimisations for all n-modes.
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing forward transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_forward_direct with storage method %d...\n",
          SO3_PROMPT,
          storage);
  }

  int m_stride = 2 * L - 1;
  int m_offset = L - 1;
  // unused: int n_stride = 2*N-1;
  int n_offset = N - 1;
  int mm_stride = 2 * L - 1;
  int mm_offset = L - 1;
  int a_stride = 2 * L - 1;
  int b_stride = L;
  int bext_stride = 2 * L - 1;
  // unused: int g_stride = 2*N-1;

  int n_start, n_stop, n_inc;

  switch (n_mode) {
  case SO3_N_MODE_ALL:
  case SO3_N_MODE_L:
    n_start = -N + 1;
    n_stop = N - 1;
    n_inc = 1;
    break;
  case SO3_N_MODE_EVEN:
    n_start = ((N - 1) % 2 == 0) ? -N + 1 : -N + 2;
    n_stop = ((N - 1) % 2 == 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_ODD:
    n_start = ((N - 1) % 2 != 0) ? -N + 1 : -N + 2;
    n_stop = ((N - 1) % 2 != 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_MAXIMUM:
    n_start = -N + 1;
    n_stop = N - 1;
    n_inc = MAX(1, 2 * N - 2);
    break;
  default:
    SO3_ERROR_GENERIC("Invalid n-mode.");
  }

  double *sqrt_tbl = calloc(2 * (L - 1) + 2, sizeof(*sqrt_tbl));
  SO3_ERROR_MEM_ALLOC_CHECK(sqrt_tbl);
  double *signs = calloc(L + 1, sizeof(*signs));
  SO3_ERROR_MEM_ALLOC_CHECK(signs);
  complex double *exps = calloc(4, sizeof(*exps));
  SO3_ERROR_MEM_ALLOC_CHECK(exps);
  complex double *expsmm = calloc(2 * L - 1, sizeof(*expsmm));
  SO3_ERROR_MEM_ALLOC_CHECK(expsmm);

  int el, m, n, mm; // mm is for m'
  // Perform precomputations.
  for (el = 0; el <= 2 * (L - 1) + 1; ++el)
    sqrt_tbl[el] = sqrt((double)el);
  for (m = 0; m <= L - 1; m += 2) {
    signs[m] = 1.0;
    signs[m + 1] = -1.0;
  }
  int i;
  for (i = 0; i < 4; ++i)
    exps[i] = cexp(I * SO3_PION2 * i);
  for (mm = -L + 1; mm <= L - 1; ++mm)
    expsmm[mm + mm_offset] = cexp(-I * mm * SSHT_PI / (2.0 * L - 1.0));

  double norm_factor = 1.0 / (2.0 * L - 1.0) / (2.0 * N - 1.0);

  // Compute Fourier transform over alpha and gamma, i.e. compute Fmn(b).
  complex double *Fmnb = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*Fmnb));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnb);
  complex double *inout = calloc((2 * L - 1) * (2 * N - 1), sizeof(*inout));
  SO3_ERROR_MEM_ALLOC_CHECK(inout);
  fftw_plan plan =
      fftw_plan_dft_2d(2 * N - 1, 2 * L - 1, inout, inout, FFTW_FORWARD, FFTW_ESTIMATE);

  int b, g;
  for (b = 0; b < L; ++b) {
    // TODO: This memcpy loop could probably be avoided by using
    // a more elaborate FFTW plan which performs the FFT directly
    // over the 1st and 3rd dimensions of f.
    for (g = 0; g < 2 * N - 1; ++g)
      memcpy(
          inout + g * a_stride,
          f + 0 + a_stride * (b + b_stride * (g)),
          a_stride * sizeof(*f));
    fftw_execute_dft(plan, inout, inout);

    // Apply spatial shift and normalisation factor
    for (n = n_start; n <= n_stop; n += n_inc) {
      int n_shift = n < 0 ? 2 * N - 1 : 0;
      for (m = -L + 1; m <= L - 1; ++m) {
        int m_shift = m < 0 ? 2 * L - 1 : 0;
        Fmnb[b + bext_stride * (m + m_offset + m_stride * (n + n_offset))] =
            inout[m + m_shift + m_stride * (n + n_shift)] * norm_factor;
      }
    }
  }
  fftw_destroy_plan(plan);

  // Extend Fmnb periodically.
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m) {
      int signmn = signs[abs(m + n) % 2];
      for (b = L; b < 2 * L - 1; ++b)
        Fmnb[b + bext_stride * (m + m_offset + m_stride * (n + n_offset))] =
            signmn * Fmnb
                         [(2 * L - 2 - b) +
                          bext_stride * (m + m_offset + m_stride * (n + n_offset))];
    }

  // Compute Fourier transform over beta, i.e. compute Fmnm'.
  complex double *Fmnm = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*Fmnm));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm);

  plan = fftw_plan_dft_1d(2 * L - 1, inout, inout, FFTW_FORWARD, FFTW_ESTIMATE);
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m) {
      memcpy(
          inout,
          Fmnb + 0 + bext_stride * (m + m_offset + m_stride * (n + n_offset)),
          bext_stride * sizeof(*Fmnb));
      fftw_execute_dft(plan, inout, inout);

      // Apply spatial shift and normalisation factor
      for (mm = -L + 1; mm <= L - 1; ++mm) {
        int mm_shift = mm < 0 ? 2 * L - 1 : 0;
        Fmnm[mm + mm_offset + mm_stride * (m + m_offset + m_stride * (n + n_offset))] =
            inout[mm + mm_shift] / (2.0 * L - 1.0);
      }
    }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m)
      for (mm = -L + 1; mm <= L - 1; ++mm)
        Fmnm[mm + mm_offset + mm_stride * (m + m_offset + m_stride * (n + n_offset))] *=
            expsmm[mm + mm_offset];

  // Compute weights.
  complex double *w = calloc(4 * L - 3, sizeof(*w));
  SO3_ERROR_MEM_ALLOC_CHECK(w);
  int w_offset = 2 * (L - 1);
  for (mm = -2 * (L - 1); mm <= 2 * (L - 1); ++mm)
    w[mm + w_offset] = so3_sampling_weight(parameters, mm);

  // Compute IFFT of w to give wr.
  complex double *wr = calloc(4 * L - 3, sizeof(*w));
  SO3_ERROR_MEM_ALLOC_CHECK(wr);
  inout = calloc(4 * L - 3, sizeof(*inout));
  SO3_ERROR_MEM_ALLOC_CHECK(inout);
  fftw_plan plan_bwd =
      fftw_plan_dft_1d(4 * L - 3, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  fftw_plan plan_fwd =
      fftw_plan_dft_1d(4 * L - 3, inout, inout, FFTW_FORWARD, FFTW_MEASURE);

  // Apply spatial shift.
  for (mm = 1; mm <= 2 * L - 2; ++mm)
    inout[mm + w_offset] = w[mm - 2 * (L - 1) - 1 + w_offset];
  for (mm = -2 * (L - 1); mm <= 0; ++mm)
    inout[mm + w_offset] = w[mm + 2 * (L - 1) + w_offset];

  fftw_execute_dft(plan_bwd, inout, inout);

  // Apply spatial shift.
  for (mm = 0; mm <= 2 * L - 2; ++mm)
    wr[mm + w_offset] = inout[mm - 2 * (L - 1) + w_offset];
  for (mm = -2 * (L - 1); mm <= -1; ++mm)
    wr[mm + w_offset] = inout[mm + 2 * (L - 1) + 1 + w_offset];

  // Compute Gmnm' by convolution implemented as product in real space.
  complex double *Fmnm_pad = calloc(4 * L - 3, sizeof(*Fmnm_pad));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm_pad);
  complex double *Gmnm = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*Gmnm));
  SO3_ERROR_MEM_ALLOC_CHECK(Gmnm);
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m) {

      // Zero-pad Fmnm'.
      for (mm = -2 * (L - 1); mm <= -L; ++mm)
        Fmnm_pad[mm + w_offset] = 0.0;
      for (mm = L; mm <= 2 * (L - 1); ++mm)
        Fmnm_pad[mm + w_offset] = 0.0;
      for (mm = -(L - 1); mm <= L - 1; ++mm)
        Fmnm_pad[mm + w_offset] = Fmnm
            [mm + mm_offset + mm_stride * (m + m_offset + m_stride * (n + n_offset))];

      // Apply spatial shift.
      for (mm = 1; mm <= 2 * L - 2; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm - 2 * (L - 1) - 1 + w_offset];
      for (mm = -2 * (L - 1); mm <= 0; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm + 2 * (L - 1) + w_offset];
      // Compute IFFT of Fmnm'.
      fftw_execute_dft(plan_bwd, inout, inout);
      // Apply spatial shift.
      for (mm = 0; mm <= 2 * L - 2; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm - 2 * (L - 1) + w_offset];
      for (mm = -2 * (L - 1); mm <= -1; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm + 2 * (L - 1) + 1 + w_offset];

      // Compute product of Fmnm' and weight in real space.
      int r;
      for (r = -2 * (L - 1); r <= 2 * (L - 1); ++r)
        Fmnm_pad[r + w_offset] *= wr[r + w_offset];

      // Apply spatial shift.
      for (mm = 1; mm <= 2 * L - 2; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm - 2 * (L - 1) - 1 + w_offset];
      for (mm = -2 * (L - 1); mm <= 0; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm + 2 * (L - 1) + w_offset];
      // Compute Gmnm' by FFT.
      fftw_execute_dft(plan_fwd, inout, inout);
      // Apply spatial shift.
      for (mm = 0; mm <= 2 * L - 2; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm - 2 * (L - 1) + w_offset];
      for (mm = -2 * (L - 1); mm <= -1; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm + 2 * (L - 1) + 1 + w_offset];

      // Extract section of Gmnm' of interest.
      for (mm = -(L - 1); mm <= L - 1; ++mm)
        Gmnm[m + m_offset + m_stride * (mm + mm_offset + mm_stride * (n + n_offset))] =
            Fmnm_pad[mm + w_offset] * 4.0 * SSHT_PI * SSHT_PI / (4.0 * L - 3.0);
    }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flmn.
  double *dl, *dl8 = NULL;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SO3_ERROR_MEM_ALLOC_CHECK(dl);
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SO3_ERROR_MEM_ALLOC_CHECK(dl8);
  }
  int dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  int dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  for (n = -N + 1; n <= N - 1; ++n)
    for (el = abs(n); el < L; ++el)
      for (m = -el; m <= el; ++m) {
        int ind;
        so3_sampling_elmn2ind(&ind, el, m, n, parameters);
        flmn[ind] = 0.0;
      }

  for (el = L0; el < L; ++el) {
    int eltmp;

    // Compute Wigner plane.
    switch (dl_method) {
    case SSHT_DL_RISBO:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_beta_risbo_eighth_table(
              dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, eltmp, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      } else {
        ssht_dl_beta_risbo_eighth_table(
            dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, el, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      }
      break;

    case SSHT_DL_TRAPANI:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, eltmp, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      } else {
        ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, el, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      }
      break;

    default:
      SO3_ERROR_GENERIC("Invalid dl method");
    }

    // Compute flmn for current el.

    switch (n_mode) {
    case SO3_N_MODE_ALL:
      n_start = MAX(-N + 1, -el);
      n_stop = MIN(N - 1, el);
      n_inc = 1;
      break;
    case SO3_N_MODE_EVEN:
      n_start = MAX(-N + 1, -el);
      n_start += (-n_start) % 2;
      n_stop = MIN(N - 1, el);
      n_stop -= n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_ODD:
      n_start = MAX(-N + 1, -el);
      n_start += 1 + n_start % 2;
      n_stop = MIN(N - 1, el);
      n_stop -= 1 - n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_MAXIMUM:
      if (el < N - 1)
        continue;
      n_start = -N + 1;
      n_stop = N - 1;
      n_inc = MAX(1, 2 * N - 2);
      break;
    case SO3_N_MODE_L:
      if (el >= N)
        continue;
      n_start = -el;
      n_stop = el;
      n_inc = MAX(1, 2 * el);
      break;
    default:
      SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    // TODO: Pull out a few multiplications into precomputations
    // or split up loops to avoid conditionals to check signs.
    for (mm = -el; mm <= el; ++mm) {
      // These signs are needed for the symmetry relations of
      // Wigner symbols.
      double elmmsign = signs[el] * signs[abs(mm)];

      for (n = n_start; n <= n_stop; n += n_inc) {
        double mmsign = mm >= 0 ? 1.0 : signs[el] * signs[abs(n)];
        double elnsign = n >= 0 ? 1.0 : elmmsign;

        // Factor which does not depend on m.
        double elnmm_factor =
            mmsign * elnsign * dl[abs(n) + dl_offset + abs(mm) * dl_stride];

        for (m = -el; m <= el; ++m) {
          mmsign = mm >= 0 ? 1.0 : signs[el] * signs[abs(m)];
          double elmsign = m >= 0 ? 1.0 : elmmsign;
          int ind;
          so3_sampling_elmn2ind(&ind, el, m, n, parameters);
          int mod = ((m - n) % 4 + 4) % 4;
          flmn[ind] += exps[mod] * elnmm_factor * mmsign * elmsign *
                       dl[abs(m) + dl_offset + abs(mm) * dl_stride] *
                       Gmnm
                           [m + m_offset +
                            m_stride * (mm + mm_offset + mm_stride * (n + n_offset))];
        }
      }
    }
  }

  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmnb);
  free(Fmnm);
  free(inout);
  free(w);
  free(wr);
  free(Fmnm_pad);
  free(Gmnm);
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(expsmm);

  if (verbosity > 0)
    printf("%sForward transform computed!\n", SO3_PROMPT);
}

/*!
 * Compute inverse Wigner transform for a real signal directly (without using
 * SSHT).
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in] flmn Harmonic coefficients for n >= 0. Note that for n = 0, these have to
 *                 respect the symmetry flm0* = (-1)^(m+n)*fl-m0, and hence fl00 has to
 * be real. \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_inverse_direct
 *                        \endlink instead for complex signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_inverse_direct_real(
    double *f, const complex double *flmn, const so3_parameters_t *parameters) {

  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Add optimisations for all n-modes.
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_inverse_direct with storage method %d...\n",
          SO3_PROMPT,
          storage);
  }

  // Iterators
  int el, m, n, mm; // mm for m'

  // Allocate memory.
  double *sqrt_tbl = calloc(2 * (L - 1) + 2, sizeof(*sqrt_tbl));
  SO3_ERROR_MEM_ALLOC_CHECK(sqrt_tbl);
  double *signs = calloc(L + 1, sizeof(*signs));
  SO3_ERROR_MEM_ALLOC_CHECK(signs);
  complex double *exps = calloc(4, sizeof(*exps));
  SO3_ERROR_MEM_ALLOC_CHECK(exps);

  // Perform precomputations.
  for (el = 0; el <= 2 * L - 1; ++el)
    sqrt_tbl[el] = sqrt((double)el);
  for (m = 0; m <= L - 1; m += 2) {
    signs[m] = 1.0;
    signs[m + 1] = -1.0;
  }
  int i;
  for (i = 0; i < 4; ++i)
    exps[i] = cexp(I * SO3_PION2 * i);

  // Compute Fmnm'
  // TODO: Currently m is fastest-varying, then n, then m'.
  // Should this order be changed to m-m'-n?
  complex double *Fmnm = calloc((2 * L - 1) * (2 * L - 1) * N, sizeof(*Fmnm));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm);
  int m_offset = L - 1;
  int m_stride = 2 * L - 1;
  int n_offset = 0;
  int n_stride = N;
  int mm_offset = L - 1;
  // unused: int mm_stride = 2*L-1;

  int n_start, n_stop, n_inc;

  double *dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SO3_ERROR_MEM_ALLOC_CHECK(dl);
  double *dl8 = NULL;
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SO3_ERROR_MEM_ALLOC_CHECK(dl8);
  }
  int dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  int dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);

  complex double *mn_factors = calloc((2 * L - 1) * N, sizeof *mn_factors);
  SO3_ERROR_MEM_ALLOC_CHECK(mn_factors);

  // TODO: SSHT starts this loop from MAX(L0, abs(spin)).
  // Can we use a similar optimisation? el can probably
  // be limited by n, but then we'd need to switch the
  // loop order, which means we'd have to recompute the
  // Wigner plane for each n. That seems wrong?
  for (el = L0; el <= L - 1; ++el) {
    int eltmp;
    // Compute Wigner plane.
    switch (dl_method) {
    case SSHT_DL_RISBO:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_beta_risbo_eighth_table(
              dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, eltmp, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      } else {
        ssht_dl_beta_risbo_eighth_table(
            dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, el, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      }
      break;

    case SSHT_DL_TRAPANI:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, eltmp, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      } else {
        ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, el, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      }
      break;

    default:
      SO3_ERROR_GENERIC("Invalid dl method");
    }

    // Compute Fmnm' contribution for current el.

    // Factor which depends only on el.
    double elfactor = (2.0 * el + 1.0) / (8.0 * SO3_PI * SO3_PI);

    switch (n_mode) {
    case SO3_N_MODE_ALL:
      n_start = 0;
      n_stop = MIN(N - 1, el);
      n_inc = 1;
      break;
    case SO3_N_MODE_EVEN:
      n_start = 0;
      n_stop = MIN(N - 1, el);
      n_stop -= n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_ODD:
      n_start = 1;
      n_stop = MIN(N - 1, el);
      n_stop -= 1 - n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_MAXIMUM:
      if (el < N - 1)
        continue;
      n_start = N - 1;
      n_stop = N - 1;
      n_inc = 1;
      break;
    case SO3_N_MODE_L:
      if (el >= N)
        continue;
      n_start = el;
      n_stop = el;
      n_inc = 1;
      break;
    default:
      SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    // Factors which do not depend on m'.
    for (n = n_start; n <= n_stop; n += n_inc)
      for (m = -el; m <= el; ++m) {
        int ind;
        so3_sampling_elmn2ind_real(&ind, el, m, n, parameters);
        int mod = ((n - m) % 4 + 4) % 4;
        mn_factors[m + m_offset + m_stride * (n + n_offset)] = flmn[ind] * exps[mod];
      }

    for (mm = 0; mm <= el; ++mm) {
      // These signs are needed for the symmetry relations of
      // Wigner symbols.
      double elmmsign = signs[el] * signs[mm];

      for (n = n_start; n <= n_stop; n += n_inc) {
        // Factor which does not depend on m.
        double elnmm_factor = elfactor * dl[n + dl_offset + mm * dl_stride];
        for (m = -el; m < 0; ++m)
          Fmnm
              [m + m_offset +
               m_stride * (n + n_offset + n_stride * (mm + mm_offset))] +=
              elnmm_factor * mn_factors[m + m_offset + m_stride * (n + n_offset)] *
              elmmsign * dl[-m + dl_offset + mm * dl_stride];
        for (m = 0; m <= el; ++m)
          Fmnm
              [m + m_offset +
               m_stride * (n + n_offset + n_stride * (mm + mm_offset))] +=
              elnmm_factor * mn_factors[m + m_offset + m_stride * (n + n_offset)] *
              dl[m + dl_offset + mm * dl_stride];
      }
    }
  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  switch (n_mode) {
  case SO3_N_MODE_ALL:
  case SO3_N_MODE_L:
    n_start = 0;
    n_stop = N - 1;
    n_inc = 1;
    break;
  case SO3_N_MODE_EVEN:
    n_start = 0;
    n_stop = ((N - 1) % 2 == 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_ODD:
    n_start = 1;
    n_stop = ((N - 1) % 2 != 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_MAXIMUM:
    n_start = N - 1;
    n_stop = N - 1;
    n_inc = 1;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid n-mode.");
  }

  // Use symmetry to compute Fmnm' for negative m'.
  for (mm = -L + 1; mm < 0; ++mm)
    for (n = n_start; n <= n_stop; n += n_inc)
      for (m = -L + 1; m <= L - 1; ++m)
        Fmnm[m + m_offset + m_stride * (n + n_offset + n_stride * (mm + mm_offset))] =
            signs[abs(m + n) % 2] *
            Fmnm
                [m + m_offset +
                 m_stride * (n + n_offset + n_stride * (-mm + mm_offset))];

  // Apply phase modulation to account for sampling offset.
  for (mm = -L + 1; mm <= L - 1; ++mm) {
    complex double mmfactor = cexp(I * mm * SO3_PI / (2.0 * L - 1.0));
    for (n = n_start; n <= n_stop; n += n_inc)
      for (m = -L + 1; m <= L - 1; ++m)
        Fmnm[m + m_offset + m_stride * (n + n_offset + n_stride * (mm + mm_offset))] *=
            mmfactor;
  }

  // Allocate space for shifted Fmnm'.
  complex double *Fmnm_shift =
      calloc((2 * L - 1) * (2 * L - 1) * N, sizeof(*Fmnm_shift));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm_shift);

  // Allocate space for function values.
  double *fext = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*fext));
  SO3_ERROR_MEM_ALLOC_CHECK(fext);

  // Set up plan before initialising array.
  // The redundant dimension needs to be the last one.
  fftw_plan plan = fftw_plan_dft_c2r_3d(
      2 * L - 1, 2 * L - 1, 2 * N - 1, Fmnm_shift, fext, FFTW_ESTIMATE);

  // Apply spatial shift.
  // This also reshapes the array to make n the inner dimension.
  for (mm = -L + 1; mm <= L - 1; ++mm) {
    int mm_shift = mm < 0 ? 2 * L - 1 : 0;
    for (n = n_start; n <= n_stop; n += n_inc) {
      for (m = -L + 1; m <= L - 1; ++m) {
        int m_shift = m < 0 ? 2 * L - 1 : 0;
        Fmnm_shift[n + n_stride * (m + m_shift + m_stride * (mm + mm_shift))] = Fmnm
            [m + m_offset + m_stride * (n + n_offset + n_stride * (mm + mm_offset))];
      }
    }
  }

  // Free Fmnm' memory.
  free(Fmnm);

  // Perform 3D FFT.
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  free(Fmnm_shift);

  // Extract f from the extended torus.
  // Again, we reshape the array in the process.
  int a, b, g;
  int a_stride = 2 * L - 1;
  // unused: int b_ext_stride = 2*L-1;
  int b_stride = L;
  int g_stride = 2 * N - 1;
  for (g = 0; g < 2 * N - 1; ++g)
    for (b = 0; b < L; ++b)
      for (a = 0; a < 2 * L - 1; ++a)
        f[a + a_stride * (b + b_stride * (g))] =
            fext[g + g_stride * (a + a_stride * (b))];

  // Free fext memory.
  free(fext);

  if (verbosity > 0)
    printf("%sInverse transform computed!\n", SO3_PROMPT);

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(mn_factors);
}

/*!
 * Compute forward Wigner transform for a real signal directly (without using
 * SSHT).
 *
 * \param[out] flmn Harmonic coefficients. If \link so3_parameters_t::n_mode n_mode
 *                  \endlink is different from \link SO3_N_MODE_ALL \endlink,
 *                  this array has to be nulled before being past to the function.
 * \param[in] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality \endlink flag
 *                        is ignored. Use \link so3_core_forward_direct
 *                        \endlink instead for complex signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_forward_direct_real(
    complex double *flmn, const double *f, const so3_parameters_t *parameters) {
  int L0, L, N;
  so3_sampling_t sampling;
  so3_storage_t storage;
  so3_n_mode_t n_mode;
  ssht_dl_method_t dl_method;
  int steerable;
  int verbosity;

  L0 = parameters->L0;
  L = parameters->L;
  N = parameters->N;
  sampling = parameters->sampling_scheme;
  storage = parameters->storage;
  // TODO: Add optimisations for all n-modes.
  n_mode = parameters->n_mode;
  dl_method = parameters->dl_method;
  verbosity = parameters->verbosity;
  steerable = parameters->steerable;

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%sComputing forward transform using MW sampling with\n", SO3_PROMPT);
    printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
    if (verbosity > 1)
      printf(
          "%sUsing routine so3_core_mw_forward_direct with storage method %d...\n",
          SO3_PROMPT,
          storage);
  }

  int m_stride = 2 * L - 1;
  int m_offset = L - 1;
  int n_offset = 0;
  int n_stride = N;
  int mm_stride = 2 * L - 1;
  int mm_offset = L - 1;
  int a_stride = 2 * L - 1;
  int b_stride = L;
  int bext_stride = 2 * L - 1;
  int g_stride = 2 * N - 1;

  int n_start, n_stop, n_inc;

  switch (n_mode) {
  case SO3_N_MODE_ALL:
  case SO3_N_MODE_L:
    n_start = 0;
    n_stop = N - 1;
    n_inc = 1;
    break;
  case SO3_N_MODE_EVEN:
    n_start = 0;
    n_stop = ((N - 1) % 2 == 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_ODD:
    n_start = 1;
    n_stop = ((N - 1) % 2 != 0) ? N - 1 : N - 2;
    n_inc = 2;
    break;
  case SO3_N_MODE_MAXIMUM:
    n_start = N - 1;
    n_stop = N - 1;
    n_inc = 1;
    break;
  default:
    SO3_ERROR_GENERIC("Invalid n-mode.");
  }

  double *sqrt_tbl = calloc(2 * (L - 1) + 2, sizeof(*sqrt_tbl));
  SO3_ERROR_MEM_ALLOC_CHECK(sqrt_tbl);
  double *signs = calloc(L + 1, sizeof(*signs));
  SO3_ERROR_MEM_ALLOC_CHECK(signs);
  complex double *exps = calloc(4, sizeof(*exps));
  SO3_ERROR_MEM_ALLOC_CHECK(exps);
  complex double *expsmm = calloc(2 * L - 1, sizeof(*expsmm));
  SO3_ERROR_MEM_ALLOC_CHECK(expsmm);

  int el, m, n, mm; // mm is for m'
  // Perform precomputations.
  for (el = 0; el <= 2 * (L - 1) + 1; ++el)
    sqrt_tbl[el] = sqrt((double)el);
  for (m = 0; m <= L - 1; m += 2) {
    signs[m] = 1.0;
    signs[m + 1] = -1.0;
  }
  int i;
  for (i = 0; i < 4; ++i)
    exps[i] = cexp(I * SO3_PION2 * i);
  for (mm = -L + 1; mm <= L - 1; ++mm)
    expsmm[mm + mm_offset] = cexp(-I * mm * SSHT_PI / (2.0 * L - 1.0));

  double norm_factor = 1.0 / (2.0 * L - 1.0) / (2.0 * N - 1.0);

  // Compute Fourier transform over alpha and gamma, i.e. compute Fmn(b).
  complex double *Fmnb = calloc((2 * L - 1) * (2 * L - 1) * N, sizeof(*Fmnb));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnb);
  double *fft_in = calloc((2 * L - 1) * (2 * N - 1), sizeof(*fft_in));
  SO3_ERROR_MEM_ALLOC_CHECK(fft_in);
  complex double *fft_out = calloc((2 * L - 1) * N, sizeof(*fft_out));
  SO3_ERROR_MEM_ALLOC_CHECK(fft_out);
  // Redundant dimension needs to be last
  fftw_plan plan =
      fftw_plan_dft_r2c_2d(2 * L - 1, 2 * N - 1, fft_in, fft_out, FFTW_ESTIMATE);

  int a, b, g;
  for (b = 0; b < L; ++b) {
    // TODO: This loop could probably be avoided by using
    // a more elaborate FFTW plan which performs the FFT directly
    // over the 1st and 3rd dimensions of f.
    // Instead, for each index in the 2nd dimension, we copy the
    // corresponding values in the 1st and 3rd dimension into a
    // new 2D array, to perform a standard 2D FFT there. While
    // we're at it, we also reshape that array such that gamma
    // is the inner dimension, as required by FFTW.
    for (a = 0; a < 2 * L - 1; ++a)
      for (g = 0; g < 2 * N - 1; ++g)
        fft_in[g + g_stride * (a)] = f[a + a_stride * (b + b_stride * (g))];

    fftw_execute(plan);

    // Apply spatial shift and normalisation factor, while
    // reshaping the dimensions once more.
    for (n = n_start; n <= n_stop; n += n_inc) {
      for (m = -L + 1; m <= L - 1; ++m) {
        int m_shift = m < 0 ? 2 * L - 1 : 0;
        Fmnb[b + bext_stride * (m + m_offset + m_stride * (n + n_offset))] =
            fft_out[n + n_stride * (m + m_shift)] * norm_factor;
      }
    }
  }
  fftw_destroy_plan(plan);

  // Extend Fmnb periodically.
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m) {
      int signmn = signs[abs(m + n) % 2];
      for (b = L; b < 2 * L - 1; ++b)
        Fmnb[b + bext_stride * (m + m_offset + m_stride * (n + n_offset))] =
            signmn * Fmnb
                         [(2 * L - 2 - b) +
                          bext_stride * (m + m_offset + m_stride * (n + n_offset))];
    }

  // Compute Fourier transform over beta, i.e. compute Fmnm'.
  complex double *Fmnm = calloc((2 * L - 1) * (2 * L - 1) * (2 * N - 1), sizeof(*Fmnm));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm);
  complex double *inout = calloc(2 * L - 1, sizeof(*inout));
  SO3_ERROR_MEM_ALLOC_CHECK(inout);

  plan = fftw_plan_dft_1d(2 * L - 1, inout, inout, FFTW_FORWARD, FFTW_ESTIMATE);
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m) {
      memcpy(
          inout,
          Fmnb + 0 + bext_stride * (m + m_offset + m_stride * (n + n_offset)),
          bext_stride * sizeof(*Fmnb));
      fftw_execute(plan);

      // Apply spatial shift and normalisation factor
      for (mm = -L + 1; mm <= L - 1; ++mm) {
        int mm_shift = mm < 0 ? 2 * L - 1 : 0;
        Fmnm[mm + mm_offset + mm_stride * (m + m_offset + m_stride * (n + n_offset))] =
            inout[mm + mm_shift] / (2.0 * L - 1.0);
      }
    }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m)
      for (mm = -L + 1; mm <= L - 1; ++mm)
        Fmnm[mm + mm_offset + mm_stride * (m + m_offset + m_stride * (n + n_offset))] *=
            expsmm[mm + mm_offset];

  // Compute weights.
  complex double *w = calloc(4 * L - 3, sizeof(*w));
  SO3_ERROR_MEM_ALLOC_CHECK(w);
  int w_offset = 2 * (L - 1);
  for (mm = -2 * (L - 1); mm <= 2 * (L - 1); ++mm)
    w[mm + w_offset] = so3_sampling_weight(parameters, mm);

  // Compute IFFT of w to give wr.
  complex double *wr = calloc(4 * L - 3, sizeof(*w));
  SO3_ERROR_MEM_ALLOC_CHECK(wr);
  inout = calloc(4 * L - 3, sizeof(*inout));
  SO3_ERROR_MEM_ALLOC_CHECK(inout);
  fftw_plan plan_bwd =
      fftw_plan_dft_1d(4 * L - 3, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  fftw_plan plan_fwd =
      fftw_plan_dft_1d(4 * L - 3, inout, inout, FFTW_FORWARD, FFTW_MEASURE);

  // Apply spatial shift.
  for (mm = 1; mm <= 2 * L - 2; ++mm)
    inout[mm + w_offset] = w[mm - 2 * (L - 1) - 1 + w_offset];
  for (mm = -2 * (L - 1); mm <= 0; ++mm)
    inout[mm + w_offset] = w[mm + 2 * (L - 1) + w_offset];

  fftw_execute_dft(plan_bwd, inout, inout);

  // Apply spatial shift.
  for (mm = 0; mm <= 2 * L - 2; ++mm)
    wr[mm + w_offset] = inout[mm - 2 * (L - 1) + w_offset];
  for (mm = -2 * (L - 1); mm <= -1; ++mm)
    wr[mm + w_offset] = inout[mm + 2 * (L - 1) + 1 + w_offset];

  // Compute Gmnm' by convolution implemented as product in real space.
  complex double *Fmnm_pad = calloc(4 * L - 3, sizeof(*Fmnm_pad));
  SO3_ERROR_MEM_ALLOC_CHECK(Fmnm_pad);
  complex double *Gmnm = calloc((2 * L - 1) * (2 * L - 1) * N, sizeof(*Gmnm));
  SO3_ERROR_MEM_ALLOC_CHECK(Gmnm);
  for (n = n_start; n <= n_stop; n += n_inc)
    for (m = -L + 1; m <= L - 1; ++m) {

      // Zero-pad Fmnm'.
      for (mm = -2 * (L - 1); mm <= -L; ++mm)
        Fmnm_pad[mm + w_offset] = 0.0;
      for (mm = L; mm <= 2 * (L - 1); ++mm)
        Fmnm_pad[mm + w_offset] = 0.0;
      for (mm = -(L - 1); mm <= L - 1; ++mm)
        Fmnm_pad[mm + w_offset] = Fmnm
            [mm + mm_offset + mm_stride * (m + m_offset + m_stride * (n + n_offset))];

      // Apply spatial shift.
      for (mm = 1; mm <= 2 * L - 2; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm - 2 * (L - 1) - 1 + w_offset];
      for (mm = -2 * (L - 1); mm <= 0; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm + 2 * (L - 1) + w_offset];
      // Compute IFFT of Fmnm'.
      fftw_execute_dft(plan_bwd, inout, inout);
      // Apply spatial shift.
      for (mm = 0; mm <= 2 * L - 2; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm - 2 * (L - 1) + w_offset];
      for (mm = -2 * (L - 1); mm <= -1; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm + 2 * (L - 1) + 1 + w_offset];

      // Compute product of Fmnm' and weight in real space.
      int r;
      for (r = -2 * (L - 1); r <= 2 * (L - 1); ++r)
        Fmnm_pad[r + w_offset] *= wr[r + w_offset];

      // Apply spatial shift.
      for (mm = 1; mm <= 2 * L - 2; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm - 2 * (L - 1) - 1 + w_offset];
      for (mm = -2 * (L - 1); mm <= 0; ++mm)
        inout[mm + w_offset] = Fmnm_pad[mm + 2 * (L - 1) + w_offset];
      // Compute Gmnm' by FFT.
      fftw_execute_dft(plan_fwd, inout, inout);
      // Apply spatial shift.
      for (mm = 0; mm <= 2 * L - 2; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm - 2 * (L - 1) + w_offset];
      for (mm = -2 * (L - 1); mm <= -1; ++mm)
        Fmnm_pad[mm + w_offset] = inout[mm + 2 * (L - 1) + 1 + w_offset];

      // Extract section of Gmnm' of interest.
      for (mm = -(L - 1); mm <= L - 1; ++mm)
        Gmnm[m + m_offset + m_stride * (mm + mm_offset + mm_stride * (n + n_offset))] =
            Fmnm_pad[mm + w_offset] * 4.0 * SSHT_PI * SSHT_PI / (4.0 * L - 3.0);
    }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flmn.
  double *dl, *dl8 = NULL;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SO3_ERROR_MEM_ALLOC_CHECK(dl);
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SO3_ERROR_MEM_ALLOC_CHECK(dl8);
  }
  int dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  int dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  for (n = 0; n <= N - 1; ++n)
    for (el = n; el < L; ++el)
      for (m = -el; m <= el; ++m) {
        int ind;
        so3_sampling_elmn2ind_real(&ind, el, m, n, parameters);
        flmn[ind] = 0.0;
      }

  for (el = L0; el < L; ++el) {
    int eltmp;

    // Compute Wigner plane.
    switch (dl_method) {
    case SSHT_DL_RISBO:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_beta_risbo_eighth_table(
              dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, eltmp, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      } else {
        ssht_dl_beta_risbo_eighth_table(
            dl8, SO3_PION2, L, SSHT_DL_QUARTER_EXTENDED, el, sqrt_tbl, signs);
        ssht_dl_beta_risbo_fill_eighth2quarter_table(
            dl, dl8, L, SSHT_DL_QUARTER, SSHT_DL_QUARTER_EXTENDED, el, signs);
      }
      break;

    case SSHT_DL_TRAPANI:
      if (el != 0 && el == L0) {
        for (eltmp = 0; eltmp <= L0; ++eltmp)
          ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, eltmp, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      } else {
        ssht_dl_halfpi_trapani_eighth_table(dl, L, SSHT_DL_QUARTER, el, sqrt_tbl);
        ssht_dl_halfpi_trapani_fill_eighth2quarter_table(
            dl, L, SSHT_DL_QUARTER, el, signs);
      }
      break;

    default:
      SO3_ERROR_GENERIC("Invalid dl method");
    }

    // Compute flmn for current el.

    switch (n_mode) {
    case SO3_N_MODE_ALL:
      n_start = 0;
      n_stop = MIN(N - 1, el);
      n_inc = 1;
      break;
    case SO3_N_MODE_EVEN:
      n_start = 0;
      n_stop = MIN(N - 1, el);
      n_stop -= n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_ODD:
      n_start = 1;
      n_stop = MIN(N - 1, el);
      n_stop -= 1 - n_stop % 2;
      n_inc = 2;
      break;
    case SO3_N_MODE_MAXIMUM:
      if (el < N - 1)
        continue;
      n_start = N - 1;
      n_stop = N - 1;
      n_inc = 1;
      break;
    case SO3_N_MODE_L:
      if (el >= N)
        continue;
      n_start = el;
      n_stop = el;
      n_inc = 1;
      break;
    default:
      SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    // TODO: Pull out a few multiplications into precomputations
    // or split up loops to avoid conditionals to check signs.
    for (mm = -el; mm <= el; ++mm) {
      // These signs are needed for the symmetry relations of
      // Wigner symbols.
      double elmmsign = signs[el] * signs[abs(mm)];

      for (n = n_start; n <= n_stop; n += n_inc) {
        double mmsign = mm >= 0 ? 1.0 : signs[el] * signs[n];

        // Factor which does not depend on m.
        double elnmm_factor = mmsign * dl[n + dl_offset + abs(mm) * dl_stride];

        for (m = -el; m <= el; ++m) {
          mmsign = mm >= 0 ? 1.0 : signs[el] * signs[abs(m)];
          double elmsign = m >= 0 ? 1.0 : elmmsign;
          int ind;
          so3_sampling_elmn2ind_real(&ind, el, m, n, parameters);
          int mod = ((m - n) % 4 + 4) % 4;
          flmn[ind] += exps[mod] * elnmm_factor * mmsign * elmsign *
                       dl[abs(m) + dl_offset + abs(mm) * dl_stride] *
                       Gmnm
                           [m + m_offset +
                            m_stride * (mm + mm_offset + mm_stride * (n + n_offset))];
        }
      }
    }
  }

  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmnb);
  free(Fmnm);
  free(inout);
  free(w);
  free(wr);
  free(Fmnm_pad);
  free(Gmnm);
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(expsmm);

  if (verbosity > 0)
    printf("%sForward transform computed!\n", SO3_PROMPT);
}

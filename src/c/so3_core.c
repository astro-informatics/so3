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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>

#include "ssht.h"

#include "so3_types.h"
#include "so3_error.h"
#include "so3_sampling.h"

#define MAX(a,b) ((a > b) ? (a) : (b))

typedef void (*inverse_complex_ssht)(complex double *, const complex double *, int, int, int, ssht_dl_method_t, int);
typedef void (*inverse_real_ssht)(double *, const complex double *, int, int, ssht_dl_method_t, int);
typedef void (*forward_complex_ssht)(complex double *, const complex double *, int, int, int, ssht_dl_method_t, int);
typedef void (*forward_real_ssht)(complex double *, const double *, int, int, ssht_dl_method_t, int);

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
    complex double *f, const complex double *flmn,
    const so3_parameters_t *parameters
) {

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
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    verbosity = parameters->verbosity;
    steerable = parameters->steerable;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_inverse_via_ssht with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }

    switch (sampling)
    {
    case SO3_SAMPLING_MW:
        fn_n_stride = L * (2*L-1);
        ssht = ssht_core_mw_lb_inverse_sov_sym;
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        ssht = ssht_core_mw_lb_inverse_sov_sym_ss;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

    // Compute fn(a,b)

    if (steerable)
    {
        // For steerable signals, we need to supersample in n/gamma,
        // in order to create a symmetric sampling.
        fftw_n = 2*N; // Each transform is over 2*N

        // We need to perform the FFT into a temporary buffer, because
        // the result will be twice as large as the output we need.
        ftemp = malloc(2*N*fn_n_stride * sizeof *ftemp);
        SO3_ERROR_MEM_ALLOC_CHECK(ftemp);

        fftw_target = ftemp;
    }
    else
    {
        fftw_n = 2*N-1; // Each transform is over 2*N-1

        fftw_target = f;
    }

    fn = calloc(fftw_n*fn_n_stride, sizeof *fn);
    SO3_ERROR_MEM_ALLOC_CHECK(fn);



    // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
    // necessary but still good practice.
    fftw_rank = 1; // We compute 1d transforms
    fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

    // We want to transform columns
    fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
    fftw_istride = fftw_ostride = fn_n_stride; // Distance between two elements of the same column

    plan = fftw_plan_many_dft(
            fftw_rank, &fftw_n, fftw_howmany,
            fn, NULL, fftw_istride, fftw_idist,
            fftw_target, NULL, fftw_ostride, fftw_odist,
            FFTW_BACKWARD, FFTW_ESTIMATE
    );

    for(n = -N+1; n <= N-1; ++n)
    {
        int ind, offset, i, el;
        int L0e = MAX(L0, abs(n)); // 'e' for 'effective'
        double factor;
        complex double *flm;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

        flm = malloc(L*L * sizeof *flm);
        SO3_ERROR_MEM_ALLOC_CHECK(flm);

        switch (storage)
        {
        case SO3_STORAGE_PADDED:
            so3_sampling_elmn2ind(&ind, 0, 0, n, parameters);
            memcpy(flm, flmn + ind, L*L * sizeof(complex double));
            break;
        case SO3_STORAGE_COMPACT:
            so3_sampling_elmn2ind(&ind, abs(n), -abs(n), n, parameters);
            memcpy(flm + n*n, flmn + ind, (L*L - n*n) * sizeof(complex double));
            for(i = 0; i < n*n; ++i)
                flm[i] = 0.0;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid storage method.");
        }

        el = L0e;
        i = offset = el*el;
        for(; el < L; ++el)
        {
            factor = sqrt((double)(2*el+1)/(16.*pow(SO3_PI, 3.)));
            for (; i < offset + 2*el+1; ++i)
                flm[i] *= factor;

            offset = i;
        }

        // The conditional applies the spatial transform, so that we store
        // the results in n-order 0, 1, 2, -2, -1
        offset = (n < 0 ? n + fftw_n : n);

        (*ssht)(
            fn + offset*fn_n_stride, flm,
            L0e, L, -n,
            dl_method,
            verbosity
        );

        if(n % 2)
            for(i = 0; i < fn_n_stride; ++i)
                fn[offset*fn_n_stride + i] = -fn[offset*fn_n_stride + i];

        if (verbosity > 0)
            printf("\n");

        free(flm);
    }


    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if (steerable)
    {
        memcpy(f, ftemp, N*fn_n_stride * sizeof(complex double));
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
    complex double *flmn, const complex double *f,
    const so3_parameters_t *parameters
) {
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
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    verbosity = parameters->verbosity;
    steerable = parameters->steerable;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing forward transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_forward_via_ssht with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }

    switch (sampling)
    {
    case SO3_SAMPLING_MW:
        fn_n_stride = L * (2*L-1);
        ssht = ssht_core_mw_lb_forward_sov_conv_sym;
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        ssht = ssht_core_mw_lb_forward_sov_conv_sym_ss;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

    if (steerable)
    {
        int g, offset;

        fn = calloc((2*N-1)*fn_n_stride, sizeof *fn);
        SO3_ERROR_MEM_ALLOC_CHECK(fn);

        for (n = -N+1; n < N; n+=2)
        {
            // The conditional applies the spatial transform, because the fn
            // are to be stored in n-order 0, 1, 2, -2, -1
            offset = (n < 0 ? n + 2*N-1 : n);

            for (g = 0; g < N; ++g)
            {
                double gamma = g * SO3_PI / N;
                for (i = 0; i < fn_n_stride; ++i)
                {
                    double weight = 2*SO3_PI/N;
                    fn[offset * fn_n_stride + i] += weight*f[g * fn_n_stride + i]*cexp(-I*n*gamma);
                }
            }
        }
    }
    else
    {
        // Make a copy of the input, because input is const
        // This could potentially be avoided by copying the input into fn and using an
        // in-place FFTW. The performance impact has to be profiled, though.
        ftemp = malloc((2*N-1)*fn_n_stride * sizeof *ftemp);
        SO3_ERROR_MEM_ALLOC_CHECK(ftemp);
        memcpy(ftemp, f, (2*N-1)*fn_n_stride * sizeof(complex double));

        fn = malloc((2*N-1)*fn_n_stride * sizeof *fn);
        SO3_ERROR_MEM_ALLOC_CHECK(fn);

        // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
        // necessary but still good practice.
        fftw_rank = 1;
        fftw_n = 2*N-1;
        fftw_howmany = fn_n_stride;
        fftw_idist = fftw_odist = 1;
        fftw_istride = fftw_ostride = fn_n_stride;

        plan = fftw_plan_many_dft(
                fftw_rank, &fftw_n, fftw_howmany,
                ftemp, NULL, fftw_istride, fftw_idist,
                fn, NULL, fftw_ostride, fftw_odist,
                FFTW_FORWARD, FFTW_ESTIMATE
        );

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        free(ftemp);

        factor = 2*SO3_PI/(double)(2*N-1);
        for(i = 0; i < (2*N-1)*fn_n_stride; ++i)
            fn[i] *= factor;
    }

    for(n = -N+1; n <= N-1; ++n)
    {
        int ind, offset, el, sign;
        int L0e = MAX(L0, abs(n)); // 'e' for 'effective'

        complex double *flm = NULL;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

        if (storage == SO3_STORAGE_COMPACT)
            flm = malloc(L*L * sizeof *flm);

        // The conditional applies the spatial transform, because the fn
        // are stored in n-order 0, 1, 2, -2, -1
        offset = (n < 0 ? n + 2*N-1 : n);

        complex double *flm_block;
        complex double *fn_block = fn + offset*fn_n_stride;

        el = L0e;
        switch (storage)
        {
        case SO3_STORAGE_PADDED:
            so3_sampling_elmn2ind(&ind, 0, 0, n, parameters);
            flm_block = flmn + ind;
            i = offset = el*el;
            break;
        case SO3_STORAGE_COMPACT:
            flm_block = flm;
            i = offset = el*el-n*n;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid storage method.");
        }

        (*ssht)(
            flm_block, fn_block,
            L0e, L, -n,
            dl_method,
            verbosity
        );

        if (storage == SO3_STORAGE_COMPACT)
        {
            so3_sampling_elmn2ind(&ind, abs(n), -abs(n), n, parameters);
            memcpy(flmn + ind, flm + n*n, (L*L - n*n) * sizeof(complex double));
        }

        if (n % 2)
            sign = -1;
        else
            sign = 1;

        for(; el < L; ++el)
        {
            factor = sign*sqrt(4.0*SO3_PI/(double)(2*el+1));
            for (; i < offset + 2*el+1; ++i)
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
 *                 respect the symmetry flm0* = (-1)^(m+n)*fl-m0, and hence fl00 has to be real.
 * \param[in]  parameters A fully populated parameters object. The \link
 *                        so3_parameters_t::reality reality\endlink flag
 *                        is ignored. Use \link so3_core_inverse_via_ssht
 *                        \endlink instead for complex signals.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_inverse_via_ssht_real(
    double *f, const complex double *flmn,
    const so3_parameters_t *parameters
) {
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
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    verbosity = parameters->verbosity;
    steerable = parameters->steerable;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_inverse_via_ssht_real with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }

    switch (sampling)
    {
    case SO3_SAMPLING_MW:
        fn_n_stride = L * (2*L-1);
        complex_ssht = ssht_core_mw_lb_inverse_sov_sym;
        real_ssht = ssht_core_mw_lb_inverse_sov_sym_real;
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        complex_ssht = ssht_core_mw_lb_inverse_sov_sym_ss;
        real_ssht = ssht_core_mw_lb_inverse_sov_sym_ss_real;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

    // Compute fn(a,b)

    // For steerable signals, we need to supersample in n/gamma,
    // in order to create a symmetric sampling.
    // Each transform is over fftw_n samples (logically; physically, fn for negative n will be omitted)
    if (steerable)
    {
        // For steerable signals, we need to supersample in n/gamma,
        // in order to create a symmetric sampling.
        fftw_n = 2*N;

        // We need to perform the FFT into a temporary buffer, because
        // the result will be twice as large as the output we need.
        ftemp = malloc(2*N*fn_n_stride * sizeof *ftemp);
        SO3_ERROR_MEM_ALLOC_CHECK(ftemp);

        fftw_target = ftemp;
    }
    else
    {
        fftw_n = 2*N-1;

        fftw_target = f;
    }

    // Only need to store for non-negative n
    fn = calloc((fftw_n/2+1)*fn_n_stride, sizeof *fn);
    SO3_ERROR_MEM_ALLOC_CHECK(fn);

    // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
    // necessary but still good practice.
    fftw_rank = 1; // We compute 1d transforms
    fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

    // We want to transform columns
    fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
    fftw_istride = fftw_ostride = fn_n_stride; // Distance between two elements of the same column

    plan = fftw_plan_many_dft_c2r(
            fftw_rank, &fftw_n, fftw_howmany,
            fn, NULL, fftw_istride, fftw_idist,
            fftw_target, NULL, fftw_ostride, fftw_odist,
            FFTW_ESTIMATE
    );

    flm = malloc(L*L * sizeof *flm);
    SO3_ERROR_MEM_ALLOC_CHECK(flm);

    for(n = 0; n <= N-1; ++n)
    {
        int ind, offset, i, el;
        int L0e = MAX(L0, abs(n)); // 'e' for 'effective'
        double factor;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

        switch (storage)
        {
        case SO3_STORAGE_PADDED:
            so3_sampling_elmn2ind_real(&ind, 0, 0, n, parameters); //L, N, SO3_STORE_NEG_FIRST_PAD);
            memcpy(flm, flmn + ind, L*L * sizeof(complex double));
            break;
        case SO3_STORAGE_COMPACT:
            so3_sampling_elmn2ind_real(&ind, n, -n, n, parameters); //L, N, SO3_STORE_NEG_FIRST_COMPACT);
            memcpy(flm + n*n, flmn + ind, (L*L - n*n) * sizeof(complex double));
            for(i = 0; i < n*n; ++i)
                flm[i] = 0.0;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid storage method.");
        }

        el = L0e;
        i = offset = el*el;
        for(; el < L; ++el)
        {
            factor = sqrt((double)(2*el+1)/(16.*pow(SO3_PI, 3.)));
            for (; i < offset + 2*el+1; ++i)
                flm[i] *= factor;

            offset = i;
        }

        if (N > 1 || n)
        {
            (*complex_ssht)(
                fn + n*fn_n_stride, flm,
                L0e, L, -n,
                dl_method,
                verbosity
            );
        }
        else
        {
            double *fn_r;

            // Create an array of real doubles for n = 0
            fn_r = calloc(fn_n_stride, sizeof *fn_r);
            SO3_ERROR_MEM_ALLOC_CHECK(fn_r);

            (*real_ssht)(
                fn_r, flm,
                L0e, L,
                dl_method,
                verbosity
            );

            for (i = 0; i < fn_n_stride; ++i)
                fn[i] = fn_r[i];
        }


        if(n % 2)
            for(i = 0; i < fn_n_stride; ++i)
                fn[n*fn_n_stride + i] = -fn[n*fn_n_stride + i];

        if (verbosity > 0)
            printf("\n");
    }

    free(flm);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if (steerable)
    {
        memcpy(f, ftemp, N*fn_n_stride * sizeof *f);
        free(ftemp);
    }

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
    complex double *flmn, const double *f,
    const so3_parameters_t *parameters
) {
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
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    steerable = parameters->steerable;
    verbosity = parameters->verbosity;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing forward transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_forward_via_ssht_real with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }

    switch (sampling)
    {
    case SO3_SAMPLING_MW:
        fn_n_stride = L * (2*L-1);
        complex_ssht = ssht_core_mw_lb_forward_sov_conv_sym;
        real_ssht = ssht_core_mw_lb_forward_sov_conv_sym_real;
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        complex_ssht = ssht_core_mw_lb_forward_sov_conv_sym_ss;
        real_ssht = ssht_core_mw_lb_forward_sov_conv_sym_ss_real;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }


    if (steerable)
    {
        int g, offset;

        fn = calloc((2*N-1)*fn_n_stride, sizeof *fn);
        SO3_ERROR_MEM_ALLOC_CHECK(fn);

        for (n = -N+1; n < N; n+=2)
        {
            // The conditional applies the spatial transform, because the fn
            // are to be stored in n-order 0, 1, 2, -2, -1
            offset = (n < 0 ? n + 2*N-1 : n);

            for (g = 0; g < N; ++g)
            {
                double gamma = g * SO3_PI / N;
                for (i = 0; i < fn_n_stride; ++i)
                {
                    double weight = 2*SO3_PI/N;
                    fn[offset * fn_n_stride + i] += weight*f[g * fn_n_stride + i]*cexp(-I*n*gamma);
                }
            }
        }
    }
    else
    {
        // Make a copy of the input, because input is const
        // This could potentially be avoided by copying the input into fn and using an
        // in-place FFTW. The performance impact has to be profiled, though.
        ftemp = malloc((2*N-1)*fn_n_stride * sizeof *ftemp);
        SO3_ERROR_MEM_ALLOC_CHECK(ftemp);
        memcpy(ftemp, f, (2*N-1)*fn_n_stride * sizeof(double));

        fn = malloc(N*fn_n_stride * sizeof *fn);
        SO3_ERROR_MEM_ALLOC_CHECK(fn);
        // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
        // necessary but still good practice.
        fftw_rank = 1; // We compute 1d transforms
        fftw_n = 2*N-1; // Each transform is over 2*N-1 (logically; physically, fn for negative n will be omitted)
        fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

        // We want to transform columns
        fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
        fftw_istride = fftw_ostride = fn_n_stride; // Distance between two elements of the same column

        plan = fftw_plan_many_dft_r2c(
                fftw_rank, &fftw_n, fftw_howmany,
                ftemp, NULL, fftw_istride, fftw_idist,
                fn, NULL, fftw_ostride, fftw_odist,
                FFTW_ESTIMATE
        );

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        free(ftemp);

        factor = 2*SO3_PI/(double)(2*N-1);
        for(i = 0; i < N*fn_n_stride; ++i)
            fn[i] *= factor;
    }

    if (storage == SO3_STORAGE_COMPACT)
        flm = malloc(L*L * sizeof *flm);

    for(n = 0; n <= N-1; ++n)
    {
        int ind, offset, el, sign;
        int L0e = MAX(L0, abs(n)); // 'e' for 'effective'

        complex double* flm_block;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

        el = L0e;
        i = offset = el*el;
        switch (storage)
        {
        case SO3_STORAGE_PADDED:
            so3_sampling_elmn2ind_real(&ind, 0, 0, n, parameters);
            flm_block = flmn + ind;
            break;
        case SO3_STORAGE_COMPACT:
            flm_block = flm;

            offset -= n*n;
            i = offset;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid storage method.");
        }


        if (N > 1)
        {
            (*complex_ssht)(
                flm_block, fn + n*fn_n_stride,
                L0e, L, -n,
                dl_method,
                verbosity
            );
        }
        else
        {
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
            (*real_ssht)(
                flm_block, fn_r,
                L0e, L,
                dl_method,
                verbosity
            );

            free(fn_r);
        }

        if (storage == SO3_STORAGE_COMPACT)
        {
            so3_sampling_elmn2ind_real(&ind, n, -n, n, parameters);
            memcpy(flmn + ind, flm + n*n, (L*L - n*n) * sizeof(complex double));
        }

        if (n % 2)
            sign = -1;
        else
            sign = 1;

        for(; el < L; ++el)
        {
            factor = sign*sqrt(4.0*SO3_PI/(double)(2*el+1));
            for (; i < offset + 2*el+1; ++i)
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

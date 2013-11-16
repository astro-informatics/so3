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

/*!
 * Compute inverse transform for MW method via SSHT.
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  flmn Harmonic coefficients.
 * \param[in]  parameters A fully populated parameters object. The reality flag
 *                        is ignored. Use so3_core_inverse_via_ssht_real instead
 *                        for real signals.
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
    int verbosity;

    // Iterator
    int n;
    // Intermediate results
    complex double *fn, *flm;
    // Stride for several arrays
    int fn_n_stride;
    // FFTW-related variables
    int fftw_rank, fftw_howmany;
    int fftw_idist, fftw_odist;
    int fftw_istride, fftw_ostride;
    int fftw_n;
    fftw_plan plan;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;
    sampling = parameters->sampling_scheme;
    storage = parameters->storage;
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    verbosity = parameters->verbosity;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_inverse_via_ssht with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }

    // Compute fn(a,b)

    switch (sampling)
    {
    case SO3_SAMPLING_MW:
        fn_n_stride = L * (2*L-1);
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

    fn = calloc((2*N-1)*fn_n_stride, sizeof *fn);
    SO3_ERROR_MEM_ALLOC_CHECK(fn);

    // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
    // necessary but still good practice.
    fftw_rank = 1; // We compute 1d transforms
    fftw_n = 2*N-1; // Each transform is over 2*N-1
    fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

    // We want to transform columns
    fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
    fftw_istride = fftw_ostride = fn_n_stride; // Distance between two elements of the same column

    plan = fftw_plan_many_dft(
            fftw_rank, &fftw_n, fftw_howmany,
            fn, NULL, fftw_istride, fftw_idist,
            f, NULL, fftw_ostride, fftw_odist,
            FFTW_BACKWARD, FFTW_ESTIMATE
    );

    flm = malloc(L*L * sizeof *flm);
    SO3_ERROR_MEM_ALLOC_CHECK(flm);

    for(n = -N+1; n <= N-1; ++n)
    {
        int ind, offset, i, el;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

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

        el = MAX(L0,abs(n));
        i = offset = el*el;
        for(; el < L; ++el)
        {
            for (; i < offset + 2*el+1; ++i)
                flm[i] *= sqrt((double)(2*el+1)/(16.*pow(SO3_PI, 3.)));

            offset = i;
        }

        // The conditional applies the spatial transform, so that we store
        // the results in n-order 0, 1, 2, -2, -1
        offset = (n < 0 ? n + 2*N-1 : n);


        switch (sampling)
        {
        case SO3_SAMPLING_MW:
            ssht_core_mw_lb_inverse_sov_sym(
                fn + offset*fn_n_stride, flm,
                L0, L, -n,
                dl_method,
                verbosity
            );
            break;
        case SO3_SAMPLING_MW_SS:
            ssht_core_mw_lb_inverse_sov_sym_ss(
                fn + offset*fn_n_stride, flm,
                L0, L, -n,
                dl_method,
                verbosity
            );
            break;
        default:
            SO3_ERROR_GENERIC("Invalid sampling scheme.");
        }

        if(n % 2)
            for(i = 0; i < fn_n_stride; ++i)
                fn[offset*fn_n_stride + i] = -fn[offset*fn_n_stride + i];

        if (verbosity > 0)
            printf("\n");
    }

    free(flm);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    free(fn);

    if (verbosity > 0)
        printf("%sInverse transform computed!\n", SO3_PROMPT);
}


/*!
 * Compute forward transform for MW method via SSHT.
 *
 * \param[out] flmn Harmonic coefficients. If n_mode is different from SO3_N_MODE_ALL,
 *                  this array has to be nulled before being past to the function.
 * \param[in] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  parameters A fully populated parameters object. The reality flag
 *                        is ignored. Use so3_core_forward_via_ssht_real instead
 *                        for real signals.
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
    int verbosity;

    // Iterator
    int i, n;
    // Intermediate results
    complex double *ftemp, *flm, *fn;
    // Stride for several arrays
    int fn_n_stride;
    // FFTW-related variables
    int fftw_rank, fftw_howmany;
    int fftw_idist, fftw_odist;
    int fftw_istride, fftw_ostride;
    int fftw_n;
    fftw_plan plan;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;
    sampling = parameters->sampling_scheme;
    storage = parameters->storage;
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    verbosity = parameters->verbosity;

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
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

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

    for(i = 0; i < (2*N-1)*fn_n_stride; ++i)
        fn[i] *= 2*SO3_PI/(double)(2*N-1);

    if (storage == SO3_STORAGE_COMPACT)
        flm = malloc(L*L * sizeof *flm);

    for(n = -N+1; n <= N-1; ++n)
    {
        int ind, offset, el, sign;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

        // The conditional applies the spatial transform, because the fn
        // are stored in n-order 0, 1, 2, -2, -1
        offset = (n < 0 ? n + 2*N-1 : n);

        switch (storage)
        {
        case SO3_STORAGE_PADDED:
            so3_sampling_elmn2ind(&ind, 0, 0, n, parameters);
            switch (sampling)
            {
            case SO3_SAMPLING_MW:
                ssht_core_mw_lb_forward_sov_conv_sym(
                    flmn + ind, fn + offset*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            case SO3_SAMPLING_MW_SS:
                ssht_core_mw_lb_forward_sov_conv_sym_ss(
                    flmn + ind, fn + offset*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            default:
                SO3_ERROR_GENERIC("Invalid sampling scheme.");
            }

            i = offset = L0*L0;
            el = L0;
            break;
        case SO3_STORAGE_COMPACT:
            switch (sampling)
            {
            case SO3_SAMPLING_MW:
                ssht_core_mw_lb_forward_sov_conv_sym(
                    flm, fn + offset*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            case SO3_SAMPLING_MW_SS:
                ssht_core_mw_lb_forward_sov_conv_sym_ss(
                    flm, fn + offset*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            default:
                SO3_ERROR_GENERIC("Invalid sampling scheme.");
            }

            so3_sampling_elmn2ind(&ind, abs(n), -abs(n), n, parameters);
            memcpy(flmn + ind, flm + n*n, (L*L - n*n) * sizeof(complex double));

            el = MAX(L0, abs(n));
            i = offset = el*el-n*n;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid storage method.");
        }

        if (n % 2)
            sign = -1;
        else
            sign = 1;

        for(; el < L; ++el)
        {
            for (; i < offset + 2*el+1; ++i)
                flmn[ind + i] *= sign*sqrt(4.0*SO3_PI/(double)(2*el+1));

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
 * Compute inverse transform for MW method via SSHT for a real signal.
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in] flmn Harmonic coefficients for n >= 0. Note that for n = 0, these have to
 *                 respect the symmetry flm0* = (-1)^(m+n)*fl-m0, and hence fl00 has to be real.
 * \param[in]  parameters A fully populated parameters object. The reality flag
 *                        is ignored. Use so3_core_inverse_via_ssht instead
 *                        for complex signals.
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
    int verbosity;

    // Iterator
    int n;
    // Intermediate results
    complex double *fn, *flm;
    // Stride for several arrays
    int fn_n_stride;
    // FFTW-related variables
    int fftw_rank, fftw_howmany;
    int fftw_idist, fftw_odist;
    int fftw_istride, fftw_ostride;
    int fftw_n;
    fftw_plan plan;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;
    sampling = parameters->sampling_scheme;
    storage = parameters->storage;
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
    verbosity = parameters->verbosity;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_inverse_via_ssht_real with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }

    // Compute fn(a,b)

    switch (sampling)
    {
    case SO3_SAMPLING_MW:
        fn_n_stride = L * (2*L-1);
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

    // Only need to store for non-negative n
    fn = calloc(N*fn_n_stride, sizeof *fn);
    SO3_ERROR_MEM_ALLOC_CHECK(fn);

    // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
    // necessary but still good practice.
    fftw_rank = 1; // We compute 1d transforms
    fftw_n = 2*N-1; // Each transform is over 2*N-1 (logically; physically, fn for negative n will be omitted)
    fftw_howmany = fn_n_stride; // We need L*(2*L-1) of these transforms

    // We want to transform columns
    fftw_idist = fftw_odist = 1; // The starts of the columns are contiguous in memory
    fftw_istride = fftw_ostride = fn_n_stride; // Distance between two elements of the same column

    plan = fftw_plan_many_dft_c2r(
            fftw_rank, &fftw_n, fftw_howmany,
            fn, NULL, fftw_istride, fftw_idist,
            f, NULL, fftw_ostride, fftw_odist,
            FFTW_ESTIMATE
    );

    flm = malloc(L*L * sizeof *flm);
    SO3_ERROR_MEM_ALLOC_CHECK(flm);

    for(n = 0; n <= N-1; ++n)
    {
        int ind, offset, i, el;

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

        el = MAX(L0,abs(n));
        i = offset = el*el;
        for(; el < L; ++el)
        {
            for (; i < offset + 2*el+1; ++i)
                flm[i] *= sqrt((double)(2*el+1)/(16.*pow(SO3_PI, 3.)));

            offset = i;
        }

        switch (sampling)
        {
        case SO3_SAMPLING_MW:
            ssht_core_mw_lb_inverse_sov_sym(
                fn + n*fn_n_stride, flm,
                L0, L, -n,
                dl_method,
                verbosity
            );
            break;
        case SO3_SAMPLING_MW_SS:
            ssht_core_mw_lb_inverse_sov_sym_ss(
                fn + n*fn_n_stride, flm,
                L0, L, -n,
                dl_method,
                verbosity
            );
            break;
        default:
            SO3_ERROR_GENERIC("Invalid sampling scheme.");
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

    free(fn);

    if (verbosity > 0)
        printf("%sInverse transform computed!\n", SO3_PROMPT);
}

/*!
 * Compute forward transform for MW method via SSHT for a real signal.
 *
 * \param[out] flmn Harmonic coefficients. If n_mode is different from SO3_N_MODE_ALL,
 *                  this array has to be nulled before being past to the function.
 * \param[in] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in]  parameters A fully populated parameters object. The reality flag
 *                        is ignored. Use so3_core_forward_via_ssht instead
 *                        for complex signals.
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
    int verbosity;

    // Iterator
    int i, n;
    // Intermediate results
    double *ftemp;
    complex double *flm, *fn;
    // Stride for several arrays
    int fn_n_stride;
    // FFTW-related variables
    int fftw_rank, fftw_howmany;
    int fftw_idist, fftw_odist;
    int fftw_istride, fftw_ostride;
    int fftw_n;
    fftw_plan plan;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;
    sampling = parameters->sampling_scheme;
    storage = parameters->storage;
    n_mode = parameters->n_mode;
    dl_method = parameters->dl_method;
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
        break;
    case SO3_SAMPLING_MW_SS:
        fn_n_stride = (L+1) * 2*L;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }

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

    for(i = 0; i < N*fn_n_stride; ++i)
        fn[i] *= 2*SO3_PI/(double)(2*N-1);

    if (storage == SO3_STORAGE_COMPACT)
        flm = malloc(L*L * sizeof *flm);

    for(n = 0; n <= N-1; ++n)
    {
        int ind, offset, el, sign;

        if ((n_mode == SO3_N_MODE_EVEN && n % 2)
            || (n_mode == SO3_N_MODE_ODD && !(n % 2))
            || (n_mode == SO3_N_MODE_MAXIMUM && abs(n) < N-1)
        ) {
            continue;
        }

        switch (storage)
        {
        case SO3_STORAGE_PADDED:
            so3_sampling_elmn2ind_real(&ind, 0, 0, n, parameters);
            switch (sampling)
            {
            case SO3_SAMPLING_MW:
                ssht_core_mw_lb_forward_sov_conv_sym(
                    flmn + ind, fn + n*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            case SO3_SAMPLING_MW_SS:
                ssht_core_mw_lb_forward_sov_conv_sym_ss(
                    flmn + ind, fn + n*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            default:
                SO3_ERROR_GENERIC("Invalid sampling scheme.");
            }

            i = offset = L0*L0;
            el = L0;
            break;
        case SO3_STORAGE_COMPACT:
            switch (sampling)
            {
            case SO3_SAMPLING_MW:
                ssht_core_mw_lb_forward_sov_conv_sym(
                    flm, fn + n*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            case SO3_SAMPLING_MW_SS:
                ssht_core_mw_lb_forward_sov_conv_sym_ss(
                    flm, fn + n*fn_n_stride,
                    L0, L, -n,
                    dl_method,
                    verbosity
                );
                break;
            default:
                SO3_ERROR_GENERIC("Invalid sampling scheme.");
            }
            so3_sampling_elmn2ind_real(&ind, n, -n, n, parameters);
            memcpy(flmn + ind, flm + n*n, (L*L - n*n) * sizeof(complex double));

            el = MAX(L0, abs(n));
            i = offset = el*el-n*n;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid storage method.");
        }

        if (n % 2)
            sign = -1;
        else
            sign = 1;

        for(; el < L; ++el)
        {
            for (; i < offset + 2*el+1; ++i)
                flmn[ind + i] *= sign*sqrt(4.0*SO3_PI/(double)(2*el+1));

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

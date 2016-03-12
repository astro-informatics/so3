// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_adjoint.c
 * Algorithms to perform adjoint Wigner transform on the rotation group SO(§).
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

#define MIN(a,b) ((a < b) ? (a) : (b))
#define MAX(a,b) ((a > b) ? (a) : (b))

typedef void (*inverse_complex_ssht)(complex double *, const complex double *, int, int, int, ssht_dl_method_t, int);
typedef void (*inverse_real_ssht)(double *, const complex double *, int, int, ssht_dl_method_t, int);
typedef void (*forward_complex_ssht)(complex double *, const complex double *, int, int, int, ssht_dl_method_t, int);
typedef void (*forward_real_ssht)(complex double *, const double *, int, int, ssht_dl_method_t, int);

void so3_adjoint_inverse_direct(
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
    if (verbosity > 0)
    {
        printf("%sComputing adjoint inverse transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_adjoint_mw_inverse_direct with storage method %d...\n"
                    , SO3_PROMPT
                    , storage);
    }    int m_stride = 2*L-1;
    int m_offset = L-1;
    // unused: int n_stride = 2*N-1;
    int n_offset = N-1;
    int mm_stride = 2*L-1;
    int mm_offset = L-1;
    int a_stride = 2*L-1;
    int b_stride = L;
    int bext_stride = 2*L-1;
    // unused: int g_stride = 2*N-1;

    int n_start, n_stop, n_inc;

    switch (n_mode)
    {
    case SO3_N_MODE_ALL:
    case SO3_N_MODE_L:
        n_start = -N+1;
        n_stop  =  N-1;
        n_inc = 1;
        break;
    case SO3_N_MODE_EVEN:
        n_start = ((N-1) % 2 == 0) ? -N+1 : -N+2;
        n_stop  = ((N-1) % 2 == 0) ?  N-1 :  N-2;
        n_inc = 2;
        break;
    case SO3_N_MODE_ODD:
        n_start = ((N-1) % 2 != 0) ? -N+1 : -N+2;
        n_stop  = ((N-1) % 2 != 0) ?  N-1 :  N-2;
        n_inc = 2;
        break;
    case SO3_N_MODE_MAXIMUM:
        n_start = -N+1;
        n_stop  =  N-1;
        n_inc = MAX(1,2*N - 2);
        break;
    default:
        SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    double *sqrt_tbl = calloc(2*(L-1)+2, sizeof(*sqrt_tbl));
    SO3_ERROR_MEM_ALLOC_CHECK(sqrt_tbl);
    double *signs = calloc(L+1, sizeof(*signs));
    SO3_ERROR_MEM_ALLOC_CHECK(signs);
    complex double *exps = calloc(4, sizeof(*exps));
    SO3_ERROR_MEM_ALLOC_CHECK(exps);
    complex double *expsmm = calloc(2*L-1, sizeof(*expsmm));
    SO3_ERROR_MEM_ALLOC_CHECK(expsmm);

    int el, m, n, mm; // mm is for m'
    // Perform precomputations.
    for (el = 0; el <= 2*(L-1)+1; ++el)
        sqrt_tbl[el] = sqrt((double)el);
    for (m = 0; m <= L-1; m += 2)
    {
        signs[m]   =  1.0;
        signs[m+1] = -1.0;
    }
    int i;
    for (i = 0; i < 4; ++i)
        exps[i] = cexp(I*SO3_PION2*i);
    for (mm = -L+1; mm <= L-1; ++mm)
        expsmm[mm + mm_offset] = cexp(-I*mm*SSHT_PI/(2.0*L-1.0));

    // Compute Fourier transform over alpha and gamma, i.e. compute Fmn(b).
    complex double *Fmnb = calloc((2*L-1)*(2*L-1)*(2*N-1), sizeof(*Fmnb));
    SO3_ERROR_MEM_ALLOC_CHECK(Fmnb);
    complex double *inout = calloc((2*L-1)*(2*N-1), sizeof(*inout));
    SO3_ERROR_MEM_ALLOC_CHECK(inout);
    fftw_plan plan = fftw_plan_dft_2d(
                        2*N-1, 2*L-1,
                        inout, inout, 
                        FFTW_FORWARD, 
                        FFTW_ESTIMATE);

    int b, g;
    for (b = 0; b < L; ++b)
    {
        // TODO: This memcpy loop could probably be avoided by using
        // a more elaborate FFTW plan which performs the FFT directly
        // over the 1st and 3rd dimensions of f.
        for (g = 0; g < 2*N-1; ++g)
            memcpy(
                inout + g*a_stride, 
                f + 0 + a_stride*(
                    b + b_stride*(
                    g)), 
                a_stride*sizeof(*f));
        fftw_execute_dft(plan, inout, inout);

        // Apply spatial shift and normalisation factor
        for (n = n_start; n <= n_stop; n += n_inc)
        {
            int n_shift = n < 0 ? 2*N-1 : 0;
            for (m = -L+1; m <= L-1; ++m)
            {
                int m_shift = m < 0 ? 2*L-1 : 0;
                Fmnb[b + bext_stride*(
                     m + m_offset + m_stride*(
                     n + n_offset))] =
                    inout[m + m_shift + m_stride*(
                          n + n_shift)];
            }
        }
    }
    fftw_destroy_plan(plan);

    // Extend Fmnb by filling it with zeroes.
    for (n = n_start; n <= n_stop; n += n_inc)
        for (m = -L+1; m <= L-1; ++m)
            for (b = L; b < 2*L-1; ++b)
                Fmnb[b + bext_stride*(
                     m + m_offset + m_stride*(
                     n + n_offset))] = 0.0;


    // Compute Fourier transform over beta, i.e. compute Fmnm'.
    complex double *Fmnm = calloc((2*L-1)*(2*L-1)*(2*N-1), sizeof(*Fmnm));
    SO3_ERROR_MEM_ALLOC_CHECK(Fmnm);

    plan = fftw_plan_dft_1d(
            2*L-1,
            inout, inout, 
            FFTW_FORWARD, 
            FFTW_ESTIMATE);
    for (n = n_start; n <= n_stop; n += n_inc)
        for (m = -L+1; m <= L-1; ++m)
        {
            memcpy(inout, 
                   Fmnb + 0 + bext_stride*(
                          m + m_offset + m_stride*(
                          n + n_offset)), 
                   bext_stride*sizeof(*Fmnb));
            fftw_execute_dft(plan, inout, inout);

            // Apply spatial shift
            for (mm = -L+1; mm <= L-1; ++mm)
            {
                int mm_shift = mm < 0 ? 2*L-1 : 0;
                Fmnm[m + m_offset + m_stride*(
                     mm + mm_offset + mm_stride*(
                     n + n_offset))] =
                    inout[mm + mm_shift];
            }
        }
    fftw_destroy_plan(plan);
    free(inout);

    // Apply phase modulation to account for sampling offset.
    for (n = n_start; n <= n_stop; n += n_inc)
        for (mm = -L+1; mm <= L-1; ++mm)
            for (m = -L+1; m <= L-1; ++m)
                Fmnm[m + m_offset + m_stride*(
                     mm + mm_offset + mm_stride*(
                     n + n_offset))] *=
                    expsmm[mm + mm_offset];

    // Compute flmn.
    double *dl, *dl8 = NULL;
    dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
    SO3_ERROR_MEM_ALLOC_CHECK(dl);
    if (dl_method == SSHT_DL_RISBO)
    {
        dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
        SO3_ERROR_MEM_ALLOC_CHECK(dl8);
    }
    int dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
    int dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
    for (n = -N+1; n <= N-1; ++n)
        for (el = abs(n); el < L; ++el)
            for (m = -el; m <= el; ++m)
            {
                int ind;
                so3_sampling_elmn2ind(&ind, el, m, n, parameters);
                flmn[ind] = 0.0;
            }

    for (el = L0; el < L; ++el)
    {
        int eltmp;

        // Compute Wigner plane.
        switch (dl_method)
        {
        case SSHT_DL_RISBO:
            if (el != 0 && el == L0)
            {
                for(eltmp = 0; eltmp <= L0; ++eltmp)
                    ssht_dl_beta_risbo_eighth_table(dl8, SO3_PION2, L,
                        SSHT_DL_QUARTER_EXTENDED,
                        eltmp, sqrt_tbl, signs);
                ssht_dl_beta_risbo_fill_eighth2quarter_table(dl,
                    dl8, L,
                    SSHT_DL_QUARTER,
                    SSHT_DL_QUARTER_EXTENDED,
                    el,
                    signs);
            }
            else
            {
                ssht_dl_beta_risbo_eighth_table(dl8, SO3_PION2, L,
                    SSHT_DL_QUARTER_EXTENDED,
                    el, sqrt_tbl, signs);
                ssht_dl_beta_risbo_fill_eighth2quarter_table(dl,
                    dl8, L,
                    SSHT_DL_QUARTER,
                    SSHT_DL_QUARTER_EXTENDED,
                    el,
                    signs);
            }
            break;

        case SSHT_DL_TRAPANI:
            if (el != 0 && el == L0)
            {
                for(eltmp = 0; eltmp <= L0; ++eltmp)
                    ssht_dl_halfpi_trapani_eighth_table(dl, L,
                        SSHT_DL_QUARTER,
                        eltmp, sqrt_tbl);
                ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
                    SSHT_DL_QUARTER,
                    el, signs);
            }
            else
            {
                ssht_dl_halfpi_trapani_eighth_table(dl, L,
                    SSHT_DL_QUARTER,
                    el, sqrt_tbl);
                ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
                    SSHT_DL_QUARTER,
                    el, signs);
            }
            break;

        default:
            SO3_ERROR_GENERIC("Invalid dl method");
        }

        // Compute flmn for current el.

        switch (n_mode)
        {
        case SO3_N_MODE_ALL:
            n_start = MAX(-N+1,-el);
            n_stop  = MIN( N-1, el);
            n_inc = 1;
            break;
        case SO3_N_MODE_EVEN:
            n_start = MAX(-N+1,-el);
            n_start += (-n_start)%2;
            n_stop  = MIN( N-1, el);
            n_stop  -=  n_stop%2;
            n_inc = 2;
            break;
        case SO3_N_MODE_ODD:
            n_start = MAX(-N+1,-el);
            n_start += 1+n_start%2;
            n_stop  = MIN( N-1, el);
            n_stop  -= 1-n_stop%2;
            n_inc = 2;
            break;
        case SO3_N_MODE_MAXIMUM:
            if (el < N-1)
                continue;
            n_start = -N+1;
            n_stop  =  N-1;
            n_inc = MAX(1,2*N-2);
            break;
        case SO3_N_MODE_L:
            if (el >= N)
                continue;
            n_start = -el;
            n_stop  =  el;
            n_inc = MAX(1,2*el);
            break;
        default:
            SO3_ERROR_GENERIC("Invalid n-mode.");
        }

        double elfactor = (2.0*el+1.0)/(8.0*SO3_PI*SO3_PI);
    
        // TODO: Pull out a few multiplications into precomputations
        // or split up loops to avoid conditionals to check signs.
        for (mm = -el; mm <= el; ++mm)
        {
            // These signs are needed for the symmetry relations of
            // Wigner symbols.
            double elmmsign = signs[el] * signs[abs(mm)];

            for (n = n_start; n <= n_stop; n += n_inc)
            {
                double mmsign = mm >= 0 ? 1.0 : signs[el] * signs[abs(n)];
                double elnsign = n >= 0 ? 1.0 : elmmsign;
                 
                // Factor which does not depend on m.
                double elnmm_factor = mmsign * elnsign * elfactor
                                      * dl[abs(n) + dl_offset + abs(mm)*dl_stride];
                
                for (m = -el; m <= el; ++m)
                {
                    mmsign = mm >= 0 ? 1.0 : signs[el] * signs[abs(m)];
                    double elmsign = m >= 0 ? 1.0 : elmmsign;
                    int ind;
                    so3_sampling_elmn2ind(&ind, el, m, n, parameters);
                    int mod = ((m-n)%4 + 4)%4;
                    flmn[ind] += 
                        exps[mod]
                        * elnmm_factor
                        * mmsign * elmsign
                        * dl[abs(m) + dl_offset + abs(mm)*dl_stride]
                        * Fmnm[m + m_offset + m_stride*(
                               mm + mm_offset + mm_stride*(
                               n + n_offset))];

                }
            }
        }
    }

    free(dl);
    if (dl_method == SSHT_DL_RISBO)
        free(dl8);
    free(Fmnb);
    free(Fmnm);
    free(sqrt_tbl);
    free(signs);
    free(exps);
    free(expsmm);

    if (verbosity > 0)
        printf("%sForward transform computed!\n", SO3_PROMPT);
}

void so3_adjoint_inverse_direct_real(
    complex double *flmn, const double *f,
    const so3_parameters_t *parameters
) { }
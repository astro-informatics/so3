// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details

/*! 
 * \file so3_core.c
 * Core algorithms to perform Wigner transform on the rotation group SO(§).
 *
 * \author Jason McEwen
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

/*!  
 * Compute inverse transform for MW method via SSHT.
 *
 * \param[out] f Function on sphere. Provide a buffer of size (2*L-1)*L*(2*N-1).
 * \param[in] flmn Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] M Azimuthal band-limit.
 * \param[in] N Orientational band-limit.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_core_mw_inverse_via_ssht(complex double *f, const complex double *flmn,
	int L, int N,
        so3_storage_t storage,
	int verbosity)
{
    int i, el, n;
    int a, b, g;
    int ind;
    complex double *fn, *flm;
    int offset;
    int fn_n_stride;
    int fftw_rank, fftw_howmany;
    int fftw_idist, fftw_odist;
    int fftw_istride, fftw_ostride;
    int fftw_n;
    fftw_plan plan;

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

    fn_n_stride = L * (2*L-1);

    flm = (complex double*)calloc(L*L, sizeof(complex double));
    // We make fn large enough so that the components with n < 0 can fit in twice.
    // This is utilized later to do a spatial shift via a simple memcpy().
    fn = (complex double*)malloc((2*N-1)*fn_n_stride * sizeof(complex double));

    // Initialize fftw_plan first. With FFTW_ESTIMATE this is technically not
    // necessary but still good practice.
    fftw_rank = 1;
    fftw_n = 2*N-1;
    fftw_howmany = fn_n_stride;
    fftw_idist = fftw_odist = 1;
    fftw_istride = fftw_ostride = fn_n_stride;

    plan = fftw_plan_many_dft(
            fftw_rank, &fftw_n, fftw_howmany,
            fn, NULL, fftw_istride, fftw_idist,
            f, NULL, fftw_ostride, fftw_odist,
            FFTW_BACKWARD, FFTW_ESTIMATE
    );

    for(n = -N+1; n < N; ++n)
    {
        so3_sampling_elmn2ind(&ind, 0, 0, n, L, N, SO3_STORE_ZERO_FIRST_PAD);
        memcpy(flm, flmn + ind, L*L * sizeof(complex double));
        
        offset = 0;
        i = 0;
        for(el = 0; el < L; ++el)
        {
            for (; i < offset + 2*el+1; ++i)
                flm[i] *= sqrt((double)(2*el+1)/(16.*pow(SO3_PI, 3.)));

            offset = i;
        }

        // The conditional applies the spatial transform, so that we store
        // the results in n-order 0, 1, 2, -2, -1
        offset = (n < 0 ? n + 2*N-1 : n);
        ssht_core_mw_inverse_sov_sym(fn + offset*fn_n_stride, flm, L, -n, SSHT_DL_TRAPANI, verbosity);
        
        if(n % 2)
            for(i = 0; i < fn_n_stride; ++i)
                fn[offset*fn_n_stride + i] = -fn[offset*fn_n_stride + i]; 

        printf("\n");
    }

    free(flm);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    free(fn);

    if (verbosity > 0)
        printf("%sInverse transform computed!\n", SO3_PROMPT);
}

    

void so3_core_mw_forward_via_ssht(complex double *flmn, const complex double *f,
	int L, int N,
        so3_storage_t storage,
	int verbosity)
{

}

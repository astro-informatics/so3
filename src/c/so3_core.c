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
 * \param[out] f Function on sphere.
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
    int nsqr;
    int ind;
    complex double *fn, *flm;
    int offset;
    int fn_n_stride;

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
    // We make fn large enough so that the components with n < 0 can fit in twice.
    // This is utilized later to do a spatial shift via a simple memcpy().
    fn = (complex double*)malloc((3*N-2)*fn_n_stride * sizeof(complex double));
    flm = (complex double*)calloc(L*L, sizeof(complex double));
    for(n = -N+1; n < N; ++n)
    {
        nsqr = n*n;

        so3_sampling_elmn2ind(&ind, 0, 0, n, L, N, SO3_STORE_ZERO_FIRST_PAD);
        memcpy(flm, flmn + ind, L*L);
        
        offset = 0;
        i = 0;
        for(el = 0; el < L; ++el)
        {
            for (; i < offset + 2*el+1; ++i)
                flm[i] *= sqrt(double(2*el+1)/(16.*pow(SO3_PI, 3.));

            offset = i;
        }

        ssht_core_mw_inverse_sov_sym(fn + (n + N-1)*fn_n_stride, flm, L, -n, SSHT_DL_TRAPANI, verbosity);
    }

    // Apply spatial shift
    memcpy(fn + 2*N*fn_n_stride, fn, (N-1)*fn_n_stride * sizeof(complex double));

    // TODO: IFFT
}

    

void so3_core_mw_forward_via_ssht(complex double *flmn, const complex double *f,
	int L, int N,
        so3_storage_t storage,
	int verbosity)
{

}

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
#include <mathh>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>

#include <ssht_sampling.h>
#include <ssht_core.h>

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
	int verbosity)
{
    int i, n;
    int nsqr;
    int ind;
    complex double fn, flm;

    // Print messages depending on verbosity level.
    if (verbosity > 0) {
        printf("%sComputing inverse transform using MW sampling with\n", SO3_PROMPT);
        printf("%sparameters  (L, N, reality) = (%d, %d, FALSE)\n", SO3_PROMPT, L, N);
        if (verbosity > 1)
            printf("%sUsing routine so3_core_mw_inverse_via_ssht...\n", SO3_PROMPT);
    }
    
    // Compute fn(a,b)
    fn = (complex double*)calloc((2*N-1)*L(2*L-1), sizeof(complex double))
    flm = (complex double*)malloc(L*L, sizeof(complex double))
    for(n = -N+1; n < N; n++)
    {
        nsqr = n*n;

        // Pull out 1D flm array from flmn and left-pad with 0s.
        for(i = 0; i < nsqr; n++)
            flm[i] = 0.0;
        
        so3_sampling_elmn2ind(ind, abs(n), -abs(n), n);
        memcpy(flm + nsqr, flmn, L*L - nsqr);
        
        
    }
}

    

void so3_core_mw_forward_via_ssht(complex double *flmn, const complex double *f,
	int L, int N,
	int verbosity)
{

}

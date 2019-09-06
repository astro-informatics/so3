// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner, Jason McEwen and Christopher Wallis
// See LICENSE.txt for license details

/*!
 * \file so3_conv.c
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

#include "so3_types.h"
#include "so3_error.h"
#include "so3_sampling.h"

/*!
 * Compute the convolution of one signal with another in harmonic space
 * h = f * g     
 * 
 * \param[out] hlmn Harmonic coefficients. Provide a buffer of size (L*L*(2*N-1).
 * \param[in]  flmn Harmonic coefficients.
 * \param[in]  glmn Harmonic coefficients.
 * \param[in]  parameters A fully populated parameters object.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_conv_harmonic_convolution(
    SO3_COMPLEX(double) * hlmn, const SO3_COMPLEX(double) * flmn,
    const SO3_COMPLEX(double) * glmn, const so3_parameters_t* parameters
)
{
    int ind, el, m, n;
    int L0, L, N;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;

    for (el = 0; el<L0; el++)
    {
        ind = 1;
        hlmn[ind];
    }

}

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

double ran2_dp(int idum);
void so3_test_gen_flmn_complex(complex double *flmn, const so3_parameters_t *parameters, int seed);
void so3_test_gen_flmn_real(complex double *flmn, const so3_parameters_t *parameters, int seed);

#define MAX(a,b) ((a > b) ? (a) : (b))


/*!
 * Generate random Wigner coefficients of a complex signal.
 *
 * \param[out] flmn Random spherical harmonic coefficients generated. Provide
 *                  enough memory for fully padded storage, i.e. (2*N-1)*L*L
 *                  elements. Unused trailing elements will be set to zero.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        L0, L, N, storage, n_mode
 *                        The reality flag is ignored. Use so3_test_gen_flmn_real
 *                        instead for real signals.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_flmn_complex(
    complex double *flmn,
    const so3_parameters_t *parameters,
    int seed)
{
    int L0, L, N;
    int i, el, m, n, ind;
    int n_start, n_stop, n_inc;
    int el_start, el_stop, el_inc;
    int m_start, m_stop, m_inc;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;

    for (i = 0; i < (2*N-1)*L*L; ++i)
        flmn[i] = 0.0;
    
    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, parameters);
    for (n = n_start; n <= n_stop; n += n_inc)
    {
        so3_sampling_el_loop_values(&el_start, &el_stop, &el_inc, n, parameters);
        for (el = el_start; el <= el_stop; el +=el_inc)
        {
            so3_sampling_m_loop_values(&m_start, &m_stop, &m_inc, el);
            for (m = m_start; m <= m_stop; m +=m_inc)
            {
                so3_sampling_elmn2ind(&ind, el, m, n, parameters);
                flmn[ind] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
            }
        }
    }
}




/*!
 * Generate random Wigner coefficients of a real signal. We only generate
 * coefficients for n >= 0, and for n = 0, we need flm0* = (-1)^(m)*fl-m0, so that
 * fl00 has to be real.
 *
 * \param[out] flmn Random spherical harmonic coefficients generated. Provide
 *                  enough memory for fully padded, complex (!) storage, i.e.
 *                  (2*N-1)*L*L elements. Unused trailing elements will be set to
 *                  zero.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        L0, L, N, storage, n_mode
 *                        The reality flag is ignored. Use so3_test_gen_flmn_complex
 *                        instead for complex signals.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_flmn_real(
    complex double *flmn,
    const so3_parameters_t *parameters,
    int seed)
{
    int L0, L, N;
    int i, el, m, n, n_start, n_stop, n_inc, ind;
    double real, imag;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;

    for (i = 0; i < (2*N-1)*L*L; ++i)
        flmn[i] = 0.0;

    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, parameters);

    for (n = n_start; n <= n_stop; n += n_inc)
    {
        if (n == 0)
        {
            // Fill fl00 with random real values
            for (el = L0; el < L; ++el)
            {
                if (parameters->n_mode == SO3_N_MODE_L && el != 0)
                    break;

                so3_sampling_elmn2ind_real(&ind, el, 0, 0, parameters);
                flmn[ind] = (2.0*ran2_dp(seed) - 1.0);
            }

            // Fill fl+-m0 with conjugated random values

            for (el = L0; el < L; ++el)
            {
                if (parameters->n_mode == SO3_N_MODE_L && el != 0)
                    break;

                for (m = 1; m <= el; ++m)
                {
                    real = (2.0*ran2_dp(seed) - 1.0);
                    imag = (2.0*ran2_dp(seed) - 1.0);
                    so3_sampling_elmn2ind_real(&ind, el, m, 0, parameters);
                    flmn[ind] = real + imag * I;
                    so3_sampling_elmn2ind_real(&ind, el, -m, 0, parameters);
                    flmn[ind] = real - imag * I;
                    if (m % 2)
                        flmn[ind] = - real + imag * I;
                    else
                        flmn[ind] = real - imag * I;
                }
            }
        }
        else
        {
            for (el = MAX(L0, n); el < L; ++el)
            {
                if (parameters->n_mode == SO3_N_MODE_L && el != n)
                    break;

                for (m = -el; m <= el; ++m)
                {

                    so3_sampling_elmn2ind_real(&ind, el, m, n, parameters);
                    flmn[ind] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
                }
            }
        }
    }
}




/*!
 * Generate uniform deviate in range [0,1) given seed. (Using double
 * precision.)
 *
 * \note Uniform deviate (Num rec 1992, chap 7.1), original routine
 * said to be 'perfect'.
 *
 * \param[in] idum Seed.
 * \retval ran_dp Generated uniform deviate.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ran2_dp(int idum)
{
    int IM1=2147483563,IM2=2147483399,IMM1=IM1-1,
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
        NTAB=32,NDIV=1+IMM1/NTAB;

    double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;
    int j,k;
    static int iv[32],iy,idum2 = 123456789;
    // N.B. in C static variables are initialised to 0 by default.

    if (idum <= 0) {
        idum= (-idum>1 ? -idum : 1); // max(-idum,1);
        idum2=idum;
        for(j=NTAB+8;j>=1;j--) {
            k=idum/IQ1;
            idum=IA1*(idum-k*IQ1)-k*IR1;
            if (idum < 0) idum=idum+IM1;
            if (j < NTAB) iv[j-1]=idum;
        }
        iy=iv[0];
    }
    k=idum/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum=idum+IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2=idum2+IM2;
    j=1+iy/NDIV;
    iy=iv[j-1]-idum2;
    iv[j-1]=idum;
    if(iy < 1)iy=iy+IMM1;
    return (AM*iy < RNMX ? AM*iy : RNMX); // min(AM*iy,RNMX);
}

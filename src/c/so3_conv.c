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

#include "so3/so3_types.h"
#include "so3/so3_error.h"
#include "so3/so3_sampling.h"
#include "so3/so3_core.h"
#include <ssht/ssht.h>

#define MAX(a,b) ((a > b) ? (a) : (b))
#define MIN(a,b) ((a < b) ? (a) : (b))


/*!
 * Compute the convolution of one signal with another in harmonic space
 * h = f (*) g     
 * 
 * \param[out] hlmn Harmonic coefficients. Provide a buffer of size (L*L*(2*N-1).
 * \param[in]  parameters of hlmn A fully populated parameters object.
 * \param[in]  flmn Harmonic coefficients.
 * \param[in]  parameters of flmn A fully populated parameters object.
 * \param[in]  glmn Harmonic coefficients.
 * \param[in]  parameters of glmn A fully populated parameters object.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_conv_harmonic_convolution(
    SO3_COMPLEX(double) * hlmn, 
    const so3_parameters_t* h_parameters,
    const SO3_COMPLEX(double) * flmn,
    const so3_parameters_t* f_parameters,
    const SO3_COMPLEX(double) * glmn, 
    const so3_parameters_t* g_parameters
)
{
    int ind_h, ind_f, ind_g, el, m, n, k;
    int n_start, n_stop, n_inc;
    int nf_start, nf_stop, nf_inc;
    int el_start, el_stop, el_inc;
    int m_start, m_stop, m_inc;

    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, h_parameters);
    for (n = n_start; n <= n_stop; n += n_inc)
    {
        so3_sampling_el_loop_values(&el_start, &el_stop, &el_inc, n, h_parameters);
        for (el = el_start; el <= el_stop; el +=el_inc)
        {
            so3_sampling_m_loop_values(&m_start, &m_stop, &m_inc, el);
            for (m = m_start; m <= m_stop; m +=m_inc)
            {
                if (h_parameters->reality) so3_sampling_elmn2ind_real(&ind_h, el, m, n, h_parameters);
                else so3_sampling_elmn2ind(&ind_h, el, m, n, h_parameters);
                hlmn[ind_h] = 0;

                so3_sampling_n_loop_values(&nf_start, &nf_stop, &nf_inc, f_parameters);
                for (k = nf_start; k <= nf_stop; k += nf_inc)
                {
                    if (so3_sampling_is_elmn_non_zero(el, n, k, g_parameters))
                    {
                        if (f_parameters->reality) so3_sampling_elmn2ind_real(&ind_f, el, m, k, f_parameters);
                        else so3_sampling_elmn2ind(&ind_f, el, m, k, f_parameters);
                        
                        if (g_parameters->reality) so3_sampling_elmn2ind_real(&ind_g, el, n, k, g_parameters);
                        else so3_sampling_elmn2ind(&ind_g, el, n, k, g_parameters);
                        
                        hlmn[ind_h] += flmn[ind_f] * conj(glmn[ind_g]);
                    }
                }
            }
        }
    }

}

/*!
 * Compute the parameters for the hlmn for the convolution of one signal with another in harmonic space
 * h = f (*) g     
 * 
 * \param[in]  parameters for glmn Harmonic coefficients.
 * \param[out] parameters for hlmn Harmonic coefficients
 * \param[in]  parameters for flmn Harmonic coefficients.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
so3_parameters_t so3_conv_get_parameters_of_convolved_lmn(
    const so3_parameters_t* f_parameters,
    const so3_parameters_t* g_parameters
)
{
    so3_parameters_t h_parameters;
    if (f_parameters->sampling_scheme != g_parameters->sampling_scheme) SO3_ERROR_GENERIC("f and g must have the same sampling_scheme");
    if (f_parameters->n_order != g_parameters->n_order) SO3_ERROR_GENERIC("f and g must have the same n_order");
    if (f_parameters->storage != g_parameters->storage) SO3_ERROR_GENERIC("f and g must have the same storage");
    if (f_parameters->n_mode != g_parameters->n_mode) SO3_ERROR_GENERIC("f and g must have the same n_mode");
    if (f_parameters->reality != g_parameters->reality) SO3_ERROR_GENERIC("f and g must have the same reality");

    h_parameters = *f_parameters;

    h_parameters.L = MIN(f_parameters->L, g_parameters->L);
    h_parameters.N = h_parameters.L;
    h_parameters.L0 = MAX(f_parameters->L0, g_parameters->L0);
    return h_parameters;
}


void so3_conv_get_parameters_of_convolved_lmn_void(
    so3_parameters_t* h_parameters,
    const so3_parameters_t* f_parameters,
    const so3_parameters_t* g_parameters
)
{
    so3_parameters_t dummy = {}, dummy1 = {}, dummy2 = {} ;
    *h_parameters = so3_conv_get_parameters_of_convolved_lmn(
        f_parameters,
        g_parameters
    );
}

/*!
 * Compute the convolution of one signal with another in real space
 * by doing the convolution in harmonic space
 * h = f (*) g     
 * 
 * \param[out] h  real space (alpha, beta, gamma)
 * \param[in]  parameters of h A fully populated parameters object.
 * \param[in]  f  real space (alpha, beta, gamma)
 * \param[in]  parameters of f A fully populated parameters object.
 * \param[in]  g  real space (alpha, beta, gamma)
 * \param[in]  parameters of g A fully populated parameters object.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_conv_convolution(
    SO3_COMPLEX(double) *h,
    const so3_parameters_t *h_parameters,
    const SO3_COMPLEX(double) *f,
    const so3_parameters_t *f_parameters,
    const SO3_COMPLEX(double) *g,
    const so3_parameters_t *g_parameters
)
{
    SO3_COMPLEX(double) *hlmn, *flmn, *glmn;
    // declare hlmn, flmn, glmn
    int hlmn_length = so3_sampling_flmn_size(h_parameters);
    hlmn = malloc(hlmn_length * sizeof *hlmn); SO3_ERROR_MEM_ALLOC_CHECK(hlmn);

    int flmn_length = so3_sampling_flmn_size(f_parameters);
    flmn = malloc(flmn_length * sizeof *flmn); SO3_ERROR_MEM_ALLOC_CHECK(flmn);

    int glmn_length = so3_sampling_flmn_size(g_parameters);
    glmn = malloc(glmn_length * sizeof *glmn); SO3_ERROR_MEM_ALLOC_CHECK(glmn);

    // harmonic transform of f and g
    so3_core_forward_direct(flmn, f, f_parameters);
    so3_core_forward_direct(glmn, g, g_parameters);

    // calculate hlmn
    so3_conv_harmonic_convolution(hlmn, h_parameters, flmn, f_parameters, glmn, g_parameters);
    
    // transform to h
    so3_core_inverse_direct(h, hlmn, h_parameters);

}

void so3_conv_s2toso3_harmonic_convolution(
    SO3_COMPLEX(double) * hlmn, 
    const so3_parameters_t* h_parameters,
    const SO3_COMPLEX(double) * flm,
    const SO3_COMPLEX(double) * glm
)
{
    int hlmn_length = so3_sampling_flmn_size(h_parameters);
    int el, m, n;
    int ind_f, ind_g;
    SO3_COMPLEX(double) psi;

    for (int i=0; i<hlmn_length; i++)
    {
        if (h_parameters->reality) so3_sampling_ind2elmn_real(&el, &m, &n, i, h_parameters);
        else so3_sampling_ind2elmn(&el, &m, &n, i, h_parameters);

        if (abs(m) <= el & abs(n) <= el)
        {
            ssht_sampling_elm2ind(&ind_f, el, m);
            ssht_sampling_elm2ind(&ind_g, el, n);
            psi = 8*SO3_PI*SO3_PI/(2*el+1);
            hlmn[i] = flm[ind_f] * conj(glm[ind_g]) * psi;
        }
        else 
        {
            hlmn[i] = 0;
        }
    }
}

// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h> 
#include "so3/so3_error.h"
#include "so3/so3_types.h"

int so3_sampling_nalpha(const so3_parameters_t *);
int so3_sampling_nbeta(const so3_parameters_t *);
int so3_sampling_ngamma(const so3_parameters_t *);

#define MAX(a,b) ((a > b) ? (a) : (b))

//============================================================================
// Sampling weights
//============================================================================

/*!
 * Compute conjugate weights for toroidal extension.
 *
 * \param[in] parameters A parameters object with (at least)
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \param[in] p Integer index to compute weight for.
 * \retval Corresponding conjugate weight.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
complex double so3_sampling_weight(
    const so3_parameters_t *parameters, 
    int p)
{
    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
    case SO3_SAMPLING_MW_SS:
        if (p == 1)
            return -I * SSHT_PION2;
        else if (p == -1)
            return I * SSHT_PION2;
        else if (p % 2 == 0)
            return 2.0 / (1.0 - p*p);
        else
            return 0.0;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
}

//============================================================================
// Sampling relations for all supported sampling schemes
//============================================================================

/*!
 * Compute size of the signal buffer f to be passed to the forward
 * and returned from the inverse transform function for a given
 * sampling scheme.
 *
 * \note
 *   Computes number of samples on rotation group, *not* over
 *   extended domain.
 * \note
 *   This includes degenerate samples on the poles which are
 *   the same for each value of phi. To get the logical number
 *   of samples (i.e. ignoring without degeneracy) use
 *   \link so3_sampling_n so3_sampling_n\endlink, instead.
 *
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::N N\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval n Number of samples stored in the signal buffers.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_f_size(const so3_parameters_t *parameters)
{
    int L, N;
    L = parameters->L;
    N = parameters->N;

    return so3_sampling_nalpha(parameters) *
           so3_sampling_nbeta(parameters) *
           so3_sampling_ngamma(parameters);
}

/*!
 * Compute total number of samples for a given sampling scheme.
 *
 * \note
 *   Computes number of samples on rotation group, *not* over
 *   extended domain.
 * \note
 *   This returns the logical number of samples (without degeneracy),
 *   and *not* the size of the signal buffer used in the transform
 *   functions. Use \link so3_sampling_f_size so3_sampling_f_size\endlink
 *   to get the actual size of the signal buffer.
 *
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::N N\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval n Number of samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_n(const so3_parameters_t *parameters)
{
    int L, N;
    L = parameters->L;
    N = parameters->N;

    // Are these actually correct?
    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
        return ((2*L-1)*(L-1) + 1)*so3_sampling_ngamma(parameters);
    case SO3_SAMPLING_MW_SS:
        return ((2*L)*(L-1) + 2)*so3_sampling_ngamma(parameters);
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
}


/*!
 * Compute number of alpha samples for a given sampling scheme.
 *
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval nalpha Number of alpha samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_nalpha(const so3_parameters_t *parameters)
{
    int L;
    L = parameters->L;

    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
        return 2*L - 1;
    case SO3_SAMPLING_MW_SS:
        return 2*L;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
}


/*!
 * Compute number of beta samples for a given sampling scheme.
 *
 * \note Computes number of samples in (0,pi], *not* over extended
 * domain.
 *
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval nbeta Number of beta samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_nbeta(const so3_parameters_t *parameters)
{
    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
        return parameters->L;
    case SO3_SAMPLING_MW_SS:
        return parameters->L + 1;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
}

/*!
 * Compute number of gamma samples for a given sampling scheme.
 *
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::N B\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval ngamma Number of gamma samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_ngamma(const so3_parameters_t *parameters)
{
    int N;
    N = parameters->N;

    if (parameters->steerable)
        return N;
    else
        return 2*N-1;
}


/*!
 * Convert alpha index to angle for a given sampling scheme.
 *
 * \note
 *  - a ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] a Alpha index.
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval alpha Alpha angle.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double so3_sampling_a2alpha(int a, const so3_parameters_t *parameters)
{
    int L;
    L = parameters->L;

    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
        return 2.0 * a * SO3_PI / (2.0*L - 1.0);
    case SO3_SAMPLING_MW_SS:
        return 2.0 * a * SO3_PI / (2.0*L);
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
}


/*!
 * Convert beta index to angle for a given sampling scheme.
 *
 * \note
 *  - b ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
 *
 * \param[in] b Beta index.
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval beta Beta angle.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double so3_sampling_b2beta(int b, const so3_parameters_t *parameters)
{
    int L;
    L = parameters->L;

    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
        return (2.0*b + 1.0) * SO3_PI / (2.0*L - 1.0);
    case SO3_SAMPLING_MW_SS:
        return 2.0 * b * SO3_PI / (2.0*L);
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
}


/*!
 * Convert gamma index to angle for a given sampling scheme.
 *
 * \note
 *  - g ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] g Gamma index.
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::N N\endlink,
 *                       \link so3_parameters_t::sampling_scheme sampling_scheme\endlink
 * \retval gamma Gamma angle.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double so3_sampling_g2gamma(int g, const so3_parameters_t *parameters)
{
    int N;
    N = parameters->N;

    if (parameters->steerable)
        return g * SO3_PI / N;
    else
        return 2.0 * g * SO3_PI / (2.0*N - 1.0);
}

//============================================================================
// Harmonic index relations
//============================================================================

/*!
 * Get storage size of flmn array for different storage methods.
 *
 * \param[in] parameters A parameters object with (at least) the following fields:
 *                       \link so3_parameters_t::L L\endlink,
 *                       \link so3_parameters_t::N N\endlink,
 *                       \link so3_parameters_t::storage storage\endlink,
 *                       \link so3_parameters_t::reality reality\endlink
 * \retval Number of coefficients to be stored.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_flmn_size(
    const so3_parameters_t *parameters
) {
    int L, N;
    L = parameters->L;
    N = parameters->N;
    switch (parameters->storage)
    {
    case SO3_STORAGE_PADDED:
        if (parameters->reality)
            return N*L*L;
        else
            return (2*N-1)*L*L;
    case SO3_STORAGE_COMPACT:
        // Both of these are based on the fact that the sum
        // over n*n from 1 to N-1 is (N-1)*N*(2*N-1)/6.
        if (parameters->reality)
            return N*(6*L*L-(N-1)*(2*N-1))/6;
        else
            return (2*N-1)*(3*L*L-N*(N-1))/3;
    default:
        SO3_ERROR_GENERIC("Invalid storage method.");
    }
}

/*!
 * Convert (el,m,n) harmonic indices to 1D index used to access flmn
 * array.
 *
 * \note Index ranges are as follows:
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - n ranges from [-el' .. el'], where el' = min{el, N}
 *  - ind ranges from [0 .. (2*N)(L**2-N(N-1)/3)-1] for compact storage methods
             and from [0 .. (2*N-1)*L**2-1] for 0-padded storage methods.
 *
 * \param[out] ind 1D index to access flmn array.
 * \param[in]  el  Harmonic index.
 * \param[in]  m   Azimuthal harmonic index.
 * \param[in]  n   Orientational harmonic index.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink,
 *                        \link so3_parameters_t::n_order n_order\endlink
 *                        <br>The \link so3_parameters_t::reality reality\endlink
 *                        flag is ignored. Use \link so3_sampling_elmn2ind_real\endlink
 *                        instead.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_elmn2ind(int *ind, int el, int m, int n, const so3_parameters_t *parameters)
{
    int L, N, offset, absn;
    L = parameters->L;
    N = parameters->N;

    // Most of the formulae here are based on the fact that the sum
    // over n*n from 1 to N-1 is (N-1)*N*(2*N-1)/6.
    switch (parameters->storage)
    {
    case SO3_STORAGE_PADDED:
        switch (parameters->n_order)
        {
        case SO3_N_ORDER_ZERO_FIRST:
            offset = ((n < 0) ? -2*n - 1 : 2*n) * L*L;
            *ind = offset + el*el + el + m;
            return;
        case SO3_N_ORDER_NEGATIVE_FIRST:
            offset = (N-1 + n) * L*L;
            *ind = offset + el*el + el + m;
            return;
        default:
            SO3_ERROR_GENERIC("Invalid n-order.");
        }
    case SO3_STORAGE_COMPACT:
        switch (parameters->n_order)
        {
        case SO3_N_ORDER_ZERO_FIRST:
            absn = abs(n);
            if (absn > el)
                SO3_ERROR_GENERIC("Tried to access component with n > l in compact storage.");
            // Initialize offset to the total storage that would be needed if N == n
            offset = (2*absn-1)*(3*L*L - absn*(absn-1))/3;
            // Advance positive n by another lm-chunk
            if (n >= 0)
                offset += L*L - n*n;
            *ind = offset + el*el - n*n + el + m;
            return;
        case SO3_N_ORDER_NEGATIVE_FIRST:
            absn = abs(n);
            if (absn > el)
                SO3_ERROR_GENERIC("Tried to access component with n > l in compact storage.");
            // Initialize offset as for padded storage, minus the correction necessary for n = 0
            offset = (N-1 + n) * L*L - (2*N - 1)*(N-1)*N/6;
            // Now correct the offset for other n due to missing padding
            if (n <= 0)
                offset += absn*(2*absn+1)*(absn+1)/6;
            else
                offset -= absn*(2*absn-1)*(absn-1)/6;
            *ind = offset + el*el - n*n + el + m;
            return;
        default:
            SO3_ERROR_GENERIC("Invalid n-order.");
        }
    default:
        SO3_ERROR_GENERIC("Invalid storage method.");
    }
}

/*!
 * Convert 1D index used to access flmn array to (el,m,n) harmonic
 * indices.
 *
 * \note Index ranges are as follows:
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - n ranges from [-el' .. el'], where el' = min{el, N}
 *  - ind ranges from [0 .. (2*N)(L**2-N(N-1)/3)-1] for compact storage methods
             and from [0 .. (2*N-1)*L**2-1] for 0-padded storage methods.
 *
 * \param[out] el  Harmonic index.
 * \param[out] m   Azimuthal harmonic index.
 * \param[out] n   Orientational harmonic index.
 * \param[in]  ind 1D index to access flm array.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink
 *                        <br>The \link so3_parameters_t::reality reality\endlink
 *                        flag is ignored. Use \link so3_sampling_ind2elmn_real \endlink
 *                        instead.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_ind2elmn(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters)
{
    int L, N, offset;
    L = parameters->L;
    N = parameters->N;

    switch (parameters->storage)
    {
    case SO3_STORAGE_PADDED:
        switch (parameters->n_order)
        {
        case SO3_N_ORDER_ZERO_FIRST:
            *n = ind/(L*L);

            if(*n % 2)
                *n = -(*n+1)/2;
            else
                *n /= 2;

            ind %= L*L;

            *el = sqrt(ind);
            *m = ind - (*el)*(*el) - *el;
            return;
        case SO3_N_ORDER_NEGATIVE_FIRST:
            *n = ind/(L*L) - (N-1);

            ind %= L*L;

            *el = sqrt(ind);
            *m = ind - (*el)*(*el) - *el;
            return;
        default:
            SO3_ERROR_GENERIC("Invalid n-order.");
        }
    case SO3_STORAGE_COMPACT:
        switch (parameters->n_order)
        {
        case SO3_N_ORDER_ZERO_FIRST:
            offset = 0;
            *n = 0;
            // TODO: Can this loop be replaced by an analytical function
            // (or two - one for positive and one for negative *n)
            while(ind + offset >= L*L)
            {
                ind -= L*L - offset;

                if (*n >= 0)
                {
                    *n = -(*n+1);
                    offset = (*n)*(*n);
                }
                else
                {
                    *n = -(*n);
                }
            }

            ind += offset;

            *el = sqrt(ind);
            *m = ind - (*el)*(*el) - *el;
            return;
        case SO3_N_ORDER_NEGATIVE_FIRST:
            *n = -N+1;
            offset = (*n)*(*n);
            // TODO: Can this loop be replaced by an analytical function
            // (or two - one for positive and one for negative *n)
            while(ind + offset >= L*L)
            {
                ind -= L*L - offset;

                (*n)++;
                offset = (*n)*(*n);
            }

            ind += offset;

            *el = sqrt(ind);
            *m = ind - (*el)*(*el) - *el;
            return;
        default:
            SO3_ERROR_GENERIC("Invalid n-order.");
        }
    default:
        SO3_ERROR_GENERIC("Invalid storage method.");
    }
}

/*!
 * Convert (el,m,n) harmonic indices to 1D index used to access flmn
 * array for a real signal.
 *
 * \note Index ranges are as follows:
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - n ranges from [0 .. el'], where el' = min{el, N}
 *  - ind ranges from [0 .. N*(L*L-(N-1)*(2*N-1)/6)-1] for compact storage methods
             and from [0 .. N*L*L-1] for 0-padded storage methods.
 *
 * \param[out] ind 1D index to access flmn array.
 * \param[in]  el  Harmonic index
 * \param[in]  m   Azimuthal harmonic index.
 * \param[in]  n   Orientational harmonic index.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink
 *                        <br>The \link so3_parameters_t::reality reality\endlink
 *                        flag is ignored. Use \link so3_sampling_elmn2ind \endlink
 *                        instead.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_elmn2ind_real(int *ind, int el, int m, int n, const so3_parameters_t *parameters)
{
    int base_ind;
    so3_parameters_t temp_params;

    // Need to make a copy, because subroutines always use
    // NEGATIVE_FIRST storage order.
    temp_params = *parameters;
    temp_params.n_order = SO3_N_ORDER_NEGATIVE_FIRST;

    // TODO: Could be optimized by computing the indices directly.
    switch(parameters->storage)
    {
    case SO3_STORAGE_PADDED:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, &temp_params);
        so3_sampling_elmn2ind(ind, el, m, n, &temp_params);
        (*ind) -= base_ind;
        return;
    case SO3_STORAGE_COMPACT:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, &temp_params);
        so3_sampling_elmn2ind(ind, el, m, n, &temp_params);
        (*ind) -= base_ind;
        return;
    default:
        SO3_ERROR_GENERIC("Invalid storage method.");
    }
}

/*!
 * Convert 1D index used to access flmn array to (el,m,n) harmonic
 * indices for a real signal.
 *
 * \note Index ranges are as follows:
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - n ranges from [0 .. el'], where el' = min{el, N}
 *  - ind ranges from [0 .. N*(L*L-(N-1)*(2*N-1)/6)-1] for compact storage methods
             and from [0 .. N*L*L-1] for 0-padded storage methods.
 *
 * \param[out] el  Harmonic index.
 * \param[out] m   Azimuthal harmonic index.
 * \param[out] n   Orientational harmonic index.
 * \param[in]  ind 1D index to access flm array.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink
 *                        <br>The \link so3_parameters_t::reality reality\endlink
 *                        flag is ignored. Use \link so3_sampling_ind2elmn
 *                        \endlink instead.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_ind2elmn_real(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters)
{
    int base_ind;
    so3_parameters_t temp_params;

    // Need to make a copy, because subroutines always use
    // NEG_FIRST storage.
    temp_params = *parameters;
    temp_params.n_order = SO3_N_ORDER_NEGATIVE_FIRST;

    // TODO: Could be optimized by computing the indices directly.
    switch(parameters->storage)
    {
    case SO3_STORAGE_PADDED:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, &temp_params);
        so3_sampling_ind2elmn(el, m, n, base_ind + ind, &temp_params);
        return;
    case SO3_STORAGE_COMPACT:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, &temp_params);
        so3_sampling_ind2elmn(el, m, n, base_ind + ind, &temp_params);
        return;
    default:
        SO3_ERROR_GENERIC("Invalid storage method.");
    }
}


/*!
 * Calculates the starting, stopping and increment of n 
 * that a applicable for a flmn.
 * This is useful for looping over the whole array
 *
 * \param[out] n_start starting n.
 * \param[out] n_stop  stopping n.
 * \param[out] n_inc   increment in n.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink
 *                        <br>The \link so3_parameters_t::reality reality\endlink
 *                        flag is ignored. Use \link so3_sampling_ind2elmn
 *                        \endlink instead.
 *  example:
 *     so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, parameters);
 *     for (n = n_start; n <= n_stop; n += n_inc) {}
 * 
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_n_loop_values(int *n_start, int *n_stop, int *n_inc, const so3_parameters_t *parameters)
{
    int L0, L, N;

    L0 = parameters->L0;
    L = parameters->L;
    N = parameters->N;

    if (parameters->reality)
    {
        switch (parameters->n_mode)
        {
        case SO3_N_MODE_ALL:
        case SO3_N_MODE_L:
            *n_start = 0;
            *n_stop  = N-1;
            *n_inc = 1;
            break;
        case SO3_N_MODE_EVEN:
            *n_start = 0;
            *n_stop  = ((N-1) % 2 == 0) ?  N-1 :  N-2;
            *n_inc = 2;
            break;
        case SO3_N_MODE_ODD:
            *n_start = 1;
            *n_stop  = ((N-1) % 2 != 0) ?  N-1 :  N-2;
            *n_inc = 2;
            break;
        case SO3_N_MODE_MAXIMUM:
            *n_start = N-1;
            *n_stop  = N-1;
            *n_inc = 1;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid n-mode.");
        }
    }
    else
    {    
        switch (parameters->n_mode)
        {
        case SO3_N_MODE_ALL:
        case SO3_N_MODE_L:
            *n_start = -N+1;
            *n_stop  =  N-1;
            *n_inc = 1;
            break;
        case SO3_N_MODE_EVEN:
            *n_start = ((N-1) % 2 == 0) ? -N+1 : -N+2;
            *n_stop  = ((N-1) % 2 == 0) ?  N-1 :  N-2;
            *n_inc = 2;
            break;
        case SO3_N_MODE_ODD:
            *n_start = ((N-1) % 2 != 0) ? -N+1 : -N+2;
            *n_stop  = ((N-1) % 2 != 0) ?  N-1 :  N-2;
            *n_inc = 2;
            break;
        case SO3_N_MODE_MAXIMUM:
            *n_start = -N+1;
            *n_stop  =  N-1;
            if (N > 1)
                *n_inc = 2*N - 2;
            else
                *n_inc = 1;
            break;
        default:
            SO3_ERROR_GENERIC("Invalid n-mode.");
        }
    }

}


/*!
 * Calculates the starting, stopping and increment of el
 * that a applicable for a flmn.
 * This is useful for looping over the whole array
 *
 * \param[out] el_start starting el.
 * \param[out] el_stop  stopping el.
 * \param[out] el_inc   increment in el.
 * \param[in]  n   The value of n the loop is being applied to.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink
 *                        <br>The \link so3_parameters_t::reality reality\endlink
 *                        flag is ignored. Use \link so3_sampling_ind2elmn
 *                        \endlink instead.
 * \retval none
 *
 * example:
 *      so3_sampling_el_loop_values(&el_start, &el_stop, &el_inc, n, &parameters);
 *      for (el = el_start; el <= el_stop; el +=el_inc){}
 * 
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_el_loop_values(int *el_start, int *el_stop, int *el_inc, const int n, const so3_parameters_t *parameters)
{
    *el_start = MAX(parameters->L0, abs(n));
    if (parameters->n_mode == SO3_N_MODE_L) *el_stop = abs(n);
    else *el_stop = parameters->L-1;
    *el_inc = 1;

}


/*!
 * Calculates the starting, stopping and increment of m
 * that a applicable for a flmn.
 * This is useful for looping over the whole array
 *
 * \param[out] m_start starting m.
 * \param[out] m_stop  stopping m.
 * \param[out] m_inc   increment in m.
 * \param[in]  el   The value of el the loop is being applied to.
 * \retval none
 *
 * example:
 *      so3_sampling_el_loop_values(&m_start, &m_stop, &m_inc, el);
 *      for (m = m_start; m <= m_stop; m +=m_inc){}
 * 
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_sampling_m_loop_values(int *m_start, int *m_stop, int *m_inc, const int el)
{
    *m_start = -el;
    *m_stop = el;
    *m_inc = 1;
}


/*!
 * Queries if integer i would be covered in a loop 
 * for (i=i_start; i<=i_stop; i += i_inc)
 *
 * \param[in]  i       integer to query
 * \param[in]  i_start starting point
 * \param[in]  i_stop  stoping point
 * \param[in]  i_inc   incrementing step

 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
bool so3_sampling_is_i_in_loop_range(const int i, const int i_start, const int i_stop, const int i_inc)
{
    if (i < i_start) return false;
    if (i > i_stop) return false;
    return (i-i_start)%i_inc == 0;
}

/*!
 * Queries if (el,m,n) harmonic indices are non zero in this setting
 *
 * \note Index ranges are as follows:
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - n ranges from [-el' .. el'], where el' = min{el, N}
 *  - ind ranges from [0 .. (2*N)(L**2-N(N-1)/3)-1] for compact storage methods
             and from [0 .. (2*N-1)*L**2-1] for 0-padded storage methods.
 *
 * \param[in]  el  Harmonic index.
 * \param[in]  m   Azimuthal harmonic index.
 * \param[in]  n   Orientational harmonic index.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link so3_parameters_t::L L\endlink,
 *                        \link so3_parameters_t::N N\endlink,
 *                        \link so3_parameters_t::storage storage\endlink,
 *                        \link so3_parameters_t::n_order n_order\endlink,
 *                        \link so3_parameters_t::reality reality\endlink
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
bool so3_sampling_is_elmn_non_zero(const int el, const int m, const int n, const so3_parameters_t *parameters)
{
    int L, N, offset, absn;
    L = parameters->L;
    N = parameters->N;

    int n_start, n_stop, n_inc;
    int el_start, el_stop, el_inc;
    int m_start, m_stop, m_inc;

    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, parameters);
    if (!(so3_sampling_is_i_in_loop_range(n, n_start, n_stop, n_inc))) return false;

    so3_sampling_el_loop_values(&el_start, &el_stop, &el_inc, n, parameters);
    if (!(so3_sampling_is_i_in_loop_range(el, el_start, el_stop, el_inc))) return false;

    so3_sampling_m_loop_values(&m_start, &m_stop, &m_inc, el);
    if (!(so3_sampling_is_i_in_loop_range(m, m_start, m_stop, m_inc))) return false;

    return true;

}

int so3_sampling_is_elmn_non_zero_return_int(const int el, const int m, const int n, const so3_parameters_t *parameters)
    {
        return (int)so3_sampling_is_elmn_non_zero(el, m, n, parameters);
    }

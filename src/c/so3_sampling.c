// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "so3_error.h"
#include "so3_types.h"

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

    // Are these actually correct?
    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
        return (2*L-1)*(L)*(2*N-1);
    case SO3_SAMPLING_MW_SS:
        return (2*L)*(L+1)*(2*N-1);
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
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
        return ((2*L-1)*(L-1) + 1)*(2*N-1);
    case SO3_SAMPLING_MW_SS:
        return ((2*L)*(L-1) + 2)*(2*N-1);
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

    // As long as we don't introduce other sampling schemes with
    // different ngamma, we could actually remove this switch
    // statement.
    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
    case SO3_SAMPLING_MW_SS:
        return 2*N-1;
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
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

    // As long as we don't introduce other sampling schemes with
    // different gamma sampling, we could actually remove this switch
    // statement.
    switch (parameters->sampling_scheme)
    {
    case SO3_SAMPLING_MW:
    case SO3_SAMPLING_MW_SS:
        return 2.0 * g * SO3_PI / (2.0*N - 1.0);
    default:
        SO3_ERROR_GENERIC("Invalid sampling scheme.");
    }
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
inline int so3_sampling_flmn_size(
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
inline void so3_sampling_elmn2ind(int *ind, int el, int m, int n, const so3_parameters_t *parameters)
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
inline void so3_sampling_ind2elmn(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters)
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
inline void so3_sampling_elmn2ind_real(int *ind, int el, int m, int n, const so3_parameters_t *parameters)
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
inline void so3_sampling_ind2elmn_real(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters)
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

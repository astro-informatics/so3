// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "so3_error.h"
#include "so3_types.h"

//============================================================================
// Sampling relations
//============================================================================

/*!
 * Compute total number of samples for McEwen and Wiaux sampling.
 *
 * /note Computes number of samples on rotation group, *not* over
 * extended domain.
 *
 * \param[in] L Harmonic band-limit.
 * \param[in] N Orientational harmonic band-limit.
 * \retval n Number of samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_mw_n(int L, int N)
{
    // Is this actually correct?
    return (2*L-1)*(L-1)*(2*N-1) + 1;
}


/*!
 * Compute number of alpha samples for McEwen and Wiaux sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval nalpha Number of alpha samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_mw_nalpha(int L)
{
    return 2*L-1;
}


/*!
 * Compute number of beta samples for McEwen and Wiaux sampling.
 *
 * /note Computes number of samples in (0,pi], *not* over extended
 * domain.
 *
 * \param[in] L Harmonic band-limit.
 * \retval nbeta Number of beta samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_mw_nbeta(int L)
{
    return L;
}


/*!
 * Compute number of gamma samples for McEwen and Wiaux sampling.
 *
 * \param[in] N Orientational harmonic band-limit.
 * \retval ngamma Number of gamma samples.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int so3_sampling_mw_ngamma(int N)
{
    return 2*N-1;
}


/*!
 * Convert alpha index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - a ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] a Alpha index.
 * \param[in] L Harmonic band-limit.
 * \retval alpha Alpha angle.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double so3_sampling_mw_a2alpha(int a, int L)
{
    return 2.0 * a * SO3_PI / (2.0*L - 1.0);
}


/*!
 * Convert beta index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - b ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
 *
 * \param[in] b Beta index.
 * \param[in] L Harmonic band-limit.
 * \retval beta Beta angle.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double so3_sampling_mw_b2beta(int b, int L)
{
    return (2.0*b + 1.0) * SO3_PI / (2.0*L - 1.0);
}


/*!
 * Convert gamma index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - g ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] g Gamma index.
 * \param[in] N Orientational harmonic band-limit.
 * \retval gamma Gamma angle.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double so3_sampling_mw_g2gamma(int g, int N)
{
    return 2.0 * g * SO3_PI / (2.0*N - 1.0);
}


//============================================================================
// Harmonic index relations
//============================================================================

/*!
 * Get storage size of flmn array for different storage methods.
 *
 * \param[in]  L        Harmonic band-limit.
 * \param[in]  N        Orientational band-limit.
 * \param[in]  storage  Indicates how flm-blocks are stored.
 * \param[in]  real     Pass in a non-zero value to indicate that the
 *                      flmn belong to a real signal.
 *
 * \retval Number of coefficients to be stored or 0 for invalid storage
 *         parameter.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
inline int so3_sampling_flmn_size(
    int L,
    int N,
    so3_storage_t storage,
    int real
) {
    switch (storage)
    {
    case SO3_STORE_ZERO_FIRST_PAD:
    case SO3_STORE_NEG_FIRST_PAD:
        if (real)
            return N*L*L;
        else
            return (2*N-1)*L*L;
    case SO3_STORE_ZERO_FIRST_COMPACT:
    case SO3_STORE_NEG_FIRST_COMPACT:
        if (real)
            return N*(6*L*L-(N-1)*(2*N-1))/6;
        else
            return (2*N-1)*(3*L*L-N*(N-1))/3;
    default:
        return 0;
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
 * \param[out] ind 1D index to access flmn array [output].
 * \param[in]  el  Harmonic index [input].
 * \param[in]  m   Azimuthal harmonic index [input].
 * \param[in]  n   Orientational harmonic index [input].
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
inline void so3_sampling_elmn2ind(int *ind, int el, int m, int n, int L, int N, so3_storage_t storage)
{
    int offset, absn;
    // Most of the formulae here are based on the fact that the sum
    // over n*n from 1 to N-1 is (N-1)*N*(2*N-1)/6.
    switch(storage)
    {
    case SO3_STORE_ZERO_FIRST_PAD:
        offset = ((n < 0) ? -2*n - 1 : 2*n) * L*L;
        *ind = offset + el*el + el + m;
        return;
    case SO3_STORE_ZERO_FIRST_COMPACT:
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
    case SO3_STORE_NEG_FIRST_PAD:
        offset = (N-1 + n) * L*L;
        *ind = offset + el*el + el + m;
        return;
    case SO3_STORE_NEG_FIRST_COMPACT:
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
        // printf("(el,m,n) = (%d,%d,%d) => ind = %d\n", el, m, n, *ind);
        // printf("offset = %d\n", offset);
        return;
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
 * \param[in]  ind 1D index to access flm array [output].
 * \param[out] el  Harmonic index [input].
 * \param[out] m   Azimuthal harmonic index [input].
 * \param[out] n   Orientational harmonic index [input].
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
inline void so3_sampling_ind2elmn(int *el, int *m, int *n, int ind, int L, int N, so3_storage_t storage)
{
    int offset;

    switch (storage)
    {
    case SO3_STORE_ZERO_FIRST_PAD:
        *n = ind/(L*L);

        if(*n % 2)
            *n = -(*n+1)/2;
        else
            *n /= 2;

        ind %= L*L;

        *el = sqrt(ind);
        *m = ind - (*el)*(*el) - *el;
        return;
    case SO3_STORE_ZERO_FIRST_COMPACT:
        offset = 0;
        *n = 0;
        // Can this loop be replaced by an analytical function
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
    case SO3_STORE_NEG_FIRST_PAD:
        *n = ind/(L*L) - (N-1);

        ind %= L*L;

        *el = sqrt(ind);
        *m = ind - (*el)*(*el) - *el;
        return;
    case SO3_STORE_NEG_FIRST_COMPACT:
        *n = -N+1;
        offset = (*n)*(*n);
        // Can this loop be replaced by an analytical function
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
 * \param[out] ind 1D index to access flmn array [output].
 * \param[in]  el  Harmonic index [input].
 * \param[in]  m   Azimuthal harmonic index [input].
 * \param[in]  n   Orientational harmonic index [input].
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
inline void so3_sampling_elmn2ind_real(int *ind, int el, int m, int n, int L, int N, so3_storage_t storage)
{
    int base_ind;

    // TODO: Could be optimized by computing the indices directly.
    switch(storage)
    {
    case SO3_STORE_ZERO_FIRST_PAD:
    case SO3_STORE_NEG_FIRST_PAD:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, L, N, SO3_STORE_NEG_FIRST_PAD);
        so3_sampling_elmn2ind(ind, el, m, n, L, N, SO3_STORE_NEG_FIRST_PAD);
        (*ind) -= base_ind;
        return;
    case SO3_STORE_ZERO_FIRST_COMPACT:
    case SO3_STORE_NEG_FIRST_COMPACT:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, L, N, SO3_STORE_NEG_FIRST_COMPACT);
        so3_sampling_elmn2ind(ind, el, m, n, L, N, SO3_STORE_NEG_FIRST_COMPACT);
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
 * \param[in]  ind 1D index to access flm array [output].
 * \param[out] el  Harmonic index [input].
 * \param[out] m   Azimuthal harmonic index [input].
 * \param[out] n   Orientational harmonic index [input].
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
inline void so3_sampling_ind2elmn_real(int *el, int *m, int *n, int ind, int L, int N, so3_storage_t storage)
{
    int base_ind;

    switch(storage)
    {
    case SO3_STORE_ZERO_FIRST_PAD:
    case SO3_STORE_NEG_FIRST_PAD:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, L, N, SO3_STORE_NEG_FIRST_PAD);
        so3_sampling_ind2elmn(el, m, n, base_ind + ind, L, N, SO3_STORE_NEG_FIRST_PAD);
        return;
    case SO3_STORE_ZERO_FIRST_COMPACT:
    case SO3_STORE_NEG_FIRST_COMPACT:
        so3_sampling_elmn2ind(&base_ind, 0, 0, 0, L, N, SO3_STORE_NEG_FIRST_COMPACT);
        so3_sampling_ind2elmn(el, m, n, base_ind + ind, L, N, SO3_STORE_NEG_FIRST_COMPACT);
        return;
    default:
        SO3_ERROR_GENERIC("Invalid storage method.");
    }
}

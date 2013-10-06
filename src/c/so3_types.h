// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*! \mainpage SO3 C documentation
 *
 * We document the C source code here.  For an example of usage,
 * see the so3_test.c program.
 * For installation instructions, see the general SSHT
 * documentation available
 * <a href="../../index.html">here</a>.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */


/*! \file so3_types.h
 *  Types used in SO3 package.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#ifndef SO3_TYPES
#define SO3_TYPES

#define SO3_PI    3.141592653589793238462643383279502884197

#define SO3_SQRT2 1.41421356237309504880168872420969807856967



#define SO3_PROMPT "[so3] "

typedef enum {
    // order of n in storage is 0, -1, 1, -2, 2, ...
    // for each n, left-pad each lm-array with 0s to L^2 values
    SO3_STORE_ZERO_FIRST_PAD,
    // order of n in storage is 0, -1, 1, -2, 2, ...
    // do not store lm for l < |n|
    SO3_STORE_ZERO_FIRST_COMPACT,
    // order of n in storage is -2, -1, 0, 1, 2, ...
    // for each n, left-pad each lm-array with 0s to L^2 values
    SO3_STORE_NEG_FIRST_PAD,
    // order of n in storage is -2, -1, 0, 1, 2, ...
    // do not store lm for l < |n|
    SO3_STORE_NEG_FIRST_COMPACT,
    // "guard" value that equals the number of usable enum values.
    // useful in loops, for instance.
    SO3_STORE_SIZE
} so3_storage_t;

typedef enum {
    // flmn are potentially non-zero for all values of n
    SO3_N_MODE_ALL,
    // flmn are only non-zero for even n
    SO3_N_MODE_EVEN,
    // flmn are only non-zero for odd n
    SO3_N_MODE_ODD,
    // flmn are only non-zero for |n| = N-1
    SO3_N_MODE_MAXIMUM,
    // "guard" value that equals the number of usable enum values.
    // useful in loops, for instance.
    SO3_N_MODE_SIZE
} so3_n_mode_t;

typedef enum {
    // McEwen and Wiaux sampling:
    // 2*L-1 samples in alpha, in [0, 2pi).
    // L samples in beta, in (0, pi].
    // 2*N-1 samples in gamma, in [0, 2pi).
    SO3_SAMPLING_MW,
    // McEwen and Wiaux symmetric sampling:
    // 2*L samples in alpha, in [0, 2pi).
    // L+1 samples in beta, in [0, pi].
    // 2*N-1 samples in gamma, in [0, 2pi).
    SO3_SAMPLING_MW_SS,
    // "guard" value that equals the number of usable enum values.
    // useful in loops, for instance.
    SO3_SAMPLING_SIZE
} so3_sampling_t;

#endif

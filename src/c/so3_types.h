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
    SO3_STORE_NEG_FIRST_COMPACT
} so3_storage_t;

typedef enum {
    // flmn are potentially non-zero for all values of n
    SO3_N_MODE_ALL,
    // flmn are only non-zero for even n
    SO3_N_MODE_EVEN,
    // flmn are only non-zero for odd n
    SO3_N_MODE_ODD,
    // flmn are only non-zero for |n| = N-1
    SO3_N_MODE_MAXIMUM
} so3_n_mode_t;

#endif


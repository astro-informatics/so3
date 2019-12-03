// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*! \mainpage SO3 C documentation
 *
 * We document the C source code here.  For an example of usage,
 * see the so3_test.c program.
 * For installation instructions, see the general SO3
 * documentation available <a href="../../index.html">here</a>.
 * Make sure to check out the documentation of so3_parameters_t, as
 * it is used across the entire API.
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

#include <ssht/ssht.h>
#ifdef __cplusplus
#include <complex>
#define SO3_COMPLEX(TYPE) std::complex<TYPE>
extern "C" {
#else
#define SO3_COMPLEX(TYPE) TYPE complex
#endif


#ifdef __cplusplus
}
#endif

#define SO3_PI    3.141592653589793238462643383279502884197
#define SO3_PION2 1.570796326794896619231321691639751442099

#define SO3_SQRT2 1.41421356237309504880168872420969807856967

#define SO3_PROMPT "[so3] "

typedef enum {
    /*! order of flm-blocks in terms of n in storage is 0, -1, 1, -2, 2, ... */
    SO3_N_ORDER_ZERO_FIRST,
    /*! order of flm-blocks in terms of n in storage is ... -2, -1, 0, 1, 2, ... */
    SO3_N_ORDER_NEGATIVE_FIRST,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    SO3_N_ORDER_SIZE
} so3_n_order_t;

typedef enum {
    /*! left-pad each flm-array with 0s to L^2 values */
    SO3_STORAGE_PADDED,
    /*! do not store lm for l < |n| */
    SO3_STORAGE_COMPACT,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    SO3_STORAGE_SIZE
} so3_storage_t;

typedef enum {
    /*! flmn are potentially non-zero for all values of n */
    SO3_N_MODE_ALL,
    /*! flmn are only non-zero for even n */
    SO3_N_MODE_EVEN,
    /*! flmn are only non-zero for odd n */
    SO3_N_MODE_ODD,
    /*! flmn are only non-zero for |n| = N-1 */
    SO3_N_MODE_MAXIMUM,
    /*! flmn are only non-zero for |n| = el */
    SO3_N_MODE_L,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    SO3_N_MODE_SIZE
} so3_n_mode_t;

typedef enum {
    /*!
     * McEwen and Wiaux sampling:
     * 2*L-1 samples in alpha, in [0, 2pi).
     * L samples in beta, in (0, pi].
     * 2*N-1 samples in gamma, in [0, 2pi).
     */
    SO3_SAMPLING_MW,
    /*!
     * McEwen and Wiaux symmetric sampling:
     * 2*L samples in alpha, in [0, 2pi).
     * L+1 samples in beta, in [0, pi].
     * 2*N-1 samples in gamma, in [0, 2pi).
     */
    SO3_SAMPLING_MW_SS,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    SO3_SAMPLING_SIZE
} so3_sampling_t;

/*!
 * A struct with all parameters that are common to several
 * functions of the API. In general only one struct needs to
 * be created and a const pointer to it is passed around.
 * \attention
 *   Make sure to use an initializer upon
 *   declaration, even if it is left empty. This
 *   ensures that all fields are initialized to
 *   zero (even if in non-static storage). This way,
 *   your code will remain backwards compatible if
 *   more fields are added to this struct in the
 *   future:
 *   \code{.c}
 *   so3_parameters_t parameters = {};
 *   \endcode
 */
typedef struct {
    /*!
     * Detail level for diagnostic console output in range [0,5].
     * \var int verbosity
     */
    int verbosity;

    /*!
     * A non-zero value indicates that the signal f is real. Not
     * all functions respect this value - instead there may be
     * separate complex and real functions. See the documentation
     * of each function for details.
     */
    int reality;

    /*!
     * Lower harmonic band-limit. flmn with l < L0 will be zero.
     * \var int L0
     */
    int L0;

    /*!
     * Upper harmonic band-limit. Only flmn with l < L will be stored
     * and considered.
     * \var int L
     */
    int L;

    /*!
     * Upper orientational band-limit. Only flmn with n < N will
     * be stored.
     * \var int N
     */
    int N;

    /*!
     * Sampling scheme to use for samples of the signal f.
     * \var so3_sampling_t sampling_scheme
     */
    so3_sampling_t sampling_scheme;

    /*!
     * Indicates the order of n-values by which individual flm-blocks
     * are stored.
     * \var so3_n_order_t n_order
     */
    so3_n_order_t n_order;

    /*!
     * Type of storage (padded or compact).
     * \var so3_storage_t storage
     */
    so3_storage_t storage;

    /*!
     * Indicates if entire blocks of flm for certain values of n are
     * zero.
     * \var so3_n_mode_t n_mode
     */
    so3_n_mode_t n_mode;

    /*!
     * Recursion method to use for computing Wigner functions.
     * \var ssht_dl_method_t dl_method
     */
    ssht_dl_method_t dl_method;

    /*!
     * A non-zero value indicates that the signal is steerable.
     */
    int steerable;
} so3_parameters_t;

#endif

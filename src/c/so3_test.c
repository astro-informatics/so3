// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Buettner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_test.c
 * Applies SO3 algorithms to perform inverse and forward Wigner
 * transforms (respectively) to check that the original
 * signal is reconstructed exactly (to numerical precision). Test is
 * performed on a random signal with harmonic coefficients uniformly
 * sampled from (-1,1).
 *
 * Usage: so3_test [L [N [seed]]], e.g. so3_test 64 4 314
 *
 * Defaults: L = 4, N = 4, seed = 1
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Buettner</a>
           <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>

#include <so3.h>

#define NREPEAT 5
#define MIN(a,b) ((a < b) ? (a) : (b))
#define MAX(a,b) ((a > b) ? (a) : (b))

double get_max_error(complex double *expected, complex double *actual, int n);
double ran2_dp(int idum);
void so3_test_gen_compact_flmn_complex(complex double *flmn, int L, int N, int seed);
void so3_test_gen_padded_zero_first_flmn_complex(complex double *flmn, int L, int N, int seed);
void so3_test_gen_padded_neg_first_flmn_complex(complex double *flmn, int L, int N, int seed);

int main(int argc, char **argv)
{
    int L, N;
    complex double *flmn_compact_orig, *flmn_compact_syn, *flmn_padded_orig, *flmn_padded_syn;
    complex double *f;
    int verbosity = 0;
    int seed;
    clock_t time_start, time_end;
    int i;

    double errors_compact_zero_first[NREPEAT];
    double errors_compact_neg_first[NREPEAT];
    double errors_padded_zero_first[NREPEAT];
    double errors_padded_neg_first[NREPEAT];
    double durations_forward_compact_zero_first[NREPEAT];
    double durations_inverse_compact_zero_first[NREPEAT];
    double durations_forward_compact_neg_first[NREPEAT];
    double durations_inverse_compact_neg_first[NREPEAT];
    double durations_forward_padded_zero_first[NREPEAT];
    double durations_inverse_padded_zero_first[NREPEAT];
    double durations_forward_padded_neg_first[NREPEAT];
    double durations_inverse_padded_neg_first[NREPEAT];

    double total_error;
    double min_time;

    // Parse command line arguments
    L = N = 32;
    seed = 1;
    if (argc > 1)
    {
        L = atoi(argv[1]);
        if (argc > 2)
            N = atoi(argv[2]);
        else
            N = L;
    }

    if (argc > 3)
        seed = atoi(argv[3]);

    flmn_padded_orig = malloc((2*N-1)*L*L * sizeof *flmn_padded_orig);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_padded_orig);
    flmn_padded_syn = malloc((2*N-1)*L*L * sizeof *flmn_padded_syn);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_padded_syn);
    flmn_compact_orig = malloc((2*N-1)*(3*L*L-N*(N-1))/3 * sizeof *flmn_compact_orig);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_compact_orig);
    flmn_compact_syn = malloc((2*N-1)*(3*L*L-N*(N-1))/3 * sizeof *flmn_compact_syn);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_compact_syn);
    f = malloc((2*L-1)*L*(2*N-1) * sizeof *f);
    SO3_ERROR_MEM_ALLOC_CHECK(f);

    // Write program name.
    printf("\n");
    printf("SO3 test program (C implementation)\n");
    printf("================================================================\n");

    for (i = 0; i < NREPEAT; ++i)
    {
        // ===========================================================================================
        // Padded storage with n-order 0, -1, 1, -2, 2, ...
        printf("Testing padded storage with n = 0 first. (run %d)\n", i+1);

        so3_test_gen_padded_zero_first_flmn_complex(flmn_padded_orig, L, N, seed);

        time_start = clock();
        so3_core_mw_inverse_via_ssht(f, flmn_padded_orig, L, N, SO3_STORE_ZERO_FIRST_PAD, verbosity);
        time_end = clock();
        durations_inverse_padded_zero_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        time_start = clock();
        so3_core_mw_forward_via_ssht(flmn_padded_syn, f, L, N, SO3_STORE_ZERO_FIRST_PAD, verbosity);
        time_end = clock();
        durations_forward_padded_zero_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        errors_padded_zero_first[i] = get_max_error(flmn_padded_orig, flmn_padded_syn, (2*N-1)*L*L);

        // ===========================================================================================
        // Padded storage with n-order ..., -2, -1, 0, 1, 2, ...
        printf("Testing padded storage with n = -N+1 first. (run %d)\n", i+1);

        so3_test_gen_padded_neg_first_flmn_complex(flmn_padded_orig, L, N, seed);

        time_start = clock();
        so3_core_mw_inverse_via_ssht(f, flmn_padded_orig, L, N, SO3_STORE_NEG_FIRST_PAD, verbosity);
        time_end = clock();
        durations_inverse_padded_neg_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        time_start = clock();
        so3_core_mw_forward_via_ssht(flmn_padded_syn, f, L, N, SO3_STORE_NEG_FIRST_PAD, verbosity);
        time_end = clock();
        durations_forward_padded_neg_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        errors_padded_neg_first[i] = get_max_error(flmn_padded_orig, flmn_padded_syn, (2*N-1)*L*L);

        // ===========================================================================================
        // Compact storage with n-order 0, -1, 1, -2, 2, ...
        printf("Testing compact storage with n = 0 first. (run %d)\n", i+1);

        so3_test_gen_compact_flmn_complex(flmn_compact_orig, L, N, seed);

        time_start = clock();
        so3_core_mw_inverse_via_ssht(f, flmn_compact_orig, L, N, SO3_STORE_ZERO_FIRST_COMPACT, verbosity);
        time_end = clock();
        durations_inverse_compact_zero_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        time_start = clock();
        so3_core_mw_forward_via_ssht(flmn_compact_syn, f, L, N, SO3_STORE_ZERO_FIRST_COMPACT, verbosity);
        time_end = clock();
        durations_forward_compact_zero_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        errors_compact_zero_first[i] = get_max_error(flmn_compact_orig, flmn_compact_syn, (2*N-1)*(3*L*L-N*(N-1))/3);

        // ===========================================================================================
        // Compact storage with n-order ..., -2, -1, 0, 1, 2, ...
        printf("Testing compact storage with n = -N+1 first. (run %d)\n", i+1);

        so3_test_gen_compact_flmn_complex(flmn_compact_orig, L, N, seed);

        time_start = clock();
        so3_core_mw_inverse_via_ssht(f, flmn_compact_orig, L, N, SO3_STORE_NEG_FIRST_COMPACT, verbosity);
        time_end = clock();
        durations_inverse_compact_neg_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        time_start = clock();
        so3_core_mw_forward_via_ssht(flmn_compact_syn, f, L, N, SO3_STORE_NEG_FIRST_COMPACT, verbosity);
        time_end = clock();
        durations_forward_compact_neg_first[i] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

        errors_compact_neg_first[i] = get_max_error(flmn_compact_orig, flmn_compact_syn, (2*N-1)*(3*L*L-N*(N-1))/3);
    }

    free(flmn_padded_orig);
    free(flmn_padded_syn);
    free(flmn_compact_orig);
    free(flmn_compact_syn);
    free(f);

    // =========================================================================
    // Summarise results

    printf("================================================================\n");
    printf("Summary\n\n");
    printf("Runs   = %40d\n", NREPEAT);
    printf("L      = %40d\n", L);
    printf("N      = %40d\n\n", N);

    printf("Padded storage with n = 0 component first.\n");
    min_time = durations_forward_padded_zero_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_forward_padded_zero_first[i]);
    printf("  Minimum time for forward transform: %fs\n", min_time);
    min_time = durations_inverse_padded_zero_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_inverse_padded_zero_first[i]);
    printf("  Minimum time for inverse transform: %fs\n", min_time);
    total_error = errors_padded_zero_first[0];
    for(i = 1; i < NREPEAT; ++i) total_error += errors_padded_zero_first[i];
    printf("  Average max error for inverse transform: %e\n", total_error/(double)NREPEAT);

    printf("Padded storage with n = -N+1 component first.\n");
    min_time = durations_forward_padded_neg_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_forward_padded_neg_first[i]);
    printf("  Minimum time for forward transform: %fs\n", min_time);
    min_time = durations_inverse_padded_neg_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_inverse_padded_neg_first[i]);
    printf("  Minimum time for inverse transform: %fs\n", min_time);
    total_error = errors_padded_neg_first[0];
    for(i = 1; i < NREPEAT; ++i) total_error += errors_padded_neg_first[i];
    printf("  Average max error for inverse transform: %e\n", total_error/(double)NREPEAT);

    printf("Compact storage with n = 0 component first.\n");
    min_time = durations_forward_compact_zero_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_forward_compact_zero_first[i]);
    printf("  Minimum time for forward transform: %fs\n", min_time);
    min_time = durations_inverse_compact_zero_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_inverse_compact_zero_first[i]);
    printf("  Minimum time for inverse transform: %fs\n", min_time);
    total_error = errors_compact_zero_first[0];
    for(i = 1; i < NREPEAT; ++i) total_error += errors_compact_zero_first[i];
    printf("  Average max error for inverse transform: %e\n", total_error/(double)NREPEAT);

    printf("Compact storage with n = -N+1 component first.\n");
    min_time = durations_forward_compact_neg_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_forward_compact_neg_first[i]);
    printf("  Minimum time for forward transform: %fs\n", min_time);
    min_time = durations_inverse_compact_neg_first[0];
    for(i = 1; i < NREPEAT; ++i) min_time = MIN(min_time, durations_inverse_compact_neg_first[i]);
    printf("  Minimum time for inverse transform: %fs\n", min_time);
    total_error = errors_compact_neg_first[0];
    for(i = 1; i < NREPEAT; ++i) total_error += errors_compact_neg_first[i];
    printf("  Average max error for inverse transform: %e\n", total_error/(double)NREPEAT);


    return 0;
}

double get_max_error(complex double *expected, complex double *actual, int n)
{
    int i;
    double error, maxError = 0;

    for (i = 0; i < n; ++i)
    {
        error = cabs(expected[i] - actual[i]);
        maxError = MAX(error, maxError);
    }

    return maxError;
}

/*!
 * Generate random spherical harmonic coefficients of a complex
 * signal for compact storage methods.
 *
 * \param[out] flmn Random spherical harmonic coefficients generated.
 * \param[in] L Harmonic band-limit.
 * \param[in] N Orientational band-limit.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Buettner</a>
           <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_compact_flmn_complex(complex double *flmn, int L, int N, int seed)
{
    int i;

    for (i = 0; i < (2*N-1)*(3*L*L-N*(N-1))/3; ++i)
        flmn[i]  = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
}

/*!
 * Generate random spherical harmonic coefficients of a complex
 * signal for storage with n = 0 components first.
 *
 * \param[out] flmn Random spherical harmonic coefficients generated.
 * \param[in] L Harmonic band-limit.
 * \param[in] N Orientational band-limit.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Buettner</a>
           <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_padded_zero_first_flmn_complex(complex double *flmn, int L, int N, int seed)
{
    int n, el, m, offset;
    for(n = -N+1; n < N; ++n)
    {
        if (n < 0)
            offset = -2*n-1;
        else
            offset = 2*n;

        for(el = 0; el < L; ++el)
            for(m = -el; m <= el; ++m)
                flmn[offset*L*L + el*el + el + m]  = (el < abs(n)) ? 0.0 : (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
    }
}

/*!
 * Generate random spherical harmonic coefficients of a complex
 * signal for storage with negative n components first.
 *
 * \param[out] flmn Random spherical harmonic coefficients generated.
 * \param[in] L Harmonic band-limit.
 * \param[in] N Orientational band-limit.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Buettner</a>
           <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_padded_neg_first_flmn_complex(complex double *flmn, int L, int N, int seed)
{
    int n, el, m;
    for(n = -N+1; n < N; ++n)
        for(el = 0; el < L; ++el)
            for(m = -el; m <= el; ++m)
                flmn[(n + N-1)*L*L + el*el + el + m]  = (el < abs(n)) ? 0.0 : (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
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
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Buettner</a>
           <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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

// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_test.c
 * Applies SO3 algorithms to perform inverse and forward Wigner
 * transforms (respectively) to check that the original
 * signal is reconstructed exactly (to numerical precision). Test is
 * performed on a random signal with harmonic coefficients uniformly
 * sampled from (-1,1).
 *
 * Usage: so3_test [L [N [L0 [seed]]], e.g. so3_test 64 4 32 314
 *
 * Defaults: L = 16, N = L, L0 = 0, seed = 1
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
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
void so3_test_gen_flmn_complex(complex double *flmn, int L0, int L, int N, so3_storage_t storage, so3_n_mode_t n_mode, int seed);
void so3_test_gen_flmn_real(complex double *flmn, int L0, int L, int N, so3_storage_t storage, so3_n_mode_t n_mode, int seed);

int main(int argc, char **argv)
{
    int L, N, L0;
    complex double *flmn_orig, *flmn_syn;
    complex double *f;
    double *f_real;
    int verbosity = 0;
    int seed;
    clock_t time_start, time_end;
    int i, storage_mode, n_mode, real;
    int flmn_size;

    const char *storage_mode_str[SO3_STORE_SIZE];
    const char *n_mode_str[SO3_N_MODE_SIZE];
    const char *reality_str[2];

    storage_mode_str[SO3_STORE_ZERO_FIRST_PAD] = "padded storage with n = 0 first";
    storage_mode_str[SO3_STORE_NEG_FIRST_PAD] = "padded storage with n = -N+1 first";
    storage_mode_str[SO3_STORE_ZERO_FIRST_COMPACT] = "compact storage with n = 0 first";
    storage_mode_str[SO3_STORE_NEG_FIRST_COMPACT] = "compact storage with n = -N+1 first";

    n_mode_str[SO3_N_MODE_ALL] = "all n";
    n_mode_str[SO3_N_MODE_EVEN] = "only even n";
    n_mode_str[SO3_N_MODE_ODD] = "only odd n";
    n_mode_str[SO3_N_MODE_MAXIMUM] = "only |n| = N-1";

    reality_str[0] = "complex";
    reality_str[1] = "real";

    // one element for each combination of storage mode, n-mode and
    // real or complex signal.
    double errors[SO3_STORE_SIZE][SO3_N_MODE_SIZE][2];
    double durations_forward[SO3_STORE_SIZE][SO3_N_MODE_SIZE][2];
    double durations_inverse[SO3_STORE_SIZE][SO3_N_MODE_SIZE][2];

    // Parse command line arguments
    L = N = 16;
    L0 = 0;
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
        L0 = atoi(argv[3]);

    if (argc > 4)
        seed = atoi(argv[4]);

    // (2*N-1)*L*L is the largest number of flmn ever needed. For more
    // compact storage modes, only part of the memory will be used.
    flmn_orig = malloc((2*N-1)*L*L * sizeof *flmn_orig);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_orig);
    flmn_syn = malloc((2*N-1)*L*L * sizeof *flmn_syn);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_syn);

    f = malloc((2*L-1)*L*(2*N-1) * sizeof *f);
    SO3_ERROR_MEM_ALLOC_CHECK(f);
    f_real = malloc((2*L-1)*L*(2*N-1) * sizeof *f_real);
    SO3_ERROR_MEM_ALLOC_CHECK(f_real);

    // Write program name.
    printf("\n");
    printf("SO3 test program (C implementation)\n");
    printf("================================================================\n");

    for (storage_mode = 0; storage_mode < SO3_STORE_SIZE; ++storage_mode)
    {
        for (n_mode = 0; n_mode < SO3_N_MODE_SIZE; ++ n_mode)
        {
            // real == 0 --> complex signal
            // real == 1 --> real signal
            for (real = 0; real < 2; ++real)
            {
                for (i = 0; i <NREPEAT; ++i)
                {
                    double duration;

                    printf("Testing %s for a %s signal. Using %s. (run %d)\n",
                           storage_mode_str[storage_mode],
                           reality_str[real],
                           n_mode_str[n_mode],
                           i+1);

                    if (real) so3_test_gen_flmn_real(flmn_orig, L0, L, N, storage_mode, n_mode, seed);
                    else      so3_test_gen_flmn_complex(flmn_orig, L0, L, N, storage_mode, n_mode, seed);

                    time_start = clock();
                    if (real) so3_core_mw_inverse_via_ssht_real(f_real, flmn_orig, L0, L, N, storage_mode, n_mode, verbosity);
                    else      so3_core_mw_inverse_via_ssht(f, flmn_orig, L0, L, N, storage_mode, n_mode, verbosity);
                    time_end = clock();

                    duration = (time_end - time_start) / (double)CLOCKS_PER_SEC;
                    if (!i || duration < durations_inverse[storage_mode][n_mode][real])
                        durations_inverse[storage_mode][n_mode][real] = duration;

                    time_start = clock();
                    if (real) so3_core_mw_forward_via_ssht_real(flmn_syn, f_real, L0, L, N, storage_mode, n_mode, verbosity);
                    else      so3_core_mw_forward_via_ssht(flmn_syn, f, L0, L, N, storage_mode, n_mode, verbosity);
                    time_end = clock();

                    duration = (time_end - time_start) / (double)CLOCKS_PER_SEC;
                    if (!i || duration < durations_forward[storage_mode][n_mode][real])
                        durations_forward[storage_mode][n_mode][real] = duration;

                    flmn_size = so3_sampling_flmn_size(L, N, storage_mode, real);
                    errors[storage_mode][n_mode][real] += get_max_error(flmn_orig, flmn_syn, flmn_size)/NREPEAT;
                }
            }
        }
    }

    free(flmn_orig);
    free(flmn_syn);
    free(f);
    free(f_real);

    // =========================================================================
    // Summarise results

    printf("================================================================\n");
    printf("Summary\n\n");
    printf("Runs   = %40d\n", NREPEAT);
    printf("L0     = %40d\n", L0);
    printf("L      = %40d\n", L);
    printf("N      = %40d\n\n", N);


    // real == 0 --> complex signal
    // real == 1 --> real signal
    for (real = 0; real < 2; ++real)
    {
        printf("Results for %s signals...\n", reality_str[real]);
        for (storage_mode = 0; storage_mode < SO3_STORE_SIZE; ++storage_mode)
        {
            printf("  ...using %s...\n", storage_mode_str[storage_mode]);
            for (n_mode = 0; n_mode < SO3_N_MODE_SIZE; ++ n_mode)
            {
                printf("    ...and %s...\n", n_mode_str[n_mode]);
                printf("      Minimum time for forward transform: %fs\n", durations_forward[storage_mode][n_mode][real]);
                printf("      Minimum time for inverse transform: %fs\n", durations_inverse[storage_mode][n_mode][real]);
                printf("      Average max errors for round-trip:  %e\n", errors[storage_mode][n_mode][real]);
            }
        }
        printf("\n");
    }

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
 * Generate random Wigner coefficients of a complex signal.
 *
 * \param[out] flmn Random spherical harmonic coefficients generated. Provide
 *                  enough memory for fully padded storage, i.e. (2*N-1)*L*L
 *                  elements. Unused trailing elements will be set to zero.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] N Orientational band-limit.
 * \param[in] storage Indicates how flm blocks are stored.
 * \param[in] n_mode Indicates which n are non-zero.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_flmn_complex(
    complex double *flmn,
    int L0, int L, int N,
    so3_storage_t storage,
    so3_n_mode_t n_mode,
    int seed)
{
    int i, el, m, n, n_start, n_stop, n_inc, ind;

    for (i = 0; i < (2*N-1)*L*L; ++i)
        flmn[i] = 0.0;

    switch (n_mode)
    {
    case SO3_N_MODE_ALL:
        n_start = -N+1;
        n_stop  =  N-1;
        n_inc = 1;
        break;
    case SO3_N_MODE_EVEN:
        n_start = ((N-1) % 2 == 0) ? -N+1 : -N+2;
        n_stop  = ((N-1) % 2 == 0) ?  N-1 :  N-2;
        n_inc = 2;
        break;
    case SO3_N_MODE_ODD:
        n_start = ((N-1) % 2 != 0) ? -N+1 : -N+2;
        n_stop  = ((N-1) % 2 != 0) ?  N-1 :  N-2;
        n_inc = 2;
        break;
    case SO3_N_MODE_MAXIMUM:
        n_start = -N+1;
        n_stop  =  N-1;
        n_inc = 2*N - 2;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    for (n = n_start; n <= n_stop; n += n_inc)
    {
        for (el = MAX(L0, abs(n)); el < L; ++el)
        {
            for (m = -el; m <= el; ++m)
            {
                so3_sampling_elmn2ind(&ind, el, m, n, L, N, storage);
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
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] N Orientational band-limit.
 * \param[in] storage Indicates how flm blocks are stored.
 * \param[in] n_mode Indicates which n are non-zero.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void so3_test_gen_flmn_real(
    complex double *flmn,
    int L0, int L, int N,
    so3_storage_t storage,
    so3_n_mode_t n_mode,
    int seed)
{
    int i, el, m, n, n_start, n_stop, n_inc, ind;
    double real, imag;

    for (i = 0; i < (2*N-1)*L*L; ++i)
        flmn[i] = 0.0;

    switch (n_mode)
    {
    case SO3_N_MODE_ALL:
        n_start = 0;
        n_stop  = N-1;
        n_inc = 1;
        break;
    case SO3_N_MODE_EVEN:
        n_start = 0;
        n_stop  = ((N-1) % 2 == 0) ?  N-1 :  N-2;
        n_inc = 2;
        break;
    case SO3_N_MODE_ODD:
        n_start = 1;
        n_stop  = ((N-1) % 2 != 0) ?  N-1 :  N-2;
        n_inc = 2;
        break;
    case SO3_N_MODE_MAXIMUM:
        n_start = N-1;
        n_stop  = N-1;
        n_inc = 1;
        break;
    default:
        SO3_ERROR_GENERIC("Invalid n-mode.");
    }

    for (n = n_start; n <= n_stop; n += n_inc)
    {
        if (n == 0)
        {
            // Fill fl00 with random real values
            for (el = L0; el < L; ++el)
            {
                so3_sampling_elmn2ind_real(&ind, el, 0, 0, L, N, storage);
                flmn[ind] = (2.0*ran2_dp(seed) - 1.0);
            }

            // Fill fl+-m0 with conjugated random values

            for (el = L0; el < L; ++el)
            {
                for (m = 1; m <= el; ++m)
                {
                    real = (2.0*ran2_dp(seed) - 1.0);
                    imag = (2.0*ran2_dp(seed) - 1.0);
                    so3_sampling_elmn2ind_real(&ind, el, m, 0, L, N, storage);
                    flmn[ind] = real + imag * I;
                    so3_sampling_elmn2ind_real(&ind, el, -m, 0, L, N, storage);
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
                for (m = -el; m <= el; ++m)
                {

                    so3_sampling_elmn2ind_real(&ind, el, m, n, L, N, storage);
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

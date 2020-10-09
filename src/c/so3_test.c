// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_test.c
 * Applies SO3 algorithms to perform inverse and forward Wigner
 * transforms (respectively) to check that the original
 * signal is reconstructed exactly (to numerical precision). Test is
 * performed on a random signal with harmonic coefficients uniformly
 * sampled from (-1,1), using a variety of options.
 *
 * \par Usage
 *   \code{.sh}
 *   so3_test [L [N [L0 [seed]]]]
 *   \endcode
 *   e.g.
 *   \code{.sh}
 *   so3_test 64 4 32 314
 *   \endcode
 *   Defaults: L = 16, N = L, L0 = 0, seed = 1
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>

#include "so3.h"
#include "so3_test_utils.h"

#define NREPEAT 2
#define MIN(a,b) ((a < b) ? (a) : (b))
#define MAX(a,b) ((a > b) ? (a) : (b))

double get_max_error(complex double *expected, complex double *actual, int n);

int main(int argc, char **argv)
{
    so3_parameters_t parameters = {};
    int L, N, L0;
    complex double *flmn_orig, *flmn_syn_direct, *flmn_syn_ssht;
    complex double *f_direct, *f_ssht;
    double *f_real_direct, *f_real_ssht;
    int seed;
    clock_t time_start, time_end, time_start_ssht, time_end_ssht, time_start_direct, time_end_direct;
    int i, sampling_scheme, n_order, storage_mode, n_mode, real, routine, steerable;
    int flmn_size;

    const char *n_order_str[SO3_N_ORDER_SIZE];
    const char *storage_mode_str[SO3_STORAGE_SIZE];
    const char *n_mode_str[SO3_N_MODE_SIZE];
    const char *reality_str[2];
    //const char *routine_str[2];
    const char *sampling_str[2];
    const char *steerable_str[2];

    n_order_str[SO3_N_ORDER_ZERO_FIRST] = "n = 0 first";
    n_order_str[SO3_N_ORDER_NEGATIVE_FIRST] = "n = -N+1 first";

    storage_mode_str[SO3_STORAGE_PADDED] = "PADDED storage";
    storage_mode_str[SO3_STORAGE_COMPACT] = "COMPACT storage";

    n_mode_str[SO3_N_MODE_ALL] = "ALL n";
    n_mode_str[SO3_N_MODE_EVEN] = "only EVEN n";
    n_mode_str[SO3_N_MODE_ODD] = "only ODD n";
    n_mode_str[SO3_N_MODE_MAXIMUM] = "only |n| = N-1";
    n_mode_str[SO3_N_MODE_L] = "only |n| = el";

    reality_str[0] = "COMPLEX";
    reality_str[1] = "REAL";

    //routine_str[0] = "using SSHT";
    //routine_str[1] = "not using SSHT";

    sampling_str[0] = "MW";
    sampling_str[1] = "MWSS";

    steerable_str[0] = "NOT STEERABLE";
    steerable_str[1] = "STEERABLE";

    // one element for each combination of sampling scheme, storage mode, n-mode and
    // real or complex signal.
    //double errors[2][SO3_SAMPLING_SIZE][SO3_N_MODE_SIZE][SO3_STORAGE_SIZE][SO3_N_MODE_SIZE][2][2];
    //double durations_forward[2][SO3_SAMPLING_SIZE][SO3_N_MODE_SIZE][SO3_STORAGE_SIZE][SO3_N_MODE_SIZE][2][2];
    //double durations_inverse[2][SO3_SAMPLING_SIZE][SO3_N_MODE_SIZE][SO3_STORAGE_SIZE][SO3_N_MODE_SIZE][2][2];
    double errors[2][SO3_SAMPLING_SIZE][SO3_N_MODE_SIZE][SO3_STORAGE_SIZE][SO3_N_MODE_SIZE][2];
    double durations_forward[2][SO3_SAMPLING_SIZE][SO3_N_MODE_SIZE][SO3_STORAGE_SIZE][SO3_N_MODE_SIZE][2];
    double durations_inverse[2][SO3_SAMPLING_SIZE][SO3_N_MODE_SIZE][SO3_STORAGE_SIZE][SO3_N_MODE_SIZE][2];

    //fftw_init_threads();
    //int nthreads = omp_get_max_threads();
    //printf("Using %d threads.\n", nthreads);
    //fftw_plan_with_nthreads(nthreads);

    // Parse command line arguments
    L = N = 4;
    int show_arrays = 0;
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

    parameters.L0 = L0;
    parameters.L = L;
    parameters.N = N;
    parameters.verbosity = 0;

    // (2*N-1)*L*L is the largest number of flmn ever needed. For more
    // compact storage modes, only part of the memory will be used.
    flmn_orig = malloc((2*N-1)*L*L * sizeof *flmn_orig);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_orig);
    flmn_syn_direct = malloc((2*N-1)*L*L * sizeof *flmn_syn_direct);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_syn_direct);
    flmn_syn_ssht = malloc((2*N-1)*L*L * sizeof *flmn_syn_ssht);
    SO3_ERROR_MEM_ALLOC_CHECK(flmn_syn_ssht);

    // We only need (2*L) * (L+1) * (2*N-1) samples for MW symmetric sampling.
    // For the usual MW sampling, only part of the memory will be used.
    f_direct = malloc((2*L)*(L+1)*(2*N-1) * sizeof *f_direct);
    SO3_ERROR_MEM_ALLOC_CHECK(f_direct);
    f_real_direct = malloc((2*L)*(L+1)*(2*N-1) * sizeof *f_real_direct);
    SO3_ERROR_MEM_ALLOC_CHECK(f_real_direct);
    f_ssht = malloc((2*L)*(L+1)*(2*N-1) * sizeof *f_ssht);
    SO3_ERROR_MEM_ALLOC_CHECK(f_ssht);
    f_real_ssht = malloc((2*L)*(L+1)*(2*N-1) * sizeof *f_real_ssht);
    SO3_ERROR_MEM_ALLOC_CHECK(f_real_ssht);

    // Write program name.
    printf("\n");
    printf("SO3 test program (C implementation)\n");
    printf("================================================================\n");

    // steerable == 0 --> don't use steerable
    // steerable == 1 --> use steerable
    for (n_mode = 0; n_mode < SO3_N_MODE_SIZE; ++ n_mode)
    // for (n_mode = 0; n_mode < 1; ++ n_mode)
    {
        parameters.n_mode = n_mode;

        for (real = 0; real < 2; ++real)
        {
            parameters.reality = real;

            // MWSS not yet supported by direct routinesSO3_SAMPLING_SIZE
            for (sampling_scheme = 0; sampling_scheme < SO3_SAMPLING_SIZE; ++sampling_scheme)
            // for (sampling_scheme = 0; sampling_scheme < 1; ++sampling_scheme)
            {
                parameters.sampling_scheme = sampling_scheme;

                // For real signals, the n_order does not matter, so skip the second option
                // in that case.
                for (n_order = 0; n_order < SO3_N_ORDER_SIZE - real; ++n_order)
                // for (n_order = 0; n_order < 1; ++n_order)
                {
                    parameters.n_order = n_order;

                    // for (storage_mode = 0; storage_mode < SO3_STORAGE_SIZE; ++storage_mode)
                    for (storage_mode = 0; storage_mode < 1; ++storage_mode)
                    {
                        parameters.storage = storage_mode;

                        for (steerable = 0; steerable < 1; ++steerable)
                        {
                            parameters.steerable = steerable;

                            durations_inverse[steerable][sampling_scheme][n_order][storage_mode][n_mode][real] = 0.0;
                            durations_forward[steerable][sampling_scheme][n_order][storage_mode][n_mode][real] = 0.0;
                            errors[steerable][sampling_scheme][n_order][storage_mode][n_mode][real] = 0.0;

                            printf("\n");
                            printf("Testing a %s signal with %s with %s sampling using %s with %s. N-Mode: %s. Running %d times: ",
                                    reality_str[real],
                                    steerable_str[steerable],
                                    sampling_str[sampling_scheme],
                                    storage_mode_str[storage_mode],
                                    n_order_str[n_order],
                                    n_mode_str[n_mode],
                                    NREPEAT);
                            printf("\n");


                            printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                            printf("Tot Variation:");
                            printf(" |-| ");
                            printf(" Inv Variation:");
                            printf(" |-| ");
                            printf(" Direct Error:  ");
                            printf("|-| ");
                            printf(" SSHT Error:");
                            printf("    |-| ");
                            printf(" Inv T Direct:");
                            printf("  |-| ");
                            printf(" For T Direct:");
                            printf("  |-| ");
                            printf(" Inv T SSHT:");
                            printf("    |-| ");
                            printf(" For T SSHT:\n");
                            printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");

                            for (i = 0; i <NREPEAT; ++i)
                            {
                                int j;
                                double duration, tot, tot_inverse, error_direct, error_ssht;
                                double inverse_duration_ssht, forward_duration_ssht;
                                double inverse_duration_direct, forward_duration_direct;



                                // Reset output array
                                for (j = 0; j < (2*N-1)*L*L; ++j)
                                {
                                    flmn_syn_direct[j] = 0.0;
                                    flmn_syn_ssht[j] = 0.0;
                                }

                                //count = 0.0;
                                tot = 0.0;
                                tot_inverse = 0.0;
                                error_direct = 0.0;
                                error_ssht = 0.0;
                                
                                //parameters.steerable = steerable;

                                if (real) so3_test_gen_flmn_real(flmn_orig, &parameters, seed);
                                else      so3_test_gen_flmn_complex(flmn_orig, &parameters, seed);

                                time_start = clock();

                                time_start_direct = clock();  
                                if (real) so3_core_inverse_direct_real(f_real_direct, flmn_orig, &parameters);
                                else      so3_core_inverse_direct(f_direct, flmn_orig, &parameters);
                                time_end_direct = clock();

                                inverse_duration_direct = (time_end_direct - time_start_direct)/ (double)CLOCKS_PER_SEC;

                                time_start_ssht = clock();
                                if (real) so3_core_inverse_via_ssht_real(f_real_ssht, flmn_orig, &parameters);
                                else      so3_core_inverse_via_ssht(f_ssht, flmn_orig, &parameters);
                                time_end_ssht = clock();

                                inverse_duration_ssht = (time_end_ssht - time_start_ssht)/ (double)CLOCKS_PER_SEC;

                                time_end = clock();

                                duration = (time_end - time_start) / (double)CLOCKS_PER_SEC;
                                if (!i || duration < durations_inverse[steerable][sampling_scheme][n_order][storage_mode][n_mode][real])
                                    durations_inverse[steerable][sampling_scheme][n_order][storage_mode][n_mode][real] = duration;

                                //parameters.steerable = 0;

                                time_start = clock();
                                
                                time_start_direct = clock();
                                if (real) so3_core_forward_direct_real(flmn_syn_direct, f_real_direct, &parameters);
                                else      so3_core_forward_direct(flmn_syn_direct, f_direct, &parameters);
                                time_end_direct = clock();

                                forward_duration_direct = (time_end_direct - time_start_direct)/ (double)CLOCKS_PER_SEC;

                                time_start_ssht = clock();   
                                if (real) so3_core_forward_via_ssht_real(flmn_syn_ssht, f_real_ssht, &parameters);
                                else      so3_core_forward_via_ssht(flmn_syn_ssht, f_ssht, &parameters);
                                time_end_ssht = clock();

                                forward_duration_ssht = (time_end_ssht - time_start_ssht)/ (double)CLOCKS_PER_SEC;

                                time_end = clock();

                                duration = (time_end - time_start) / (double)CLOCKS_PER_SEC;
                                if (!i || duration < durations_forward[steerable][sampling_scheme][n_order][storage_mode][n_mode][real])
                                    durations_forward[steerable][sampling_scheme][n_order][storage_mode][n_mode][real] = duration;

                                flmn_size = so3_sampling_flmn_size(&parameters);
                                errors[steerable][sampling_scheme][n_order][storage_mode][n_mode][real] += get_max_error(flmn_orig, flmn_syn_ssht, flmn_size)/NREPEAT;
                                //printf("%.10e,",*flmn_syn);
                                //printf("%.10e",*flmn_orig);
                                //printf(" ---");

                            
                                //for (j = 0; j < (2*N-1)*L*L; ++j)
                                //    tot = 0.0;
                                //    tot_inverse = 0.0;
                                //    error_direct = 0.0;
                                //    error_ssht = 0.0;
                                //    tot += cabs(flmn_syn_direct[j] - flmn_syn_ssht[j]);
                                //    if (real) tot_inverse += cabs(f_real_direct[j] - f_real_ssht[j]);
                                //    else tot_inverse += cabs(f_direct[j] - f_ssht[j]);

                                tot = get_max_error(flmn_syn_direct, flmn_syn_ssht, flmn_size);

                                if (real) tot_inverse = get_max_error(f_real_direct, f_real_ssht, flmn_size);
                                else tot_inverse = get_max_error(f_direct, f_ssht, flmn_size);

                                error_direct = get_max_error(flmn_orig, flmn_syn_direct, flmn_size);
                                error_ssht = get_max_error(flmn_orig, flmn_syn_ssht, flmn_size);

                                printf("%.7e", tot);
                                printf("  |-|  ");
                                printf("%.7e", tot_inverse);
                                printf("  |-|  ");
                                printf("%.7e", error_direct);
                                printf("  |-|  ");
                                printf("%.7e", error_ssht);
                                printf("  |-|  ");
                                printf("%.7e", inverse_duration_direct);
                                printf("  |-|  ");
                                printf("%.7e", forward_duration_direct);
                                printf("  |-|  ");
                                printf("%.7e", inverse_duration_ssht);
                                printf("  |-|  ");
                                printf("%.7e\n", forward_duration_ssht);
                                //for (j = 0; j < (2*N-1)*L*L; ++j)
                                  //  if (real) tot = flmn_syn[j];
                                    //else tot = cabs(flmn_syn[j]);

                                   // if (tot != 0.0) count = count + tot;
                                //printf("%.30e", count);

                            }

                            if (show_arrays == 1)
                            {    
                                int c1, c2, c3;
                                c1 = c2 = c3 = 0;

                                printf("\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                printf("\n");
                                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SSHT flmn real component ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                for (int k = 0; k < (2*N-1)*L*L; ++k)
                                {
                                    if (k < 10)
                                    {
                                        printf("|| %d   ||",k);  
                                    }
                                    if (k > 9 && k < 100)
                                    {
                                        printf("|| %d  ||",k);  
                                    }
                                    if (k > 99)
                                    {
                                        printf("|| %d ||",k);  
                                    }

                                    if (creal(flmn_syn_ssht[k]) >= 0.0)
                                    {
                                        if (creal(flmn_syn_ssht[k]) == 0.0) printf(" --------- ");
                                        else printf(" %.3e ", creal(flmn_syn_ssht[k]));
                                    }
                                    else
                                    {
                                        printf("%.3e ", creal(flmn_syn_ssht[k]));
                                    }
                                    ++c1;

                                    if (c1 == 7)
                                    {
                                        printf("\n");
                                        c1 = 0;
                                    }
                                }

                                printf("\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                printf("\n");
                                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DIRECT flmn real component ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                for (int h = 0; h < (2*N-1)*L*L; ++h)
                                {
                                    if (h < 10)
                                    {
                                        printf("|| %d   ||",h);  
                                    }
                                    if (h > 9 && h < 100)
                                    {
                                        printf("|| %d  ||",h);  
                                    }
                                    if (h > 99)
                                    {
                                        printf("|| %d ||",h);  
                                    }

                                    if (creal(flmn_syn_direct[h]) >= 0.0)
                                    {
                                        if (creal(flmn_syn_direct[h]) == 0.0) printf(" --------- ");
                                        else printf(" %.3e ", creal(flmn_syn_direct[h]));
                                    }
                                    else
                                    {
                                        printf("%.3e ", creal(flmn_syn_direct[h]));
                                    }
                                
                                    ++c2;

                                    if (c2 == 7)
                                    {
                                        printf("\n");
                                        c2 = 0;
                                    }
                                }

                                printf("\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                printf("\n");
                                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRUE flmn real component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                for (int r = 0; r < (2*N-1)*L*L; ++r)
                                {
                                    if (r < 10)
                                    {
                                        printf("|| %d   ||",r);  
                                    }
                                    if (r > 9 && r < 100)
                                    {
                                        printf("|| %d  ||",r);  
                                    }
                                    if (r > 99)
                                    {
                                        printf("|| %d ||",r);  
                                    }

                                    if (creal(flmn_orig[r]) >= 0.0)
                                    {
                                        if (creal(flmn_orig[r]) == 0.0) printf(" --------- ");
                                        else printf(" %.3e ", creal(flmn_orig[r]));
                                    }
                                    else
                                    {
                                        printf("%.3e ", creal(flmn_orig[r]));
                                    }
                                    ++c3;

                                    if (c3 == 7)
                                    {
                                        printf("\n");
                                        c3 = 0;
                                    }
                                }
                            }

                            if (show_arrays == 2)
                            {    
                                int c1, c2, c3;
                                c1 = c2 = c3 = 0;

                                printf("\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                printf("\n");
                                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SSHT real part of f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                for (int k = 0; k < (2*N-1)*L*L; ++k)
                                {
                                    if (k < 10)
                                    {
                                        printf("|| %d   ||",k);  
                                    }
                                    if (k > 9 && k < 100)
                                    {
                                        printf("|| %d  ||",k);  
                                    }
                                    if (k > 99)
                                    {
                                        printf("|| %d ||",k);  
                                    }

                                    if (creal(f_ssht[k]) >= 0.0)
                                    {
                                        if (creal(f_ssht[k]) == 0.0) printf(" --------- ");
                                        else printf(" %.3e ", creal(f_ssht[k]));
                                    }
                                    else
                                    {
                                        printf("%.3e ", creal(f_ssht[k]));
                                    }
                                    ++c1;

                                    if (c1 == 7)
                                    {
                                        printf("\n");
                                        c1 = 0;
                                    }
                                }

                                printf("\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                printf("\n");
                                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DIRECT real part of f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                for (int h = 0; h < (2*N-1)*L*L; ++h)
                                {
                                    if (h < 10)
                                    {
                                        printf("|| %d   ||",h);  
                                    }
                                    if (h > 9 && h < 100)
                                    {
                                        printf("|| %d  ||",h);  
                                    }
                                    if (h > 99)
                                    {
                                        printf("|| %d ||",h);  
                                    }

                                    if (creal(f_direct[h]) >= 0.0)
                                    {
                                        if (creal(f_direct[h]) == 0.0) printf(" --------- ");
                                        else printf(" %.3e ", creal(f_direct[h]));
                                    }
                                    else
                                    {
                                        printf("%.3e ", creal(f_direct[h]));
                                    }
                                
                                    ++c2;

                                    if (c2 == 7)
                                    {
                                        printf("\n");
                                        c2 = 0;
                                    }
                                }

                                printf("\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                printf("\n");
                                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Difference between real part of  f's ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                                printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                                for (int r = 0; r < (2*N-1)*L*L; ++r)
                                {
                                    if (r < 10)
                                    {
                                        printf("|| %d   ||",r);  
                                    }
                                    if (r > 9 && r < 100)
                                    {
                                        printf("|| %d  ||",r);  
                                    }
                                    if (r > 99)
                                    {
                                        printf("|| %d ||",r);  
                                    }

                                    if (creal(f_direct[r] - f_ssht[r]) >= 0.0)
                                    {
                                        if ((creal(f_direct[r]) - creal(f_ssht[r])) <= 0.0000000000001) printf(" --------- ");
                                        else printf(" %.3e ", (creal(f_direct[r]) - creal(f_ssht[r])));
                                    }
                                    else
                                    {
                                        if ((creal(f_direct[r]) - creal(f_ssht[r])) >= -0.0000000000001) printf(" --------- ");
                                        else printf("%.3e ", (creal(f_direct[r]) - creal(f_ssht[r])));
                                    }
                                    ++c3;

                                    if (c3 == 7)
                                    {
                                        printf("\n");
                                        c3 = 0;
                                    }
                                }
                            }
                             
                                    
                            
                            printf("\n");
                            printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
                        }
                    }
                }
            }
        }
    }


    free(flmn_orig);
    free(flmn_syn_ssht);
    free(flmn_syn_direct);
    free(f_direct);
    free(f_ssht);
    free(f_real_direct);
    free(f_real_ssht);

    // =========================================================================
    // Summarise results

    printf("================================================================\n");
    printf("Summary\n\n");
    printf("Runs   = %40d\n", NREPEAT);
    printf("L0     = %40d\n", L0);
    printf("L      = %40d\n", L);
    printf("N      = %40d\n\n", N);

    for (steerable = 1; steerable < 2; ++steerable)
    {

        printf("Results for %s...\n", steerable_str[steerable]);
        // real == 0 --> complex signal
        // real == 1 --> real signal

        for (real = 0; real < 2; ++real)
        {
            //printf("  ...with %s signals...\n", reality_str[real]);

            // MWSS not yet supported by direct routines
            for (sampling_scheme = 0; sampling_scheme < 1; ++sampling_scheme)
            {
                //printf("    ...using %s sampling...\n", sampling_str[sampling_scheme]);

                //for (storage_mode = 0; storage_mode < SO3_STORAGE_SIZE; ++storage_mode)
                for (storage_mode = 0; storage_mode < 1; ++storage_mode)
                {
                    //printf("      ...with %s...\n", storage_mode_str[storage_mode]);

                    // For real signals, the n_order does not matter, so skip the second option
                    // in that case.
                    for (n_order = 0; n_order < SO3_N_ORDER_SIZE - real; ++n_order)
                    {
                        //printf("        ...using %s...\n", n_order_str[n_order]);

                        for (n_mode = 0; n_mode < SO3_N_MODE_SIZE; ++ n_mode)
                        {
                            //printf("          ...and %s...\n", n_mode_str[n_mode]);
                            //printf("            Minimum time for forward transform: %fs\n", durations_forward[steerable][sampling_scheme][n_order][storage_mode][n_mode][real]);
                            //printf("            Minimum time for inverse transform: %fs\n", durations_inverse[steerable][sampling_scheme][n_order][storage_mode][n_mode][real]);
                            //printf("            Average max errors for round-trip:  %e\n", errors[steerable][sampling_scheme][n_order][storage_mode][n_mode][real]);
                        }
                    }
                    printf("\n");
                }
            }
        }
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

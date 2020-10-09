#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "so3_types.h"
#include "so3_sampling.h"
#include "so3_types.h"
#include "so3_sampling.h"
#include "so3_conv.h"
#include "so3_error.h"
#include "so3_test_utils.h"

#define TINY 1E-9

static void test_so3_harmonic_convolution();
static void test_sampling_elmn2ind();
static void test_sampling_ind2elmn();
static void test_sampling_elmn2ind_real();
static void test_sampling_ind2elmn_real();
static void test_so3_sampling_n_loop_values();
static void test_so3_conv_get_parameters_of_convolved_lmn();

int main() {
    test_so3_sampling_n_loop_values();
    test_so3_conv_get_parameters_of_convolved_lmn();
    test_so3_harmonic_convolution();
    test_sampling_elmn2ind();
    test_sampling_ind2elmn();
    test_sampling_elmn2ind_real();
    test_sampling_ind2elmn_real();
    printf("All unit tests passed!\n");
    return 0;
}

void test_so3_sampling_n_loop_values()
{
    so3_parameters_t parameters = {};
    int n_start_old, n_stop_old, n_inc_old;
    int n_start, n_stop, n_inc;
    int N=4;
    parameters.N = N;
    for (int n_mode = 0; n_mode < SO3_N_MODE_SIZE; ++ n_mode)
    // for (n_mode = 0; n_mode < 1; ++ n_mode)
    {
        parameters.n_mode = n_mode;
        for (int real = 0; real < 2; ++real)
        {
            parameters.reality = real;

            // MWSS not yet supported by direct routinesSO3_SAMPLING_SIZE
            for (int sampling_scheme = 0; sampling_scheme < SO3_SAMPLING_SIZE; ++sampling_scheme)
            // for (sampling_scheme = 0; sampling_scheme < 1; ++sampling_scheme)
            {
                parameters.sampling_scheme = sampling_scheme;

                // For real signals, the n_order does not matter, so skip the second option
                // in that case.
                for (int n_order = 0; n_order < SO3_N_ORDER_SIZE - real; ++n_order)
                // for (n_order = 0; n_order < 1; ++n_order)
                {
                    parameters.n_order = n_order;

                    // for (storage_mode = 0; storage_mode < SO3_STORAGE_SIZE; ++storage_mode)
                    for (int storage_mode = 0; storage_mode < 1; ++storage_mode)
                    {
                        parameters.storage = storage_mode;
                        // printf("storage_mode: %d\n", storage_mode);

                        for (int steerable = 0; steerable < 1; ++steerable)
                        {
                            if (parameters.reality)
                            {
                                switch (parameters.n_mode)
                                {
                                case SO3_N_MODE_ALL:
                                case SO3_N_MODE_L:
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
                            }
                            else
                            {                            
                                switch (parameters.n_mode)
                                {
                                    case SO3_N_MODE_ALL:
                                    case SO3_N_MODE_L:
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
                                        if (N > 1)
                                            n_inc = 2*N - 2;
                                        else
                                            n_inc = 1;
                                        break;
                                    default:
                                        SO3_ERROR_GENERIC("Invalid n-mode.");
                                }
                            }
                            so3_sampling_n_loop_values(&n_start_old, &n_stop_old, &n_inc_old, &parameters);
                            assert (n_start == n_start_old);
                            assert (n_stop == n_stop_old);
                            assert (n_inc == n_inc_old);
                        }
                    }
                }
            }
        }
    }
}

void test_so3_conv_get_parameters_of_convolved_lmn()
{
    so3_parameters_t h_parameters = {}, f_parameters = {}, g_parameters = {};
    int L1=3, L01=1, N1=2;
    int L2=4, L02=0, N2=4;
    int L3=3, L03=1, N3=3;
    int seed=1;

    f_parameters.L0 = L01;
    f_parameters.L = L1;
    f_parameters.N = N1;
    f_parameters.verbosity = 1;
    f_parameters.n_mode = 0;
    f_parameters.reality = 0;
    f_parameters.sampling_scheme = 0;
    f_parameters.n_order = 0;
    f_parameters.storage = 0;
    f_parameters.steerable = 0;

    g_parameters.L = L2;
    g_parameters.L0 = L02;
    g_parameters.N = N2;
    g_parameters.verbosity = 1;
    g_parameters.n_mode = 0;
    g_parameters.reality = 0;
    g_parameters.sampling_scheme = 0;
    g_parameters.n_order = 0;
    g_parameters.storage = 0;
    g_parameters.steerable = 0;


    h_parameters = so3_conv_get_parameters_of_convolved_lmn(&f_parameters, &g_parameters);

    assert (h_parameters.L == L3);
    assert (h_parameters.L0 == L03);
    assert (h_parameters.N == N3);
}

void test_so3_harmonic_convolution()
{
    SO3_COMPLEX(double) *hlmn, *flmn, *glmn, h_value;
    so3_parameters_t h_parameters = {}, f_parameters = {}, g_parameters = {};
    int L=2, N=2;
    int el, m, n;
    int seed=1;
    int n_start, n_stop, n_inc;
    int ind_f, ind_g;

    f_parameters.L0 = 0;
    f_parameters.L = L;
    f_parameters.N = N;
    f_parameters.verbosity = 1;
    f_parameters.n_mode = 0;
    f_parameters.reality = 0;
    f_parameters.sampling_scheme = 0;
    f_parameters.n_order = 0;
    f_parameters.storage = 0;
    f_parameters.steerable = 0;

    g_parameters = f_parameters;
    h_parameters = so3_conv_get_parameters_of_convolved_lmn(&f_parameters, &g_parameters);

    int flmn_length = so3_sampling_flmn_size(&f_parameters);
    flmn = malloc(flmn_length * sizeof *flmn); SO3_ERROR_MEM_ALLOC_CHECK(flmn);

    int glmn_length = so3_sampling_flmn_size(&g_parameters);
    glmn = malloc(glmn_length * sizeof *glmn); SO3_ERROR_MEM_ALLOC_CHECK(glmn);

    int hlmn_length = so3_sampling_flmn_size(&h_parameters);
    hlmn = malloc(hlmn_length * sizeof *hlmn); SO3_ERROR_MEM_ALLOC_CHECK(hlmn);

    if (f_parameters.reality) 
    {
        so3_test_gen_flmn_real(flmn, &f_parameters, seed);
        so3_test_gen_flmn_real(glmn, &f_parameters, seed);
    }
    else
    {
        so3_test_gen_flmn_complex(flmn, &f_parameters, seed);
        so3_test_gen_flmn_complex(glmn, &f_parameters, seed);
    }

    so3_conv_harmonic_convolution(hlmn, &h_parameters, flmn, &f_parameters, glmn, &g_parameters);

    for (int i=0; i<flmn_length; i++)
    {
        h_value = 0;
        if (h_parameters.reality) so3_sampling_ind2elmn_real(&el, &m, &n, i, &h_parameters);
        else so3_sampling_ind2elmn(&el, &m, &n, i, &h_parameters);
        so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, &f_parameters);
        for (int k=n_start; k<=n_stop; k+=n_inc)
        {
            if (so3_sampling_is_elmn_non_zero(el, n, k, &g_parameters))
            {

                if (f_parameters.reality) so3_sampling_elmn2ind_real(&ind_f, el, m, k, &f_parameters);
                else so3_sampling_elmn2ind(&ind_f, el, m, k, &f_parameters);
                if (g_parameters.reality) so3_sampling_elmn2ind_real(&ind_g, el, n, k, &g_parameters);
                else so3_sampling_elmn2ind(&ind_g, el, n, k, &g_parameters);

                h_value += flmn[ind_f] * conj(glmn[ind_g]);
            }
        }
        assert ( creal(h_value - hlmn[i]) < TINY);
        assert ( cimag(h_value - hlmn[i]) < TINY);
    }

}

void test_sampling_elmn2ind()
{
    so3_parameters_t parameters = {};
    int ind;

    // Test padded storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |      -1       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.

    parameters.L = 3;
    parameters.N = -1;
    parameters.n_order = SO3_N_ORDER_ZERO_FIRST;
    parameters.storage = SO3_STORAGE_PADDED;

    so3_sampling_elmn2ind(&ind, 0, 0, 0, &parameters);
    assert( ind == 0 &&
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, &parameters);
    assert( ind == 7 &&
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, &parameters);
    assert( ind == 16 &&
            "Element (2,1,-1) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, &parameters);
    assert( ind == 25 &&
            "Element (2,1,1) in wrong position for L = 3." );

    // Test compact storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |    -1     |     1     |
    // el | 0 |     1     |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.

    parameters.L = 3;
    parameters.N = -1;
    parameters.n_order = SO3_N_ORDER_ZERO_FIRST;
    parameters.storage = SO3_STORAGE_COMPACT;

    so3_sampling_elmn2ind(&ind, 0, 0, 0, &parameters);
    assert( ind == 0 &&
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, &parameters);
    assert( ind == 7 &&
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, &parameters);
    assert( ind == 15 &&
            "Element (2,1,-1) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, &parameters);
    assert( ind == 23 &&
            "Element (2,1,1) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -2, &parameters);
    assert( ind == 28 &&
            "Element (2,1,-2) in wrong position for L = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 2, &parameters);
    assert( ind == 33 &&
            "Element (2,1,2) in wrong position for L = 3." );

    // Test padded storage with n-order ..., -2, -1, 0, 1, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |      -1       |       0       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | _ |-1 | 0 | 1 | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |

    parameters.L = 3;
    parameters.N = 3;
    parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    parameters.storage = SO3_STORAGE_PADDED;

    so3_sampling_elmn2ind(&ind, 0, 0, 0, &parameters);
    assert( ind == 18 &&
            "Element (0,0,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, -2, -2, &parameters);
    assert( ind == 4 &&
            "Element (2,-2,-2) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, &parameters);
    assert( ind == 25 &&
            "Element (2,1,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, &parameters);
    assert( ind == 16 &&
            "Element (2,1,-1) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, &parameters);
    assert( ind == 34 &&
            "Element (2,1,1) in wrong position for L = N = 3." );

    // Test compact storage with n-order ..., -2, -1, 0, 1, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |    -1     |       0       |     1     |
    // el |     1     | 0 |     1     |     1     |
    // m  |-1 | 0 | 1 | 0 |-1 | 0 | 1 |-1 | 0 | 1 |

    parameters.L = 3;
    parameters.N = 3;
    parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    parameters.storage = SO3_STORAGE_COMPACT;

    so3_sampling_elmn2ind(&ind, 0, 0, 0, &parameters);
    assert( ind == 13 &&
            "Element (0,0,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, -2, -2, &parameters);
    assert( ind == 0 &&
            "Element (2,-2,-2) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 0, &parameters);
    assert( ind == 20 &&
            "Element (2,1,0) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -1, &parameters);
    assert( ind == 11 &&
            "Element (2,1,-1) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 1, &parameters);
    assert( ind == 28 &&
            "Element (2,1,1) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, -2, &parameters);
    assert( ind == 3 &&
            "Element (2,1,-2) in wrong position for L = N = 3." );

    so3_sampling_elmn2ind(&ind, 2, 1, 2, &parameters);
    assert( ind == 33 &&
            "Element (2,1,2) in wrong position for L = N = 3." );
}

void test_sampling_ind2elmn()
{
    so3_parameters_t parameters = {};
    int el, m, n;

    // Test padded storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |      -1       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.

    parameters.L = 3;
    parameters.N = -1;
    parameters.n_order = SO3_N_ORDER_ZERO_FIRST;
    parameters.storage = SO3_STORAGE_PADDED;

    so3_sampling_ind2elmn(&el, &m, &n, 0, &parameters);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 7, &parameters);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 16, &parameters);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 16 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 25, &parameters);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 25 yields wrong indices (el,m,n)." );

    // Test compact storage with n-order 0, -1, 1, -2, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |    -1     |     1     |
    // el | 0 |     1     |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.

    parameters.L = 3;
    parameters.N = -1;
    parameters.n_order = SO3_N_ORDER_ZERO_FIRST;
    parameters.storage = SO3_STORAGE_COMPACT;

    so3_sampling_ind2elmn(&el, &m, &n, 0, &parameters);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 7, &parameters);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 15, &parameters);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 15 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 23, &parameters);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 23 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 28, &parameters);
    assert( el == 2 && m == 1 && n == -2 &&
            "Index 28 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 33, &parameters);
    assert( el == 2 && m == 1 && n == 2 &&
            "Index 33 yields wrong indices (el,m,n)." );

    // Test padded storage with n-order ..., -2, -1, 0, 1, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |      -1       |       0       |       1       |
    // el | 0 |     1     | 0 |     1     | 0 |     1     |
    // m? | _ |-1 | 0 | 1 | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |

    parameters.L = 3;
    parameters.N = 3;
    parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    parameters.storage = SO3_STORAGE_PADDED;

    so3_sampling_ind2elmn(&el, &m, &n, 18, &parameters);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 18 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 4, &parameters);
    assert( el == 2 && m == -2 && n == -2 &&
            "Index 4 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 25, &parameters);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 25 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 16, &parameters);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 16 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 34, &parameters);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 34 yields wrong indices (el,m,n)." );

    // Test compact storage with n-order ..., -2, -1, 0, 1, 2, ...
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |    -1     |       0       |     1     |
    // el |     1     | 0 |     1     |     1     |
    // m  |-1 | 0 | 1 | 0 |-1 | 0 | 1 |-1 | 0 | 1 |

    parameters.L = 3;
    parameters.N = 3;
    parameters.n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    parameters.storage = SO3_STORAGE_COMPACT;

    so3_sampling_ind2elmn(&el, &m, &n, 13, &parameters);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 13 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 0, &parameters);
    assert( el == 2 && m == -2 && n == -2 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 20, &parameters);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 20 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 11, &parameters);
    assert( el == 2 && m == 1 && n == -1 &&
            "Index 11 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 28, &parameters);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 28 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 3, &parameters);
    assert( el == 2 && m == 1 && n == -2 &&
            "Index 3 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn(&el, &m, &n, 33, &parameters);
    assert( el == 2 && m == 1 && n == 2 &&
            "Index 33 yields wrong indices (el,m,n)." );
}

void test_sampling_elmn2ind_real()
{
    so3_parameters_t parameters = {};
    int ind;

    // Test padded storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |       1       |
    // el | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |

    parameters.L = 3;
    parameters.N = 3;
    parameters.storage = SO3_STORAGE_PADDED;

    so3_sampling_elmn2ind_real(&ind, 0, 0, 0, &parameters);
    assert( ind == 0 &&
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 0, &parameters);
    assert( ind == 7 &&
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 1, &parameters);
    assert( ind == 16 &&
            "Element (2,1,1) in wrong position for L = 3." );

    // Test compact storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |     1     |
    // el | 0 |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |

    parameters.L = 3;
    parameters.N = 3;
    parameters.storage = SO3_STORAGE_COMPACT;

    so3_sampling_elmn2ind_real(&ind, 0, 0, 0, &parameters);
    assert( ind == 0 &&
            "Element (0,0,0) should be first." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 0, &parameters);
    assert( ind == 7 &&
            "Element (2,1,0) in wrong position for L = 3." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 1, &parameters);
    assert( ind == 15 &&
            "Element (2,1,1) in wrong position for L = 3." );

    so3_sampling_elmn2ind_real(&ind, 2, 1, 2, &parameters);
    assert( ind == 20 &&
            "Element (2,1,2) in wrong position for L = 3." );
}

void test_sampling_ind2elmn_real()
{
    so3_parameters_t parameters = {};
    int el, m, n;

    // Test padded storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell in the m? row corresponds to one coefficient in memory.
    // It contains an underscore for invalid el-m-n combinations
    // (which will be zeroes in memory) or the actual m-index for values
    // that will really be stored.
    //
    // n  |       0       |       1       |
    // el | 0 |     1     | 0 |     1     |
    // m? | 0 |-1 | 0 | 1 | _ |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.

    parameters.L = 3;
    parameters.N = 3;
    parameters.storage = SO3_STORAGE_PADDED;

    so3_sampling_ind2elmn_real(&el, &m, &n, 0, &parameters);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 7, &parameters);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 16, &parameters);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 16 yields wrong indices (el,m,n)." );

    // Test compact storage for a real signal (n-order 0, 1, 2, ...)
    // That is, for L = N = 2, the flmn are layed out as follows.
    // Each cell of the m-row corresponds to one coefficient in memory.
    //
    // n  |       0       |     1     |
    // el | 0 |     1     |     1     |
    // m  | 0 |-1 | 0 | 1 |-1 | 0 | 1 |
    //
    // We pass in a nonsensical N to make sure implementation does
    // not depend on it.

    parameters.L = 3;
    parameters.N = 3;
    parameters.storage = SO3_STORAGE_COMPACT;

    so3_sampling_ind2elmn_real(&el, &m, &n, 0, &parameters);
    assert( el == 0 && m == 0 && n == 0 &&
            "Index 0 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 7, &parameters);
    assert( el == 2 && m == 1 && n == 0 &&
            "Index 7 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 15, &parameters);
    assert( el == 2 && m == 1 && n == 1 &&
            "Index 15 yields wrong indices (el,m,n)." );

    so3_sampling_ind2elmn_real(&el, &m, &n, 20, &parameters);
    assert( el == 2 && m == 1 && n == 2 &&
            "Index 23 yields wrong indices (el,m,n)." );
}

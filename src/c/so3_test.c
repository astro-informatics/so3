// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details

// TODO: Write test demo.

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include <so3.h>

void print_max_error(complex double *expected, complex double *actual, int n);

int main(int argc, char **argv)
{
    int L, N;
    complex double *flmn, *f, *flmnout;
    int i;

    srand(1618);

    L = N = 3;
    flmn = malloc((2*N-1)*L*L * sizeof *flmn);
    f = malloc((2*L-1)*L*(2*N-1) * sizeof *f);
    flmnout = malloc((2*N-1)*L*L * sizeof *flmnout);

    flmn[0]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[1]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[2]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[3]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[4]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[5]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[6]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[7]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[8]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[9]  = 0.0;
    flmn[10] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[11] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[12] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[13] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[14] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[15] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[16] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[17] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[18] = 0.0;
    flmn[19] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[20] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[21] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[22] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[23] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[24] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[25] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[26] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[27] = 0.0;
    flmn[28] = 0.0;
    flmn[29] = 0.0;
    flmn[30] = 0.0;
    flmn[31] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[32] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[33] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[34] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[35] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[36] = 0.0;
    flmn[37] = 0.0;
    flmn[38] = 0.0;
    flmn[39] = 0.0;
    flmn[40] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[41] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[42] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[43] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[44] = (complex double)rand()/(complex double)RAND_MAX;

    printf("Testing padded storage with n = 0 first.\n");

    so3_core_mw_inverse_via_ssht(f, flmn, L, N, SO3_STORE_ZERO_FIRST_PAD, 0);
    so3_core_mw_forward_via_ssht(flmnout, f, L, N, SO3_STORE_ZERO_FIRST_PAD, 0);

    print_max_error(flmn, flmnout, (2*N-1)*L*L);

    flmn[0]  = 0.0;
    flmn[1]  = 0.0;
    flmn[2]  = 0.0;
    flmn[3]  = 0.0;
    flmn[4]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[5]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[6]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[7]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[8]  = (complex double)rand()/(complex double)RAND_MAX;
    flmn[9]  = 0.0;
    flmn[10] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[11] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[12] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[13] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[14] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[15] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[16] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[17] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[18] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[19] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[20] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[21] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[22] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[23] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[24] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[25] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[26] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[27] = 0.0;
    flmn[28] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[29] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[30] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[31] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[32] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[33] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[34] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[35] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[36] = 0.0;
    flmn[37] = 0.0;
    flmn[38] = 0.0;
    flmn[39] = 0.0;
    flmn[40] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[41] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[42] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[43] = (complex double)rand()/(complex double)RAND_MAX;
    flmn[44] = (complex double)rand()/(complex double)RAND_MAX;

    printf("Testing padded storage with n = -N+1 first.\n");

    so3_core_mw_inverse_via_ssht(f, flmn, L, N, SO3_STORE_NEG_FIRST_PAD, 0);
    so3_core_mw_forward_via_ssht(flmnout, f, L, N, SO3_STORE_NEG_FIRST_PAD, 0);

    print_max_error(flmn, flmnout, (2*N-1)*L*L);

    flmn = realloc(flmn, (2*N-1)*(3*L*L-N*(N-1))/3 * sizeof *flmn);
    flmnout = realloc(flmnout, (2*N-1)*(3*L*L-N*(N-1))/3 * sizeof *flmnout);

    for (i = 0; i < (2*N-1)*(3*L*L-N*(N-1))/3; ++i)
        flmn[i]  = (complex double)rand()/(complex double)RAND_MAX;

    printf("Testing compact storage with n = 0 first.\n");

    so3_core_mw_inverse_via_ssht(f, flmn, L, N, SO3_STORE_ZERO_FIRST_COMPACT, 0);
    so3_core_mw_forward_via_ssht(flmnout, f, L, N, SO3_STORE_ZERO_FIRST_COMPACT, 0);

    print_max_error(flmn, flmnout, (2*N-1)*(3*L*L-N*(N-1))/3);

    printf("Testing compact storage with n = -N+1 first.\n");

    so3_core_mw_inverse_via_ssht(f, flmn, L, N, SO3_STORE_NEG_FIRST_COMPACT, 0);
    so3_core_mw_forward_via_ssht(flmnout, f, L, N, SO3_STORE_NEG_FIRST_COMPACT, 0);

    print_max_error(flmn, flmnout, (2*N-1)*(3*L*L-N*(N-1))/3);

    free(flmn);
    free(f);
    free(flmnout);

    return 0;
}

void print_max_error(complex double *expected, complex double *actual, int n)
{
    int i;
    double error, maxError = 0;

    for (i = 0; i < n; ++i)
    {
        error = cabs(expected[i] - actual[i]);
        if (error > maxError)
            maxError = error;
    }

    printf("Maximum error is %e\n", maxError);
}

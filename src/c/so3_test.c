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

    {
        int n, el, m, offset;
        for(n = -N+1; n < N; ++n)
        {
            if (n < 0)
                offset = -2*n-1;
            else
                offset = 2*n;

            for(el = 0; el < L; ++el)
            {
                for(m = -el; m <= el; ++m)
                {
                    flmn[offset*L*L + el*el + el + m]  = (el < abs(n)) ? 0.0 : (complex double)rand()/(complex double)RAND_MAX;
                }
            }
        }
    }

    printf("Testing padded storage with n = 0 first.\n");

    so3_core_mw_inverse_via_ssht(f, flmn, L, N, SO3_STORE_ZERO_FIRST_PAD, 0);
    so3_core_mw_forward_via_ssht(flmnout, f, L, N, SO3_STORE_ZERO_FIRST_PAD, 0);

    print_max_error(flmn, flmnout, (2*N-1)*L*L);

    {
        int n, el, m;
        for(n = -N+1; n < N; ++n)
        {
            for(el = 0; el < L; ++el)
            {
                for(m = -el; m <= el; ++m)
                {
                    flmn[(n + N-1)*L*L + el*el + el + m]  = (el < abs(n)) ? 0.0 : (complex double)rand()/(complex double)RAND_MAX;
                }
            }
        }
    }

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

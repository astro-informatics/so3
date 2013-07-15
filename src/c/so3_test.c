// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details

// TODO: Write test demo.

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include <so3.h>

int main(int argc, char **argv)
{
    int L, N;
    complex double *flmn, *f, *flmnout;

    L = N = 3;
    flmn = (complex double*)calloc((2*N-1)*L*L, sizeof(complex double));
    flmn[0]  = 0.585267750979777;
    flmn[1]  = 0.223811939491137;
    flmn[2]  = 0.751267059305653;
    flmn[3]  = 0.255095115459269;
    flmn[4]  = 0.505957051665142;
    flmn[5]  = 0.699076722656686;
    flmn[6]  = 0.890903252535798;
    flmn[7]  = 0.959291425205444;
    flmn[8]  = 0.547215529963803;
    flmn[9]  = 0.0;
    flmn[10] = 0.276025076998578;
    flmn[11] = 0.679702676853675;
    flmn[12] = 0.655098003973841;
    flmn[13] = 0.162611735194631;
    flmn[14] = 0.118997681558377;
    flmn[15] = 0.498364051982143;
    flmn[16] = 0.959743958516081;
    flmn[17] = 0.340385726666133;
    flmn[18] = 0.0;
    flmn[19] = 0.138624442828679;
    flmn[20] = 0.149294005559057;
    flmn[21] = 0.257508254123736;
    flmn[22] = 0.840717255983663;
    flmn[23] = 0.254282178971531;
    flmn[24] = 0.814284826068816;
    flmn[25] = 0.243524968724989;
    flmn[26] = 0.929263623187228;
    flmn[27] = 0.0;
    flmn[28] = 0.0;
    flmn[29] = 0.0;
    flmn[30] = 0.0;
    flmn[31] = 0.489764395788231;
    flmn[32] = 0.445586200710899;
    flmn[33] = 0.646313010111265;
    flmn[34] = 0.709364830858073;
    flmn[35] = 0.754686681982361;
    flmn[36] = 0.0;
    flmn[37] = 0.0;
    flmn[38] = 0.0;
    flmn[39] = 0.0;
    flmn[40] = 0.349983765984809;
    flmn[41] = 0.196595250431208;
    flmn[42] = 0.251083857976031;
    flmn[43] = 0.616044676146639;
    flmn[44] = 0.473288848902729;

    f = (complex double*)malloc((2*L-1)*L*(2*N-1) * sizeof(complex double));
    so3_core_mw_inverse_via_ssht(f, flmn, L, N, SO3_STORE_ZERO_FIRST_PAD, 2);

    flmnout = (complex double*)malloc((2*N-1)*L*L * sizeof(complex double));

    so3_core_mw_forward_via_ssht(flmnout, f, L, N, SO3_STORE_ZERO_FIRST_PAD, 2);

    {
        int i;
        double error, maxError = 0;

        for (i = 0; i < (2*N-1)*L*L; ++i)
        {
            error = cabs(flmn[i] - flmnout[i]);
            if (error > maxError)
                maxError = error;
        }

        printf("Maximum error is %e\n", maxError);
    }

    //for(i = 0; i < (2*L-1)*L*(2*N-1); ++i)
        //printf("%f + i%f\n", creal(f[i]), cimag(f[i]));

    return 0;
}

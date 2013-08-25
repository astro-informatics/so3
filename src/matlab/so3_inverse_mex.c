// SO3 package to perform Wigner transforms
// Copyright (C) 2013 Martin Buettner and Jason McEwen
// See LICENSE.txt for license details


#include <so3.h>
#include "so3_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute inverse transform.
 *
 * Usage:
 *   [f] = ...
 *     so3_inverse_mex(flmn, L, N, order, storage);
 *
 * \author Martin Buettner
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    int i, iin, iout, a, b, g;

    int flmn_m, flmn_n;
    double *flmn_real, *flmn_imag;
    complex double *flmn;
    int L, N;
    int len;
    char order[SO3_STRING_LEN], storage[SO3_STRING_LEN];

    so3_storage_t method;

    complex double *f;
    double *f_real, *f_imag;
    mwSize ndim = 3;
    mwSize dims[ndim];
    int nalpha, nbeta, ngamma;

    /* Check number of arguments. */
    if (nrhs != 5)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:nrhs",
                          "Require five inputs.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidOutput:nlhs",
                          "Require one output.");

    /* Parse harmonic coefficients flmn. */
    iin = 0;
    flmn_m = mxGetM(prhs[iin]);
    flmn_n = mxGetN(prhs[iin]);
    if (flmn_m != 1 && flmn_n != 1)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnVector",
                          "Harmonic coefficients must be contained in vector.");
    flmn = malloc(flmn_m * flmn_n * sizeof(*flmn));
    flmn_real = mxGetPr(prhs[iin]);
    for (i = 0; i < flmn_m*flmn_n; ++i)
        flmn[i] = flmn_real[i];
    if (mxIsComplex(prhs[iin]))
    {
        flmn_imag = mxGetPi(prhs[iin]);
        for (i = 0; i < flmn_m*flmn_n; ++i)
            flmn[i] += I * flmn_imag[i];
    }

    /* Parse harmonic band-limit L. */
    iin = 1;
    if (!mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1)
    {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:harmonicBandLimit",
                          "Harmonic band-limit must be integer.");
    }
    L = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)L ||
        L <= 0)
    {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:harmonicBandLimitNonInt",
                          "Harmonic band-limit must be positive integer.");
    }

    /* Parse orientational band-limit N. */
    iin = 2;
    if (!mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1)
    {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:orientationalBandLimit",
                          "Orientational band-limit must be integer.");
    }
    N = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)N ||
        N <= 0)
    {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:orientationalBandLimitNonInt",
                          "Orientational band-limit must be positive integer.");
    }

    /* Parse storage order. */
    iin = 3;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:orderChar",
                          "Storage order must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:orderTooLong",
                          "Storage order exceeds string length.");
    mxGetString(prhs[iin], order, len);

    /* Parse storage type. */
    iin = 4;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storageChar",
                          "Storage type must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storageTooLong",
                          "Storage type exceeds string length.");
    mxGetString(prhs[iin], storage, len);

    if (strcmp(storage, SO3_STORAGE_PADDED) == 0)
    {
        if (flmn_m * flmn_n != (2*N-1)*L*L)
            mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnSize",
                              "Invalid number of harmonic coefficients.");

        if (strcmp(order, SO3_ORDER_ZEROFIRST) == 0)
            method = SO3_STORE_ZERO_FIRST_PAD;
        else if (strcmp(order, SO3_ORDER_NEGFIRST) == 0)
            method = SO3_STORE_NEG_FIRST_PAD;
        else
            mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:order",
                              "Invalid storage order.");
    }
    else if (strcmp(storage, SO3_STORAGE_COMPACT) == 0)
    {
        if (flmn_m * flmn_n != (2*N-1)*(3*L*L-N*(N-1))/3)
            mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnSize",
                              "Invalid number of harmonic coefficients.");

        if (strcmp(order, SO3_ORDER_ZEROFIRST) == 0)
            method = SO3_STORE_ZERO_FIRST_COMPACT;
        else if (strcmp(order, SO3_ORDER_NEGFIRST) == 0)
            method = SO3_STORE_NEG_FIRST_COMPACT;
        else
            mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:order",
                              "Invalid storage order.");
    }
    else
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storage",
                          "Invalid storage type.");

    /* Compute inverse transform. */
    nalpha = so3_sampling_mw_nalpha(L);
    nbeta = so3_sampling_mw_nbeta(L);
    ngamma = so3_sampling_mw_ngamma(N);

    f = malloc(nalpha * nbeta * ngamma * sizeof(*f));
    so3_core_mw_inverse_via_ssht(
        f, flmn,
        0, L, N,
        method,
        SO3_N_MODE_ALL,
        0
    );

    // Copy result to output argument

    dims[0] = ngamma;
    dims[1] = nbeta;
    dims[2] = nalpha;

    iout = 0;
    plhs[iout] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    f_real = mxGetPr(plhs[iout]);
    f_imag = mxGetPi(plhs[iout]);
    for(g = 0; g < ngamma; ++g)
    {
        for(b = 0; b < nbeta; ++b)
        {
            for(a = 0; a < nalpha; ++a)
            {
                f_real[a*ngamma*nbeta + b*ngamma + g] = creal(f[g*nalpha*nbeta + b*nalpha + a]);
                f_imag[a*ngamma*nbeta + b*ngamma + g] = cimag(f[g*nalpha*nbeta + b*nalpha + a]);
            }
        }
    }

    /* Free memory. */
    free(flmn);
    free(f);
}

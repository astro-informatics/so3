// SO3 package to perform Wigner transforms
// Copyright (C) 2013 Martin Buettner and Jason McEwen
// See LICENSE.txt for license details


#include <so3.h>
#include "so3_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute forward transform.
 *
 * Usage:
 *   [flmn] = ...
 *     so3_forward_mex(f, L, N, order, storage);
 *
 * \author Martin Buettner and Jason McEwen
 */
 void mexFunction( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
 {
    int i, iin, iout, a, b, g;

    mwSize ndim;
    const mwSize *dims;
    int f_na, f_nb, f_ng;
    double *f_real, *f_imag;
    complex double *f;
    int L, N;
    int len;
    char order[SO3_STRING_LEN], storage[SO3_STRING_LEN];
    int nalpha, nbeta, ngamma;

    so3_storage_t method;

    int flmn_size;
    complex double *flmn;
    double *flmn_real, *flmn_imag;

    /* Check number of arguments. */
    if (nrhs != 5)
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:nrhs",
                          "Require five inputs.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidOutput:nlhs",
                          "Require one output.");

    /* Parse function samples f. */
    iin = 0;
    if (mxGetNumberOfDimensions(prhs[iin]) != 3)
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:fVector",
                          "Function samples must be contained in a 3d-array.");
    dims = mxGetDimensions(prhs[iin]);
    f_na = dims[2];
    f_nb = dims[1];
    f_ng = dims[0];

    f = malloc(f_ng * f_nb * f_na * sizeof(*f));
    f_real = mxGetPr(prhs[iin]);
    for(g = 0; g < f_ng; ++g)
    {
        for(b = 0; b < f_nb; ++b)
        {
            for(a = 0; a < f_na; ++a)
            {
                f[g*f_na*f_nb + b*f_na + a] = f_real[a*f_ng*f_nb + b*f_ng + g];
            }
        }
    }
    if (mxIsComplex(prhs[iin]))
    {
        f_imag = mxGetPi(prhs[iin]);
        for(g = 0; g < f_ng; ++g)
        {
            for(b = 0; b < f_nb; ++b)
            {
                for(a = 0; a < f_na; ++a)
                {
                   f[g*f_na*f_nb + b*f_na + a] += I * f_imag[a*f_ng*f_nb + b*f_ng + g];
                }
            }
        }
    }

    /* Parse harmonic band-limit L. */
    iin = 1;
    if (!mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1)
    {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:harmonicBandLimit",
                          "Harmonic band-limit must be integer.");
    }
    L = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)L ||
        L <= 0)
    {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:harmonicBandLimitNonInt",
                          "Harmonic band-limit must be positive integer.");
    }

    /* Parse orientational band-limit N. */
    iin = 2;
    if (!mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1)
    {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:orientationalBandLimit",
                          "Orientational band-limit must be integer.");
    }
    N = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)N ||
        N <= 0)
    {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:orientationalBandLimitNonInt",
                          "Orientational band-limit must be positive integer.");
    }

    /* Parse storage order. */
    iin = 3;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:orderChar",
                          "Storage order must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:orderTooLong",
                          "Storage order exceeds string length.");
    mxGetString(prhs[iin], order, len);

    /* Parse storage type. */
    iin = 4;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:storageChar",
                          "Storage type must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:storageTooLong",
                          "Storage type exceeds string length.");
    mxGetString(prhs[iin], storage, len);

    nalpha = so3_sampling_mw_nalpha(L);
    nbeta = so3_sampling_mw_nbeta(L);
    ngamma = so3_sampling_mw_ngamma(N);

    if (f_na != nalpha || f_nb != nbeta || f_ng != ngamma)
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:fSize",
                          "Invalid dimension sizes of function samples.");

    if (strcmp(storage, SO3_STORAGE_PADDED) == 0)
    {
        flmn_size = (2*N-1)*L*L;

        if (strcmp(order, SO3_ORDER_ZEROFIRST) == 0)
            method = SO3_STORE_ZERO_FIRST_PAD;
        else if (strcmp(order, SO3_ORDER_NEGFIRST) == 0)
            method = SO3_STORE_NEG_FIRST_PAD;
        else
            mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:order",
                              "Invalid storage order.");
    }
    else if (strcmp(storage, SO3_STORAGE_COMPACT) == 0)
    {
        flmn_size = (2*N-1)*(3*L*L-N*(N-1))/3;

        if (strcmp(order, SO3_ORDER_ZEROFIRST) == 0)
            method = SO3_STORE_ZERO_FIRST_COMPACT;
        else if (strcmp(order, SO3_ORDER_NEGFIRST) == 0)
            method = SO3_STORE_NEG_FIRST_COMPACT;
        else
            mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:order",
                              "Invalid storage order.");
    }
    else
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:storage",
                          "Invalid storage type.");

    /* Compute inverse transform. */

    flmn = malloc(flmn_size * sizeof(*flmn));
    so3_core_mw_forward_via_ssht(
        flmn, f,
        L, N,
        method, 0
    );

    // Copy result to output argument

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(flmn_size, 1, mxCOMPLEX);
    flmn_real = mxGetPr(plhs[iout]);
    flmn_imag = mxGetPi(plhs[iout]);
    for(i = 0; i < flmn_size; ++i)
    {
        flmn_real[i] = creal(flmn[i]);
        flmn_imag[i] = cimag(flmn[i]);
    }

    /* Free memory. */
    free(f);
    free(flmn);
}

// SO3 package to perform Wigner transforms
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#include "ssht.h"
#include <so3.h>
#include "so3_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute inverse transform.
 *
 * Usage:
 *   [f] = ...
 *     so3_inverse_mex(flmn, L0, L, N, order, storage, n_mode, dl_method, reality, sampling_scheme);
 *
 * \author Martin Büttner
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
    int i, iin, iout, a, b, g;

    int flmn_m, flmn_n, flmn_size;
    double *flmn_real, *flmn_imag;
    complex double *flmn;
    int L0, L, N;
    int len;
    char order_str[SO3_STRING_LEN], storage_str[SO3_STRING_LEN], n_mode_str[SO3_STRING_LEN], dl_method_str[SO3_STRING_LEN], sampling_str[SO3_STRING_LEN];

    so3_n_order_t n_order;
    so3_storage_t storage_method;
    so3_n_mode_t n_mode;
    ssht_dl_method_t dl_method;
    so3_sampling_t sampling_scheme;

    so3_parameters_t parameters = {};

    int reality;

    complex double *f;
    double *f_real, *f_imag;
    double *fr;
    mwSize ndim = 3;
    mwSize dims[ndim];
    int nalpha, nbeta, ngamma;

    /* Check number of arguments. */
    if (nrhs != 10)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:nrhs",
                          "Require ten inputs.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidOutput:nlhs",
                          "Require one output.");

    /* Parse reality. */
    iin = 8;
    if( !mxIsLogicalScalar(prhs[iin]) )
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:reality",
                          "Reality flag must be logical.");
    reality = mxIsLogicalScalarTrue(prhs[iin]);

    /* Parse harmonic coefficients flmn. */
    iin = 0;
    flmn_m = mxGetM(prhs[iin]);
    flmn_n = mxGetN(prhs[iin]);
    if (flmn_m != 1 && flmn_n != 1)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnVector",
                          "Harmonic coefficients must be contained in vector.");
    flmn_size = flmn_m * flmn_n;
    flmn = malloc(flmn_size * sizeof(*flmn));
    flmn_real = mxGetPr(prhs[iin]);
    for (i = 0; i < flmn_m*flmn_n; ++i)
        flmn[i] = flmn_real[i];
    if (mxIsComplex(prhs[iin]))
    {
        flmn_imag = mxGetPi(prhs[iin]);
        for (i = 0; i < flmn_m*flmn_n; ++i)
            flmn[i] += I * flmn_imag[i];
    }

    /* Parse lower harmonic band-limit L0. */
    iin = 1;
    if (!mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1)
    {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:lowerHarmonicBandLimit",
                          "Lower harmonic band-limit must be integer.");
    }
    L0 = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)L0 ||
        L0 < 0)
    {
        mexErrMsgIdAndTxt("so3_forward_mex:InvalidInput:lowerHarmonicBandLimitNonInt",
                          "Lower harmonic band-limit must be non-negative integer.");
    }

    /* Parse harmonic band-limit L. */
    iin = 2;
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
    iin = 3;
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
    iin = 4;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:orderChar",
                          "Storage order must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:orderTooLong",
                          "Storage order exceeds string length.");
    mxGetString(prhs[iin], order_str, len);

    /* Parse storage type. */
    iin = 5;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storageChar",
                          "Storage type must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storageTooLong",
                          "Storage type exceeds string length.");
    mxGetString(prhs[iin], storage_str, len);

    if (strcmp(storage_str, SO3_STORAGE_PADDED_STR) == 0)
    {
        storage_method = SO3_STORAGE_PADDED;

        if ((reality && flmn_size != N*L*L)
            || (!reality && flmn_size != (2*N-1)*L*L))
            mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnSize",
                              "Invalid number of harmonic coefficients.");
    }
    else if (strcmp(storage_str, SO3_STORAGE_COMPACT_STR) == 0)
    {
        storage_method = SO3_STORAGE_COMPACT;

        if ((reality && flmn_size != N*(6*L*L-(N-1)*(2*N-1))/6)
            || (!reality && flmn_size != (2*N-1)*(3*L*L-N*(N-1))/3))
            mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnSize",
                              "Invalid number of harmonic coefficients.");
    }
    else
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storage",
                          "Invalid storage type.");

    if (strcmp(order_str, SO3_N_ORDER_ZERO_FIRST_STR) == 0)
        n_order = SO3_N_ORDER_ZERO_FIRST;
    else if (strcmp(order_str, SO3_N_ORDER_NEGATIVE_FIRST_STR) == 0)
        n_order = SO3_N_ORDER_NEGATIVE_FIRST;
    else
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:order",
                          "Invalid storage order.");

    /* Parse n-mode. */
    iin = 6;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:nModeChar",
                          "Storage type must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:nModeTooLong",
                          "n-mode exceeds string length.");
    mxGetString(prhs[iin], n_mode_str, len);

    if (strcmp(n_mode_str, SO3_N_MODE_ALL_STR) == 0)
        n_mode = SO3_N_MODE_ALL;
    else if (strcmp(n_mode_str, SO3_N_MODE_L_STR) == 0)
        n_mode = SO3_N_MODE_L;
    else if (strcmp(n_mode_str, SO3_N_MODE_EVEN_STR) == 0)
        n_mode = SO3_N_MODE_EVEN;
    else if (strcmp(n_mode_str, SO3_N_MODE_ODD_STR) == 0)
        n_mode = SO3_N_MODE_ODD;
    else if (strcmp(n_mode_str, SO3_N_MODE_MAXIMUM_STR) == 0)
        n_mode = SO3_N_MODE_MAXIMUM;
    else
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:nMode",
                          "Invalid n-mode.");

    /* Parse Wigner recursion method. */
    iin = 7;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:dlMethodChar",
                          "Wigner recursion method must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:dlMethodTooLong",
                          "Wigner recursion method exceeds string length.");
    mxGetString(prhs[iin], dl_method_str, len);

    if (strcmp(dl_method_str, SSHT_RECURSION_RISBO) == 0)
        dl_method = SSHT_DL_RISBO;
    else if (strcmp(dl_method_str, SSHT_RECURSION_TRAPANI) == 0)
        dl_method = SSHT_DL_TRAPANI;
    else
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:storage",
                          "Invalid Wigner recursion method.");

    /* Parse sampling scheme method. */
    iin = 9;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:samplingSchemeChar",
                          "Sampling scheme must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:samplingSchemeTooLong",
                          "Sampling scheme exceeds string length.");
    mxGetString(prhs[iin], sampling_str, len);

    parameters.L0 = L0;
    parameters.L = L;
    parameters.N = N;
    parameters.n_order = n_order;
    parameters.storage = storage_method;
    parameters.n_mode = n_mode;
    parameters.dl_method = dl_method;
    parameters.verbosity = 0;
    parameters.reality = reality;

    if (flmn_size != so3_sampling_flmn_size(&parameters))
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:flmnSize",
                          "Invalid number of harmonic coefficients.");

    if (strcmp(sampling_str, SO3_SAMPLING_MW_STR) == 0)
    {
        sampling_scheme = SO3_SAMPLING_MW;
        parameters.sampling_scheme = sampling_scheme;

        nalpha = so3_sampling_nalpha(&parameters);
        nbeta = so3_sampling_nbeta(&parameters);
        ngamma = so3_sampling_ngamma(&parameters);
    }
    else if (strcmp(sampling_str, SO3_SAMPLING_MW_SS_STR) == 0)
    {
        sampling_scheme = SO3_SAMPLING_MW_SS;
        parameters.sampling_scheme = sampling_scheme;

        nalpha = so3_sampling_nalpha(&parameters);
        nbeta = so3_sampling_nbeta(&parameters);
        ngamma = so3_sampling_ngamma(&parameters);
    }
    else
        mexErrMsgIdAndTxt("so3_inverse_mex:InvalidInput:samplingScheme",
                          "Invalid sampling scheme.");

    /* Compute inverse transform. */

    if (reality)
    {
        fr = malloc(nalpha * nbeta * ngamma * sizeof(*fr));
        so3_core_inverse_via_ssht_real(
            fr, flmn,
            &parameters
        );
    }
    else
    {
        f = malloc(nalpha * nbeta * ngamma * sizeof(*f));
        so3_core_inverse_via_ssht(
            f, flmn,
            &parameters
        );
    }

    // Copy result to output argument

    dims[0] = ngamma;
    dims[1] = nbeta;
    dims[2] = nalpha;

    if (reality)
    {
        iout = 0;
        plhs[iout] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
        f_real = mxGetPr(plhs[iout]);
        for(g = 0; g < ngamma; ++g)
        {
            for(b = 0; b < nbeta; ++b)
            {
                for(a = 0; a < nalpha; ++a)
                {
                    f_real[a*ngamma*nbeta + b*ngamma + g] = fr[g*nalpha*nbeta + b*nalpha + a];
                }
            }
        }
    }
    else
    {
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
    }

    /* Free memory. */
    free(flmn);
    if (reality)
    {
        free(fr);
    }
    else
    {
        free(f);
    }
}

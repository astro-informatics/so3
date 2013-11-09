// SO3 package to perform Wigner transforms
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#include <string.h>

#include <so3.h>
#include "so3_mex.h"
#include "mex.h"


/**
 * Compute sampling positions.
 *
 * Usage:
 *   [n, nalpha, nbeta, ngamma, alphas, betas, gammas] = so3_sampling_mex(L, N, sampling_scheme)
 *
 * \author Martin Büttner
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int L, N, nalpha, nbeta, ngamma, n, a, b, g, iin;
    double *alphas, *betas, *gammas;
    int iout = 0;
    int len;
    char sampling_str[SO3_STRING_LEN];
    so3_parameters_t parameters = {};

    /* Check number of arguments. */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:nrhs",
                "Require three inputs.");
    }
    if(nlhs!=7) {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidOutput:nlhs",
                "Require seven outputs.");
    }

    /* Parse harmonic band-limit L. */
    iin = 0;
    if( !mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1 )
    {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:bandLimit",
                          "Harmonic band-limit must be integer.");
    }
    L = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:bandLimitNonInt",
                          "Harmonic band-limit must be positive integer.");
    }

    /* Parse orientational band-limit N. */
    iin = 1;
    if( !mxIsDouble(prhs[iin]) ||
        mxIsComplex(prhs[iin]) ||
        mxGetNumberOfElements(prhs[iin])!=1 )
    {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:bandLimit",
                          "Harmonic band-limit must be integer.");
    }
    N = (int)mxGetScalar(prhs[iin]);
    if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:bandLimitNonInt",
                          "Orientational band-limit must be positive integer.");
    }
    if (N > L)
    {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:NgreaterL",
                          "Orientational band-limit must not be greater than harmonic band-limit.");
    }

    /* Parse sampling scheme method. */
    iin = 2;
    if( !mxIsChar(prhs[iin]) ) {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:samplingSchemeChar",
                          "Sampling scheme must be string.");
    }
    len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
    if (len >= SO3_STRING_LEN)
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:samplingSchemeTooLong",
                          "Sampling scheme exceeds string length.");
    mxGetString(prhs[iin], sampling_str, len);

    parameters.L = L;
    parameters.N = N;

    if (strcmp(sampling_str, SO3_SAMPLING_MW_STR) == 0)
    {
        parameters.sampling_scheme = SO3_SAMPLING_MW;
    }
    else if (strcmp(sampling_str, SO3_SAMPLING_MW_SS_STR) == 0)
    {
        parameters.sampling_scheme = SO3_SAMPLING_MW_SS;
    }
    else
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:samplingScheme",
                          "Invalid sampling scheme.");

    n = so3_sampling_n(&parameters);
    nalpha = so3_sampling_nalpha(&parameters);
    nbeta = so3_sampling_nbeta(&parameters);
    ngamma = so3_sampling_ngamma(&parameters);

    /* Compute sample positions. */
    plhs[iout++] = mxCreateDoubleScalar(n);
    plhs[iout++] = mxCreateDoubleScalar(nalpha);
    plhs[iout++] = mxCreateDoubleScalar(nbeta);
    plhs[iout++] = mxCreateDoubleScalar(ngamma);

    plhs[iout] = mxCreateDoubleMatrix(nalpha, 1, mxREAL);
    alphas = mxGetPr(plhs[iout++]);
    for (a = 0; a < nalpha; a++)
        alphas[a] = so3_sampling_a2alpha(a, &parameters);
    plhs[iout] = mxCreateDoubleMatrix(nbeta, 1, mxREAL);
    betas = mxGetPr(plhs[iout++]);
    for (b = 0; b < nbeta; b++)
        betas[b] = so3_sampling_b2beta(b, &parameters);
    plhs[iout] = mxCreateDoubleMatrix(ngamma, 1, mxREAL);
    gammas = mxGetPr(plhs[iout++]);
    for (g = 0; g < ngamma; g++)
        gammas[g] = so3_sampling_g2gamma(g, &parameters);
}

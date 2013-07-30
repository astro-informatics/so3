// SO3 package to perform Wigner transforms
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details


#include <so3.h>
#include "so3_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute sampling positions.
 *
 * Usage: 
 *   [n, nalpha, nbeta, ngamma, alphas, betas, gammas] = so3_sampling_mex(L, N)
 *
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int L, N, nalpha, nbeta, ngamma, n, a, b, g, iin;
    double *alphas, *betas, *gammas;
    int iout = 0;

    /* Check number of arguments. */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("so3_sampling_mex:InvalidInput:nrhs",
                "Require two inputs.");
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

    /* Compute sample positions. */
    nalpha = so3_sampling_mw_nalpha(L);
    nbeta = so3_sampling_mw_nbeta(L);
    ngamma = so3_sampling_mw_ngamma(N);
    n = so3_sampling_mw_n(L, N);
    plhs[iout++] = mxCreateDoubleScalar(n);
    plhs[iout++] = mxCreateDoubleScalar(nalpha);
    plhs[iout++] = mxCreateDoubleScalar(nbeta);
    plhs[iout++] = mxCreateDoubleScalar(ngamma);

    plhs[iout] = mxCreateDoubleMatrix(nalpha, 1, mxREAL);
    alphas = mxGetPr(plhs[iout++]);
    for (a = 0; a < nalpha; a++)
        alphas[a] = so3_sampling_mw_a2alpha(a, L);
    plhs[iout] = mxCreateDoubleMatrix(nbeta, 1, mxREAL);
    betas = mxGetPr(plhs[iout++]);
    for (b = 0; b < nbeta; b++)
        betas[b] = so3_sampling_mw_b2beta(b, L);
    plhs[iout] = mxCreateDoubleMatrix(ngamma, 1, mxREAL);
    gammas = mxGetPr(plhs[iout++]);
    for (g = 0; g < ngamma; g++)
        gammas[g] = so3_sampling_mw_g2gamma(g, N);

}

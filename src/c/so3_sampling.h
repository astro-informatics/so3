// SO3 package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SO3_SAMPLING
#define SO3_SAMPLING

#include <complex.h>

complex double so3_sampling_weight_mw(int p);
double so3_sampling_weight_dh(double theta_t, int L);
void so3_sampling_gl_thetas_weights(double *thetas, double *weights, int L);

double so3_sampling_mw_t2theta(int t, int L);
double so3_sampling_mw_p2phi(int p, int L);
int so3_sampling_mw_n(int L);
int so3_sampling_mw_ntheta(int L);
int so3_sampling_mw_nphi(int L);

double so3_sampling_mw_ss_t2theta(int t, int L);
double so3_sampling_mw_ss_p2phi(int p, int L);
int so3_sampling_mw_ss_n(int L);
int so3_sampling_mw_ss_ntheta(int L);
int so3_sampling_mw_ss_nphi(int L);

double so3_sampling_dh_t2theta(int t, int L);
double so3_sampling_dh_p2phi(int p, int L);
int so3_sampling_dh_n(int L);
int so3_sampling_dh_ntheta(int L);
int so3_sampling_dh_nphi(int L);

double so3_sampling_gl_p2phi(int p, int L);
int so3_sampling_gl_n(int L);
int so3_sampling_gl_ntheta(int L);
int so3_sampling_gl_nphi(int L);

extern inline void so3_sampling_elm2ind(int *ind, int el, int m);
extern inline void so3_sampling_ind2elm(int *el, int *m, int ind);

#endif

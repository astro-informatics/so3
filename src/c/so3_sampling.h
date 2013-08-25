// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Buettner and Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_SAMPLING
#define SO3_SAMPLING

#include "so3_types.h"

int so3_sampling_mw_n(int L, int N);
int so3_sampling_mw_nalpha(int L);
int so3_sampling_mw_nbeta(int L);
int so3_sampling_mw_ngamma(int N);

double so3_sampling_mw_a2alpha(int a, int L);
double so3_sampling_mw_b2beta(int b, int L);
double so3_sampling_mw_g2gamma(int g, int N);

// Note, if this is compiled using C99-standard then the "extern" belongs in the
// .c file instead.
extern inline int so3_sampling_flmn_size(int L, int N, so3_storage_t storage);
extern inline void so3_sampling_elmn2ind(int *ind, int el, int m, int n, int L, int N, so3_storage_t storage);
extern inline void so3_sampling_ind2elmn(int *el, int *m, int *n, int ind, int L, int N, so3_storage_t storage);

#endif

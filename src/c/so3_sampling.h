// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_SAMPLING
#define SO3_SAMPLING

#include "so3_types.h"

complex double so3_sampling_weight(const so3_parameters_t *parameters, int p);

int so3_sampling_f_size(const so3_parameters_t *parameters);
int so3_sampling_n(const so3_parameters_t *parameters);
int so3_sampling_nalpha(const so3_parameters_t *parameters);
int so3_sampling_nbeta(const so3_parameters_t *parameters);
int so3_sampling_ngamma(const so3_parameters_t *parameters);

double so3_sampling_a2alpha(int a, const so3_parameters_t *parameters);
double so3_sampling_b2beta(int b, const so3_parameters_t *parameters);
double so3_sampling_g2gamma(int g, const so3_parameters_t *parameters);

// Note, if this is compiled using C99-standard then the "extern" belongs in the
// .c file instead.
extern inline int so3_sampling_flmn_size(const so3_parameters_t *parameters);
extern inline void so3_sampling_elmn2ind(int *ind, int el, int m, int n, const so3_parameters_t *parameters);
extern inline void so3_sampling_ind2elmn(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters);
extern inline void so3_sampling_elmn2ind_real(int *ind, int el, int m, int n, const so3_parameters_t *parameters);
extern inline void so3_sampling_ind2elmn_real(int *el, int *m, int *n, int ind, const so3_parameters_t *parameters);

#endif

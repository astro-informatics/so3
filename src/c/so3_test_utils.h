// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_TEST_UTILS
#define SO3_TEST_UTILS

#include "so3_types.h"
#include <ssht/ssht.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif
void so3_test_gen_flmn_complex(
    complex double *flmn,
    const so3_parameters_t *parameters,
    int seed);

void so3_test_gen_flmn_real(
    complex double *flmn,
    const so3_parameters_t *parameters,
    int seed);

double ran2_dp(int idum);

#ifdef __cplusplus
}
#endif
#endif

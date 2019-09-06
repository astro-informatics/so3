// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner, Jason McEwen and Christopher Wallis
// See LICENSE.txt for license details

#ifndef SO3_CONV
#define SO3_CONV

#include "so3_types.h"
#include "ssht.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif
void so3_conv_harmonic_convolution(
    SO3_COMPLEX(double) * hlmn, const SO3_COMPLEX(double) * flmn,
    const SO3_COMPLEX(double) * glmn, const so3_parameters_t* parameters
    );
#ifdef __cplusplus
}
#endif
#endif

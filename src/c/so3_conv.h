// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner, Jason McEwen and Christopher Wallis
// See LICENSE.txt for license details

#ifndef SO3_CONV
#define SO3_CONV

#include "so3_types.h"
#include <ssht/ssht.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif
void so3_conv_convolution(
    SO3_COMPLEX(double) *h,
    const so3_parameters_t *h_parameter,
    const SO3_COMPLEX(double) *f,
    const so3_parameters_t *f_parameter,
    const SO3_COMPLEX(double) *g,
    const so3_parameters_t *g_parameter
);

void so3_conv_harmonic_convolution(
    SO3_COMPLEX(double) * hlmn, 
    const so3_parameters_t* h_parameters,
    const SO3_COMPLEX(double) * flmn,
    const so3_parameters_t* f_parameters,
    const SO3_COMPLEX(double) * glmn, 
    const so3_parameters_t* g_parameters
);

so3_parameters_t so3_conv_get_parameters_of_convolved_lmn(
    const so3_parameters_t* f_parameters,
    const so3_parameters_t* g_parameters
);

void so3_conv_get_parameters_of_convolved_lmn_void(
    so3_parameters_t* h_parameters,
    const so3_parameters_t* f_parameters,
    const so3_parameters_t* g_parameters
);

void so3_conv_s2toso3_harmonic_convolution(
    SO3_COMPLEX(double) * hlmn, 
    const so3_parameters_t* h_parameters,
    const SO3_COMPLEX(double) * flm,
    const SO3_COMPLEX(double) * glm
);


#ifdef __cplusplus
}
#endif
#endif

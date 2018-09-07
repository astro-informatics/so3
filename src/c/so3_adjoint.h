// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_ADJOINT
#define SO3_ADJOINT

#include "ssht.h"
#include <complex.h>

// void so3_adjoint_inverse_via_ssht(
//     complex double *f, const complex double *flmn,
//     const so3_parameters_t *parameters
// );

// void so3_adjoint_forward_via_ssht(
//     complex double *flmn, const complex double *f,
//     const so3_parameters_t *parameters
// );

// void so3_adjoint_inverse_via_ssht_real(
//     double *f, const complex double *flmn,
//     const so3_parameters_t *parameters
// );

// void so3_adjoint_forward_via_ssht_real(
//     complex double *flmn, const double *f,
//     const so3_parameters_t *parameters
// );



void so3_adjoint_inverse_direct(
    complex double *flmn, const complex double *f,
    const so3_parameters_t *parameters
);

void so3_adjoint_forward_direct(
    complex double *f, const complex double *flmn,
    const so3_parameters_t *parameters
);

void so3_adjoint_inverse_direct_real(
    complex double *flmn, const double *f,
    const so3_parameters_t *parameters
);

void so3_adjoint_forward_direct_real(
    double *f, const complex double *flmn,
    const so3_parameters_t *parameters
);



#endif

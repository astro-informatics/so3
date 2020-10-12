// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_CORE
#define SO3_CORE

#include "so3_types.h"
#include <ssht/ssht.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif
void so3_core_inverse_via_ssht(
    SO3_COMPLEX(double) * f, const SO3_COMPLEX(double) * flmn,
    const so3_parameters_t* parameters);

void so3_core_inverse_via_ssht_real(
    double* f, const SO3_COMPLEX(double) * flmn,
    const so3_parameters_t* parameters);

void so3_core_forward_via_ssht(
    SO3_COMPLEX(double) * flmn, const SO3_COMPLEX(double) * f,
    const so3_parameters_t* parameters);

void so3_core_forward_via_ssht_real(
    SO3_COMPLEX(double) * flmn, const double* f,
    const so3_parameters_t* parameters);

void so3_core_inverse_direct(
    SO3_COMPLEX(double) * f, const SO3_COMPLEX(double) * flmn,
    const so3_parameters_t* parameters);

void so3_core_inverse_direct_real(
    double* f, const SO3_COMPLEX(double) * flmn,
    const so3_parameters_t* parameters);

void so3_core_forward_direct(
    SO3_COMPLEX(double) * flmn, const SO3_COMPLEX(double) * f,
    const so3_parameters_t* parameters);

void so3_core_forward_direct_real(
    SO3_COMPLEX(double) * flmn, const double* f,
    const so3_parameters_t* parameters);

#ifdef __cplusplus
}
#endif
#endif

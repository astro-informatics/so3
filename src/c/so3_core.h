// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_CORE
#define SO3_CORE

#include "ssht.h"
#include <complex.h>

void so3_core_mw_inverse_via_ssht(
    complex double *f, const complex double *flmn,
    int L0, int L, int N,
    so3_sampling_t sampling,
    so3_storage_t storage,
    so3_n_mode_t n_mode,
    ssht_dl_method_t dl_method,
    int use_dft,
    int verbosity
);

void so3_core_mw_forward_via_ssht(
    complex double *flmn, const complex double *f,
    int L0, int L, int N,
    so3_sampling_t sampling,
    so3_storage_t storage,
    so3_n_mode_t n_mode,
    ssht_dl_method_t dl_method,
    int verbosity
);

void so3_core_mw_inverse_via_ssht_real(
    double *f, const complex double *flmn,
    int L0, int L, int N,
    so3_sampling_t sampling,
    so3_storage_t storage,
    so3_n_mode_t n_mode,
    ssht_dl_method_t dl_method,
    int use_dft,
    int verbosity
);

void so3_core_mw_forward_via_ssht_real(
    complex double *flmn, const double *f,
    int L0, int L, int N,
    so3_sampling_t sampling,
    so3_storage_t storage,
    so3_n_mode_t n_mode,
    ssht_dl_method_t dl_method,
    int verbosity
);


#endif

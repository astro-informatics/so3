// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details

#ifndef SO3_CORE
#define SO3_CORE

#include <complex.h>

void so3_core_mw_inverse_via_ssht(complex double *f, const complex double *flmn,
	int L, int N,
        so3_storage_t storage,
	int verbosity);

void so3_core_mw_forward_via_ssht(complex double *flmn, const complex double *f,
	int L, int N,
        so3_storage_t storage,
	int verbosity);
	

#endif

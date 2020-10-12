#ifndef SSHT_TEST_UTILITIES
#define SSHT_TEST_UTILITIES
#include <complex.h>

#include "so3/so3_types.h"

double ran2_dp(int idum);
int max(int a, int b);
void gen_flmn_real(complex double *flmn, const so3_parameters_t *parameters, int seed);
void gen_flmn_complex(
    complex double *flmn, const so3_parameters_t *parameters, int seed);
#endif

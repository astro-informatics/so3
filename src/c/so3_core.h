// SO3 package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SO3_CORE
#define SO3_CORE

#include <complex.h>
#include <ssht.h>






/*! Set of n indices for which Wigner coefficients are non-zero 
 (arises from steerable kernel functions in spherical convolution). */
typedef enum {SO3_NSET_ALL = 0, 
  SO3_NSET_ODD, 
  SO3_NSET_EVEN} so3_nset_t;

/*! Size of array used to store Wigner coefficients flmn. */
typedef enum {SO3_WIGNERSIZE_FULL = 0, 
  SO3_WIGNERSIZE_HALFM, 
  SO3_WIGNERSIZE_HALFN, 
  SO3_WIGNERSIZE_HALFMN} so3_wignersize_t;










void so3_core_mw_inverse_sov_sym(complex double *f, complex double *flm, 
                                 int L, int spin, 
                                 ssht_dl_method_t dl_method, 
                                 int verbosity);
void so3_core_mw_inverse_sov_sym_real(double *f, complex double *flm, 
                                      int L, 
                                      ssht_dl_method_t dl_method, 
                                      int verbosity);
void so3_core_mw_forward_sov_conv_sym(complex double *flm, complex double *f, 
                                      int L, int spin, 
                                      ssht_dl_method_t dl_method,
                                      int verbosity);
void so3_core_mw_forward_sov_conv_sym_real(complex double *flm, double *f, 
					    int L, 
					    ssht_dl_method_t dl_method, 
					    int verbosity);
void so3_core_mw_inverse_sov_sym_pole(complex double *f, 
				       complex double *f_sp, double *phi_sp,
				       complex double *flm, 
				       int L, int spin, 
				       ssht_dl_method_t dl_method,
				       int verbosity);
void so3_core_mw_inverse_sov_sym_real_pole(double *f, 
					    double *f_sp,
					    complex double *flm, 
					    int L, 
					    ssht_dl_method_t dl_method, 
					    int verbosity);
void so3_core_mw_forward_sov_conv_sym_pole(complex double *flm, complex double *f,
					    complex double f_sp, double phi_sp,
					    int L, int spin, 
					    ssht_dl_method_t dl_method,
					    int verbosity);
void so3_core_mw_forward_sov_conv_sym_real_pole(complex double *flm, 
						 double *f, 
						 double f_sp,
						 int L, 
						 ssht_dl_method_t dl_method, 
						 int verbosity);
// Note that mw direct algoritms are for testing purposes only.
void so3_core_mwdirect_inverse(complex double *f, complex double *flm, 
				 int L, int spin, int verbosity);
void so3_core_mwdirect_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);






#endif

// SO3 package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SO3_CORE
#define SO3_CORE

#include <complex.h>
#include <ssht.h>

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


void so3_core_mw_inverse_sov_sym_ss(complex double *f, complex double *flm, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method, 
				     int verbosity);
void so3_core_mw_inverse_sov_sym_ss_real(double *f, complex double *flm, 
					  int L, 
					  ssht_dl_method_t dl_method, 
					  int verbosity);
void so3_core_mw_forward_sov_conv_sym_ss(complex double *flm, complex double *f, 
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void so3_core_mw_forward_sov_conv_sym_ss_real(complex double *flm, double *f, 
					       int L, 
					       ssht_dl_method_t dl_method, 
					       int verbosity);
void so3_core_mw_inverse_sov_sym_ss_pole(complex double *f, 
					  complex double *f_np, double *phi_np,
					  complex double *f_sp, double *phi_sp,
					  complex double *flm, 
					  int L, int spin,
					  ssht_dl_method_t dl_method,
					  int verbosity);
void so3_core_mw_inverse_sov_sym_ss_real_pole(double *f, 
					       double *f_np,
					       double *f_sp,
					       complex double *flm, 
					       int L, 
					       ssht_dl_method_t dl_method,
					       int verbosity);
void so3_core_mw_forward_sov_conv_sym_ss_pole(complex double *flm, complex double *f,
					       complex double f_np, double phi_np,
					       complex double f_sp, double phi_sp,
					       int L, int spin, 
					       ssht_dl_method_t dl_method,
					       int verbosity);
void so3_core_mw_forward_sov_conv_sym_ss_real_pole(complex double *flm, 
						    double *f, 
						    double f_np,
						    double f_sp,
						    int L, 
						    ssht_dl_method_t dl_method,
						    int verbosity);
// Note that mw direct algoritms are for testing purposes only.
void so3_core_mwdirect_inverse_ss(complex double *f, complex double *flm, 
				   int L, int spin, int verbosity);


void so3_core_gl_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);
void so3_core_gl_inverse_sov_real(double *f, complex double *flm, 
				   int L, int verbosity);
void so3_core_gl_forward_sov(complex double *flm, complex double *f, 
			      int L, int spin, int verbosity);
void so3_core_gl_forward_sov_real(complex double *flm, double *f, 
				   int L, int verbosity);


void so3_core_dh_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);
void so3_core_dh_inverse_sov_real(double *f, complex double *flm, 
				   int L, int verbosity);
void so3_core_dh_forward_sov(complex double *flm, complex double *f, 
			      int L, int spin, int verbosity);
void so3_core_dh_forward_sov_real(complex double *flm, double *f, 
				   int L, int verbosity);


#endif

# cython: language_level=3

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#

cdef extern from "so3.h":

    int so3_sampling_f_size(so3_parameters_t *params)
    int so3_sampling_nalpha(const so3_parameters_t* parameters);
    int so3_sampling_nbeta(const so3_parameters_t* parameters);
    int so3_sampling_ngamma(const so3_parameters_t* parameters);
    int so3_sampling_flmn_size(const so3_parameters_t* parameters);
    void so3_sampling_elmn2ind(int* ind, int el, int m, int n, const so3_parameters_t* parameters);
    void so3_sampling_ind2elmn(int* el, int* m, int* n, int ind, const so3_parameters_t* parameters);
    void so3_sampling_elmn2ind_real(int* ind, int el, int m, int n, const so3_parameters_t* parameters);
    void so3_sampling_ind2elmn_real(int* el, int* m, int* n, int ind, const so3_parameters_t* parameters);
    void so3_sampling_n_loop_values(int *n_start, int *n_stop, int *n_inc, const so3_parameters_t *parameters);
    void so3_sampling_el_loop_values(int *el_start, int *el_stop, int *el_inc, const int n, const so3_parameters_t *parameters);
    void so3_sampling_m_loop_values(int *m_start, int *m_stop, int *m_inc, const int el);
    int so3_sampling_is_elmn_non_zero_return_int(const int el, const int m, const int n, const so3_parameters_t *parameters);

    void so3_core_inverse_direct(double complex * f, const double complex * flmn, const so3_parameters_t* parameters)

    void so3_core_inverse_direct_real(double* f, const double complex * flmn,const so3_parameters_t* parameters)

    void so3_core_forward_direct(double complex * flmn, const double complex * f, const so3_parameters_t* parameters)

    void so3_core_forward_direct_real(
        double complex * flmn, const double* f,
        const so3_parameters_t* parameters)


    void so3_conv_convolution(
        double complex *h,
        const so3_parameters_t *h_parameter,
        const double complex *f,
        const so3_parameters_t *f_parameter,
        const double complex *g,
        const so3_parameters_t *g_parameter
    )

    void so3_conv_harmonic_convolution(
        double complex * hlmn, 
        const so3_parameters_t* h_parameters,
        const double complex * flmn,
        const so3_parameters_t* f_parameters,
        const double complex * glmn, 
        const so3_parameters_t* g_parameters
    )

    void so3_conv_get_parameters_of_convolved_lmn_void(
        so3_parameters_t* h_parameters,
        const so3_parameters_t* f_parameters,
        const so3_parameters_t* g_parameters
    )

    void so3_conv_s2toso3_harmonic_convolution(
        double complex  * hlmn, 
        const so3_parameters_t* h_parameters,
        const double complex * flm,
        const double complex * glm, 
    )

    ctypedef struct so3_parameters_t:
        int verbosity
        int reality
        int L0
        int L
        int N
        so3_sampling_t sampling_scheme
        so3_n_order_t n_order
        so3_storage_t storage
        so3_n_mode_t n_mode
        ssht_dl_method_t dl_method
        int steerable 

    ctypedef enum so3_sampling_t:
            SO3_SAMPLING_MW, SO3_SAMPLING_MW_SS, SO3_SAMPLING_SIZE

    ctypedef enum so3_n_order_t:
        SO3_N_ORDER_ZERO_FIRST, SO3_N_ORDER_NEGATIVE_FIRST, SO3_N_ORDER_SIZE

    ctypedef enum so3_storage_t:
        SO3_STORAGE_PADDED, SO3_STORAGE_COMPACT, SO3_STORAGE_SIZE

    ctypedef enum so3_n_mode_t:
        SO3_N_MODE_ALL, SO3_N_MODE_EVEN, SO3_N_MODE_ODD, SO3_N_MODE_MAXIMUM, SO3_N_MODE_L, SO3_N_MODE_SIZE

    ctypedef enum ssht_dl_method_t:
        SSHT_DL_RISBO, SSHT_DL_TRAPANI

# funcitons to include

def create_parameter_dict(
        int L,
        int N,
        int L0=0,
        int verbosity=0,
        int reality=0,
        so3_sampling_t sampling_scheme=SO3_SAMPLING_MW,
        so3_n_order_t n_order=SO3_N_ORDER_NEGATIVE_FIRST,
        so3_storage_t storage=SO3_STORAGE_PADDED,
        so3_n_mode_t n_mode=SO3_N_MODE_ALL,
        ssht_dl_method_t dl_method=SSHT_DL_RISBO,
        int steerable=0
        ):
    """function to create params dict from input"""
    cdef so3_parameters_t parameters = {}
    parameters.L = L
    parameters.N = N
    parameters.L0 = L0
    parameters.verbosity = verbosity
    parameters.reality = reality
    parameters.sampling_scheme = sampling_scheme
    parameters.n_order = n_order
    parameters.storage = storage
    parameters.n_mode = n_mode
    parameters.dl_method = dl_method
    parameters.steerable = steerable

    return parameters

cdef so3_parameters_t create_parameter_struct(dict parameters_dict):
    cdef so3_parameters_t parameters = {}
    parameters.L = parameters_dict['L']
    parameters.N = parameters_dict['N']
    parameters.L0 = parameters_dict['L0']
    parameters.verbosity = parameters_dict['verbosity']
    parameters.reality = parameters_dict['reality']
    parameters.sampling_scheme = parameters_dict['sampling_scheme']
    parameters.n_order = parameters_dict['n_order']
    parameters.storage = parameters_dict['storage']
    parameters.n_mode = parameters_dict['n_mode']
    parameters.dl_method = parameters_dict['dl_method']
    parameters.steerable = parameters_dict['steerable']

    return parameters


# all of the sampling.h functions! (except index to angle)

def f_size(dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    return so3_sampling_f_size(&parameters)

def n_alpha(dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    return so3_sampling_nalpha(&parameters)

def n_beta(dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    return so3_sampling_nbeta(&parameters)

def n_gamma(dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    return so3_sampling_ngamma(&parameters)

def flmn_size(dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    return so3_sampling_flmn_size(&parameters)

def elmn2ind(int el, int m, int n, dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    cdef int ind
    if parameters_dict['reality'] == 1:
        so3_sampling_elmn2ind_real(&ind, el, m, n, &parameters)
    else:
        so3_sampling_elmn2ind(&ind, el, m, n, &parameters)
    return ind

def ind2elmn(int ind, dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    cdef int el, m, n
    if parameters_dict['reality'] == 1:
        so3_sampling_ind2elmn_real(&el, &m, &n, ind, &parameters)
    else:
        so3_sampling_ind2elmn(&el, &m, &n, ind, &parameters)
    return (el, m, n)

def get_loop_n_values(dict parameters_dict):
    cdef int n_start, n_stop, n_inc
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, &parameters)
    return n_start, n_stop, n_inc

def loop_over_n(dict parameters_dict):
    cdef int n_start, n_stop, n_inc, n
    n_start, n_stop, n_inc = get_loop_n_values(parameters_dict)
    for n in range(n_start, n_stop+1, n_inc):
        yield n

def get_loop_el_values(int n, dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)
    cdef int el_start, el_stop, el_inc
    so3_sampling_el_loop_values(&el_start, &el_stop, &el_inc, n, &parameters)
    return el_start, el_stop, el_inc

def loop_over_el(int n, dict parameters_dict):
    cdef int el_start, el_stop, el_inc, el
    el_start, el_stop, el_inc = get_loop_el_values(n, parameters_dict)
    for el in range(el_start, el_stop+1, el_inc):
        yield el

def get_loop_m_values(int el):
    cdef int m_start, m_stop, m_inc
    so3_sampling_m_loop_values(&m_start, &m_stop, &m_inc, el)
    return m_start, m_stop, m_inc

def loop_over_m(int el):
    cdef int m_start, m_stop, m_inc, m
    m_start, m_stop, m_inc = get_loop_m_values(el)
    for m in range(m_start, m_stop+1, m_inc):
        yield m

def is_elmn_non_zero(int el, int m, int n, parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)

    return bool(so3_sampling_is_elmn_non_zero_return_int(el, m, n, &parameters))


# forward and inverse for MW and MWSS for complex functions
def inverse(np.ndarray[ double complex, ndim=1, mode="c"] flmn not None, dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)

    if parameters_dict['reality']:
        f_length = f_size(parameters_dict)
        f = np.zeros([f_length,], dtype=float)
        so3_core_inverse_direct_real(<double *> np.PyArray_DATA(f), <const double complex*> np.PyArray_DATA(flmn), &parameters)
    else:
        f_length = f_size(parameters_dict)
        f = np.zeros([f_length,], dtype=complex)
        so3_core_inverse_direct(<double complex*> np.PyArray_DATA(f), <const double complex*> np.PyArray_DATA(flmn), &parameters)

    return f

def forward(np.ndarray[ double complex, ndim=1, mode="c"] f not None, dict parameters_dict):
    cdef so3_parameters_t parameters=create_parameter_struct(parameters_dict)

    if parameters_dict['reality']:
        flmn_length = flmn_size(parameters_dict)
        flmn = np.zeros([flmn_length,], dtype=float)
        so3_core_forward_direct_real(<double complex*> np.PyArray_DATA(flmn), <const double *> np.PyArray_DATA(f), &parameters)
    else:
        flmn_length = flmn_size(parameters_dict)
        flmn = np.zeros([flmn_length,], dtype=complex)
        so3_core_forward_direct(<double complex*> np.PyArray_DATA(flmn), <const double complex*> np.PyArray_DATA(f), &parameters)

    return flmn

# convolution both in real and harmonic space and helper params function

def get_convolved_parameter_dict(dict f_parameters_dict, dict g_parameters_dict):
    cdef so3_parameters_t f_parameters=create_parameter_struct(f_parameters_dict)
    cdef so3_parameters_t g_parameters=create_parameter_struct(f_parameters_dict)
    cdef so3_parameters_t h_parameters
    
    so3_conv_get_parameters_of_convolved_lmn_void(
        &h_parameters,
        &f_parameters,
        &g_parameters
    )    
    return h_parameters

def convolve(
    np.ndarray[ double complex, ndim=1, mode="c"] f not None, 
    dict f_parameters_dict,
    np.ndarray[ double complex, ndim=1, mode="c"] g not None, 
    dict g_parameters_dict
    ):
    cdef so3_parameters_t f_parameters=create_parameter_struct(f_parameters_dict)
    cdef so3_parameters_t g_parameters=create_parameter_struct(f_parameters_dict)
    cdef so3_parameters_t h_parameters

    so3_conv_get_parameters_of_convolved_lmn_void(
        &h_parameters,
        &f_parameters,
        &g_parameters
    )    

    h_length = f_size(<dict>h_parameters)
    h = np.zeros([h_length,], dtype=complex)

    so3_conv_convolution(
        <double complex *> np.PyArray_DATA(h),
        &h_parameters,
        <const double complex *> np.PyArray_DATA(f),
        &f_parameters,
        <const double complex *> np.PyArray_DATA(g),
        &g_parameters
    )
    return h, h_parameters

def convolve_harmonic(    
    np.ndarray[ double complex, ndim=1, mode="c"] flmn not None, 
    dict f_parameters_dict,
    np.ndarray[ double complex, ndim=1, mode="c"] glmn not None, 
    dict g_parameters_dict
    ):

    cdef so3_parameters_t f_parameters=create_parameter_struct(f_parameters_dict)
    cdef so3_parameters_t g_parameters=create_parameter_struct(f_parameters_dict)
    cdef so3_parameters_t h_parameters

    so3_conv_get_parameters_of_convolved_lmn_void(
        &h_parameters,
        &f_parameters,
        &g_parameters
    )    

    hlmn_length = flmn_size(<dict>h_parameters)
    hlmn = np.zeros([hlmn_length,], dtype=complex)

    so3_conv_convolution(
        <double complex *> np.PyArray_DATA(hlmn),
        &h_parameters,
        <const double complex *> np.PyArray_DATA(flmn),
        &f_parameters,
        <const double complex *> np.PyArray_DATA(glmn),
        &g_parameters
    )
    return hlmn, h_parameters

def s2toso3_harmonic_convolution(
    dict h_parameters_dict,
    np.ndarray[ double complex, ndim=1, mode="c"] flm not None,
    np.ndarray[ double complex, ndim=1, mode="c"] glm not None):

    cdef so3_parameters_t h_parameters=create_parameter_struct(h_parameters_dict)

    hlmn_length = flmn_size(<dict>h_parameters)
    hlmn = np.zeros([hlmn_length,], dtype=complex)

    so3_conv_s2toso3_harmonic_convolution(
    <double complex *> np.PyArray_DATA(hlmn), 
    &h_parameters,
    <const double complex *> np.PyArray_DATA(flm),
    <const double complex *> np.PyArray_DATA(glm)
    )
    return hlmn

def test_func():
    return "hello"
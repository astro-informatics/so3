# cython: language_level=3

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# if you want to use the Numpy-C-API from Cython
np.import_array()

#----------------------------------------------------------------------------------------------------#

cdef extern from "so3/so3.h":

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

    void so3_core_inverse_via_ssht(double complex * f, const double complex * flmn, const so3_parameters_t* parameters)

    void so3_core_inverse_via_ssht_real(double* f, const double complex * flmn,const so3_parameters_t* parameters)

    void so3_core_forward_via_ssht(double complex * flmn, const double complex * f, const so3_parameters_t* parameters)

    void so3_core_forward_via_ssht_real(
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

class SO3Parameters:
    def __init__(self,
        int L=32,
        int N=4,
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
        self.L = L
        self.N = N
        self.L0 = L0
        self.verbosity = verbosity
        self.reality = reality
        self.sampling_scheme = sampling_scheme
        self.n_order = n_order
        self.storage = storage
        self.n_mode = n_mode
        self.dl_method = dl_method
        self.steerable = steerable    

    def from_dict(self, dict parameters_dict):
        self.L = parameters_dict['L']
        self.N = parameters_dict['N']
        self.L0 = parameters_dict['L0']
        self.verbosity = parameters_dict['verbosity']
        self.reality = parameters_dict['reality']
        self.sampling_scheme = parameters_dict['sampling_scheme']
        self.n_order = parameters_dict['n_order']
        self.storage = parameters_dict['storage']
        self.n_mode = parameters_dict['n_mode']
        self.dl_method = parameters_dict['dl_method']
        self.steerable = parameters_dict['steerable']
        return self


def create_parameter_dict(
        int L,
        int N,
        int L0=0,
        int verbosity=0,
        int reality=0,
        str sampling_scheme_str="SO3_SAMPLING_MW",
        str n_order_str="SO3_N_ORDER_NEGATIVE_FIRST",
        str storage_str="SO3_STORAGE_PADDED",
        str n_mode_str="SO3_N_MODE_ALL",
        str dl_method_str="SSHT_DL_RISBO",
        int steerable=0
        ):
    """function to create params class from input
    this is just here to prevent some breaking changes
    """
    cdef so3_sampling_t sampling_scheme
    cdef so3_n_order_t n_order
    cdef so3_storage_t storage
    cdef so3_n_mode_t n_mode
    cdef ssht_dl_method_t dl_method
    
    if sampling_scheme_str == "SO3_SAMPLING_MW":
        sampling_scheme = SO3_SAMPLING_MW
    else:
        sampling_scheme = SO3_SAMPLING_MW_SS
    
    if n_order_str == "SO3_N_ORDER_NEGATIVE_FIRST":
        n_order = SO3_N_ORDER_NEGATIVE_FIRST
    else:
        n_order = SO3_N_ORDER_ZERO_FIRST
    
    if storage_str=="SO3_STORAGE_PADDED":
        storage = SO3_STORAGE_PADDED
    else:
        storage = SO3_STORAGE_COMPACT
    
    if n_mode_str=="SO3_N_MODE_ALL":
        n_mode = SO3_N_MODE_ALL
    elif n_mode_str=="SO3_N_MODE_EVEN":
        n_mode = SO3_N_MODE_EVEN
    elif n_mode_str=="SO3_N_MODE_ODD":
        n_mode = SO3_N_MODE_ODD
    elif n_mode_str=="SO3_N_MODE_MAXIMUM":
        n_mode = SO3_N_MODE_MAXIMUM
    else: 
        n_mode = SO3_N_MODE_L
    
    if dl_method_str=="SSHT_DL_RISBO":
        dl_method = SSHT_DL_RISBO
    else:
        dl_method = SSHT_DL_TRAPANI
    
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

    return SO3Parameters(
    L = L,
    N = N,
    L0 = L0,
    verbosity = verbosity,
    reality = reality,
    sampling_scheme = sampling_scheme,
    n_order = n_order,
    storage = storage,
    n_mode = n_mode,
    dl_method = dl_method,
    steerable = steerable
    )

cdef so3_parameters_t create_parameter_struct(so3_parameters):
    cdef so3_parameters_t parameters = {}
    parameters.L = so3_parameters.L
    parameters.N = so3_parameters.N
    parameters.L0 = so3_parameters.L0
    parameters.verbosity = so3_parameters.verbosity
    parameters.reality = so3_parameters.reality
    parameters.sampling_scheme = so3_parameters.sampling_scheme
    parameters.n_order = so3_parameters.n_order
    parameters.storage = so3_parameters.storage
    parameters.n_mode = so3_parameters.n_mode
    parameters.dl_method = so3_parameters.dl_method
    parameters.steerable = so3_parameters.steerable

    return parameters


# all of the sampling.h functions! (except index to angle)

def f_size(so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    return so3_sampling_f_size(&parameters)

def n_alpha(so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    return so3_sampling_nalpha(&parameters)

def n_beta(so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    return so3_sampling_nbeta(&parameters)

def n_gamma(so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    return so3_sampling_ngamma(&parameters)

def flmn_size(so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    return so3_sampling_flmn_size(&parameters)

def elmn2ind(int el, int m, int n, so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    cdef int ind
    if so3_parameters.reality == 1:
        so3_sampling_elmn2ind_real(&ind, el, m, n, &parameters)
    else:
        so3_sampling_elmn2ind(&ind, el, m, n, &parameters)
    return ind

def ind2elmn(int ind, so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    cdef int el, m, n
    if so3_parameters.reality == 1:
        so3_sampling_ind2elmn_real(&el, &m, &n, ind, &parameters)
    else:
        so3_sampling_ind2elmn(&el, &m, &n, ind, &parameters)
    return (el, m, n)

def get_loop_n_values(so3_parameters not None):
    cdef int n_start, n_stop, n_inc
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    so3_sampling_n_loop_values(&n_start, &n_stop, &n_inc, &parameters)
    return n_start, n_stop, n_inc

def loop_over_n(so3_parameters not None):
    cdef int n_start, n_stop, n_inc, n
    n_start, n_stop, n_inc = get_loop_n_values(so3_parameters)
    for n in range(n_start, n_stop+1, n_inc):
        yield n

def get_loop_el_values(int n, so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)
    cdef int el_start, el_stop, el_inc
    so3_sampling_el_loop_values(&el_start, &el_stop, &el_inc, n, &parameters)
    return el_start, el_stop, el_inc

def loop_over_el(int n, so3_parameters not None):
    cdef int el_start, el_stop, el_inc, el
    el_start, el_stop, el_inc = get_loop_el_values(n, so3_parameters)
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

def is_elmn_non_zero(int el, int m, int n, so3_parameters):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)

    return bool(so3_sampling_is_elmn_non_zero_return_int(el, m, n, &parameters))


# forward and inverse for MW and MWSS for complex functions
def inverse(np.ndarray[ double complex, ndim=1, mode="c"] flmn not None, so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)

    if so3_parameters.reality:
        f_length = f_size(so3_parameters)
        f = np.zeros([f_length,], dtype=float)
        so3_core_inverse_via_ssht_real(<double *> np.PyArray_DATA(f), <const double complex*> np.PyArray_DATA(flmn), &parameters)
    else:
        f_length = f_size(so3_parameters)
        f = np.zeros([f_length,], dtype=complex)
        so3_core_inverse_via_ssht(<double complex*> np.PyArray_DATA(f), <const double complex*> np.PyArray_DATA(flmn), &parameters)

    return f

def forward(np.ndarray[ double complex, ndim=1, mode="c"] f not None, so3_parameters not None):
    cdef so3_parameters_t parameters=create_parameter_struct(so3_parameters)

    if so3_parameters.reality:
        flmn_length = flmn_size(so3_parameters)
        flmn = np.zeros([flmn_length,], dtype=float)
        so3_core_forward_via_ssht_real(<double complex*> np.PyArray_DATA(flmn), <const double *> np.PyArray_DATA(f), &parameters)
    else:
        flmn_length = flmn_size(so3_parameters)
        flmn = np.zeros([flmn_length,], dtype=complex)
        so3_core_forward_via_ssht(<double complex*> np.PyArray_DATA(flmn), <const double complex*> np.PyArray_DATA(f), &parameters)

    return flmn

# convolution both in real and harmonic space and helper params function

def get_convolved_parameter_dict(dict f_so3_parameters_dict, dict g_so3_parameters_dict):
    """To avoid *some* breaking changes """
    print("This function is depreaciated use SO3Parameters class to define parameters instead")
    f_so3_parameters = SO3Parameters().from_dict(f_so3_parameters_dict)
    g_so3_parameters = SO3Parameters().from_dict(g_so3_parameters_dict)
    return get_convolved_parameters(f_so3_parameters, g_so3_parameters)

def get_convolved_parameters(f_so3_parameters, g_so3_parameters):
    cdef so3_parameters_t f_parameters=create_parameter_struct(f_so3_parameters)
    cdef so3_parameters_t g_parameters=create_parameter_struct(g_so3_parameters)
    cdef so3_parameters_t h_parameters
    
    so3_conv_get_parameters_of_convolved_lmn_void(
        &h_parameters,
        &f_parameters,
        &g_parameters
    )    
    
    return SO3Parameters().from_dict(h_parameters)

def convolve(
    np.ndarray[ double complex, ndim=1, mode="c"] f not None, 
    f_parameters,
    np.ndarray[ double complex, ndim=1, mode="c"] g not None, 
    g_parameters
    ):

    h_parameters = get_convolved_parameters(f_parameters, g_parameters)

    cdef so3_parameters_t f_parameters_struct=create_parameter_struct(f_parameters)
    cdef so3_parameters_t g_parameters_struct=create_parameter_struct(g_parameters)
    cdef so3_parameters_t h_parameters_struct=create_parameter_struct(h_parameters)

    h_length = f_size(h_parameters)
    h = np.zeros([h_length,], dtype=complex)

    so3_conv_convolution(
        <double complex *> np.PyArray_DATA(h),
        &h_parameters_struct,
        <const double complex *> np.PyArray_DATA(f),
        &f_parameters_struct,
        <const double complex *> np.PyArray_DATA(g),
        &g_parameters_struct
    )
    return h, h_parameters

def convolve_harmonic(    
    np.ndarray[ double complex, ndim=1, mode="c"] flmn not None, 
    f_parameters,
    np.ndarray[ double complex, ndim=1, mode="c"] glmn not None, 
    g_parameters
    ):


    h_parameters = get_convolved_parameters(f_parameters, g_parameters)

    cdef so3_parameters_t f_parameters_struct=create_parameter_struct(f_parameters)
    cdef so3_parameters_t g_parameters_struct=create_parameter_struct(g_parameters)
    cdef so3_parameters_t h_parameters_struct=create_parameter_struct(h_parameters)

    hlmn_length = flmn_size(h_parameters)
    hlmn = np.zeros([hlmn_length,], dtype=complex)

    so3_conv_convolution(
        <double complex *> np.PyArray_DATA(hlmn),
        &h_parameters_struct,
        <const double complex *> np.PyArray_DATA(flmn),
        &f_parameters_struct,
        <const double complex *> np.PyArray_DATA(glmn),
        &g_parameters_struct
    )
    return hlmn, h_parameters

def s2toso3_harmonic_convolution(
    h_so3_parameters,
    np.ndarray[ double complex, ndim=1, mode="c"] flm not None,
    np.ndarray[ double complex, ndim=1, mode="c"] glm not None):

    cdef so3_parameters_t h_parameters=create_parameter_struct(h_so3_parameters)

    hlmn_length = flmn_size(h_so3_parameters)
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

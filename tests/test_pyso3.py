import numpy as np
from pytest import approx

import so3


def test_pyso3():

    f_p = so3.create_parameter_dict(16, 8)
    f_length = so3.f_size(f_p)
    rng = np.random.default_rng()
    f_before = rng.normal(size=(f_length, 2)) @ [1, 1j]
    flmn = so3.forward(f_before, f_p)
    f_before = so3.inverse(flmn, f_p)

    flmn = so3.forward(f_before, f_p)

    f = so3.inverse(flmn, f_p)
    flmn_new = so3.forward(f, f_p)

    assert flmn_new == approx(flmn)
    assert f_before == approx(f)


def test_convolve_smoke_test():
    rng = np.random.default_rng()
    g_p = so3.create_parameter_dict(16, 8)
    g_length = so3.f_size(g_p)
    g_before = rng.normal(size=(g_length, 2)) @ [1, 1j]
    f_before = rng.normal(size=(g_length, 2)) @ [1, 1j]
    f = so3.inverse(so3.forward(f_before, g_p), g_p)
    g = so3.inverse(so3.forward(g_before, g_p), g_p)
    so3.convolve(f, g_p, g, g_p)

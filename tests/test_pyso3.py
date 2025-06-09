import numpy as np
import pytest
from pytest import approx

import so3


@pytest.mark.parametrize("L, N, L0", [(16, 8, 0), (16, 8, 1), (4, 1, 0), (3, 2, 0)])
@pytest.mark.parametrize("reality", [False, True])
def test_pyso3(L, N, L0, reality):
    f_p = so3.create_parameter_dict(L, N, L0, reality=reality * 1)
    f_length = so3.f_size(f_p)
    rng = np.random.default_rng()
    multiplier = [1, 1] if reality else [1, 1j]
    f_before = rng.normal(size=(f_length, 2)) @ multiplier
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

import numpy as np

from plotter.deltar_utils import delta_r


def test_delta_r_pure_eta_separation():
    np.testing.assert_allclose(delta_r(0.0, 0.0, 1.5, 0.0), 1.5)


def test_delta_r_pure_phi_separation():
    np.testing.assert_allclose(delta_r(0.0, 0.0, 0.0, 1.0), 1.0)


def test_delta_r_pythagorean():
    # phi difference must stay within [0, pi] to be physical (no wraparound).
    np.testing.assert_allclose(delta_r(0.0, 0.0, 4.0, 3.0), 5.0)


def test_delta_r_phi_wraps_around_2pi():
    # phi = -pi+0.1 and phi = pi-0.1 are 0.2 rad apart in reality, not ~2pi.
    dr = delta_r(0.0, -np.pi + 0.1, 0.0, np.pi - 0.1)
    np.testing.assert_allclose(dr, 0.2, atol=1e-9)

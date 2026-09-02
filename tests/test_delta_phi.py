import numpy as np

from plotter.delta_phi_utils import compute_delta_phi_deg


def test_delta_phi_basic():
    np.testing.assert_allclose(compute_delta_phi_deg(0.0, np.pi / 2), 90.0)
    np.testing.assert_allclose(compute_delta_phi_deg(0.0, 0.0), 0.0)


def test_delta_phi_wraps_around_2pi():
    # phi values on opposite sides of the +/-pi branch cut are physically
    # close together (small angle), not ~360 degrees apart.
    phi1 = -np.pi + 0.1
    phi2 = np.pi - 0.1
    dphi = compute_delta_phi_deg(phi1, phi2)
    np.testing.assert_allclose(dphi, np.degrees(0.2), atol=1e-9)


def test_delta_phi_never_exceeds_180():
    phi1 = np.linspace(-np.pi, np.pi, 37)
    phi2 = np.linspace(-np.pi, np.pi, 37)[::-1]
    dphi = compute_delta_phi_deg(phi1, phi2)
    assert np.all(dphi <= 180.0 + 1e-9)
    assert np.all(dphi >= 0.0)

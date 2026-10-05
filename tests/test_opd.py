import numpy as np
import pytest
from scipy.integrate import quad

from thermal_opd.opd import steady_opd
from thermal_opd.steady import axial_integral, axial_response, steady_temperature


@pytest.mark.parametrize("lam", [40.0, 300.0, 6000.0])
def test_axial_integral_matches_independent_quadrature(lam):
    args = (lam, 1000.0, 1.38, 0.005, 15.0, 5.0)
    expected, _ = quad(
        lambda z: axial_response(z, *args)[0], 0.0, 0.005,
        epsabs=1e-14, epsrel=1e-11,
    )

    assert np.isclose(axial_integral(*args), expected, rtol=1e-10, atol=0.0)


def test_high_mode_integral_matches_insulated_face_energy_identity():
    # Integrating k*Z''=k*lam^2*Z with front flux C and insulated rear gives
    # integral Z dz = C/(k*lam^2) when both faces have no environmental loss.
    lam, source, conductivity = 150000.0, 1000.0, 1.38
    expected = source / (conductivity * lam**2)
    actual = axial_integral(lam, source, conductivity, 0.005, 0.0, 0.0)

    assert np.isclose(actual, expected, rtol=1e-12, atol=0.0)


def test_opd_matches_independent_temperature_integral_on_broadcast_grid():
    r = np.array([0.0, 0.002])
    phi = np.array([0.0, np.pi / 2])
    coefficients = {0: np.array([1000.0]), 2: np.array([200.0])}
    roots = {0: np.array([40.0]), 2: np.array([300.0])}
    args = (coefficients, roots, 1.38, 0.005, 15.0, 5.0)
    dn_dT = 0.96e-5
    actual = steady_opd(r[:, None], phi[None, :], *args, dn_dT)
    expected = np.empty((len(r), len(phi)))
    for i, ri in enumerate(r):
        for j, phij in enumerate(phi):
            integral, _ = quad(
                lambda z: steady_temperature(ri, phij, z, *args),
                0.0, 0.005, epsabs=1e-13, epsrel=1e-11,
            )
            expected[i, j] = dn_dT * integral

    assert np.allclose(actual, expected, rtol=1e-10, atol=0.0)


def test_opd_preserves_thermo_optic_coefficient_sign():
    args = ({0: np.array([1000.0])}, {0: np.array([40.0])}, 1.38, 0.005, 15.0, 5.0)
    positive = steady_opd(0.0, 0.0, *args, 1e-5)
    negative = steady_opd(0.0, 0.0, *args, -1e-5)

    assert positive > 0
    assert np.isclose(negative, -positive, rtol=1e-12, atol=0.0)

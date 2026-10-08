import numpy as np
import pytest

from thermal_opd import convergence


def test_field_metrics_weight_temperature_and_remove_only_opd_piston():
    # Replacing weighted errors by ordinary means or retaining piston must fail.
    error = convergence.field_errors(
        np.array([[2., 6.]]), np.array([[2., 4.]]), np.array([1., 3.]),
        np.array([11., 14., 17.]), np.array([1., 3., 5.]),
    )
    assert error['temperature_relative_l2'] == pytest.approx([np.sqrt(12 / 52)])
    assert error['temperature_relative_linf'] == pytest.approx([0.5])
    assert error['opd_piston_removed_rms_over_reference_pv'] == pytest.approx(np.sqrt(2/3)/4)
    assert error['opd_piston_removed_linf_over_reference_pv'] == pytest.approx(0.25)
    assert error['opd_relative_pv'] == pytest.approx(0.5)


def test_joint_acceptance_checks_every_field_and_plane():
    good = {'temperature_relative_l2': [0., 0.0001, 0.],
            'temperature_relative_linf': [0., 0., 0.],
            'opd_piston_removed_rms_over_reference_pv': 0.,
            'opd_piston_removed_linf_over_reference_pv': 0., 'opd_relative_pv': 0.}
    assert convergence.within_tolerance(good, 0.001)
    assert not convergence.within_tolerance({**good, 'temperature_relative_linf': [0., 0.002, 0.]}, 0.001)
    assert not convergence.within_tolerance({**good, 'opd_relative_pv': np.nan}, 0.001)


def test_continuous_pv_recovers_quadratic_extrema_between_grid_points():
    result = convergence.continuous_aperture_pv(
        lambda x, y: (x - 0.13)**2 + 2*(y + 0.09)**2,
        radius=1., thickness=0.2, refraction_angle=0., aperture_radius=0.9,
        num_seed_radius=7, num_seed_phi=24,
    )
    assert result['minimum'] == pytest.approx(0., abs=1e-10)
    # Independently maximize the boundary's one-variable polynomial by dense scan.
    phi = np.linspace(0, 2*np.pi, 200001)
    want = np.max((0.9*np.cos(phi)-0.13)**2 + 2*(0.9*np.sin(phi)+0.09)**2)
    assert result['pv'] == pytest.approx(want, rel=1e-8)


@pytest.mark.parametrize('angle', [np.arctan(0.3), -np.arctan(0.3)])
def test_continuous_pv_enforces_rear_transmission_for_signed_angle(angle):
    # R=1,L=1,d=+-0.3 clips y at +-0.7, leaving the other end at +-0.9.
    result = convergence.continuous_aperture_pv(
        lambda x, y: y, radius=1., thickness=1., refraction_angle=angle,
        aperture_radius=0.9, num_seed_radius=7, num_seed_phi=24,
    )
    assert result['pv'] == pytest.approx(1.6, abs=1e-8)
    assert result['maximum'] == pytest.approx(0.7 if angle > 0 else 0.9, abs=1e-8)

"""Compose existing modal functions, with numerical inputs supplied by scripts."""

from dataclasses import dataclass
import numpy as np

from thermal_opd.eigenmodes import radial_eigenvalues
from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.projection import fourier_bessel_coefficients


@dataclass(frozen=True)
class ModalResolution:
    M: int
    N: int
    projection_r: int
    projection_phi: int
    ray_points: int


def solve_modal_source(parameters, angle_deg, resolution):
    """Return source coefficients and radial roots; no temperature interpolation."""
    p, s = parameters, resolution
    if s.M < 0 or s.M % 2 or min(s.N, s.ray_points) < 1 or s.projection_r < 2 or s.projection_phi <= 2*s.M:
        raise ValueError('Require even nonnegative M and adequately sampled projection')
    r = np.linspace(0, p.radius, s.projection_r)
    phi = np.linspace(0, 2*np.pi, s.projection_phi, endpoint=False)
    q = oblique_gaussian_heat_flux(r[:, None]*np.cos(phi), r[:, None]*np.sin(phi),
                                   p.beam_radius, np.deg2rad(angle_deg), p.absorbed_power)
    # Exactly normal incidence is axisymmetric; numerical cosine sums should
    # not create artificial angular modes at roundoff amplitude.
    orders = (0,) if angle_deg == 0 else range(0, s.M + 1, 2)
    roots = {m: radial_eigenvalues(m, p.radius, p.conductivity, p.h_eff, s.N) for m in orders}
    coefficients = {m: fourier_bessel_coefficients(q, r, phi, m, values) for m, values in roots.items()}
    return coefficients, roots

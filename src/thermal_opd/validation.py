"""Steady modal boundary checks and independent surface heat-loss integrals.

These checks establish consistency with the retained source modes. They do
not establish convergence to the prescribed Gaussian or agreement with FEM.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv

from thermal_opd.projection import reconstruct_cosine_field
from thermal_opd.steady import steady_temperature, steady_temperature_gradient


def steady_boundary_residuals(
    r, phi, side_z, coefficients, eigenvalues,
    radius, conductivity, thickness, h_front, h_back, h_side,
):
    """Return sampled maximum absolute front/rear/side residuals (W/m²).

    r, phi and side_z are 1D evaluation coordinates in m, rad and m.
    phi spans a uniform full period without a repeated endpoint. Front z=0:
    h_front*T-k*dT/dz=q_modal; rear z=L: h_back*T+k*dT/dz=0;
    side r=R: h_side*T+k*dT/dr=0. T is the temperature rise.
    """
    r, phi, side_z = np.asarray(r), np.asarray(phi), np.asarray(side_z)
    args = (coefficients, eigenvalues, conductivity, thickness, h_front, h_back)
    front = steady_temperature(r[:, None], phi[None, :], 0.0, *args)
    rear = steady_temperature(r[:, None], phi[None, :], thickness, *args)
    _, front_dz = steady_temperature_gradient(r[:, None], phi[None, :], 0.0, *args)
    _, rear_dz = steady_temperature_gradient(r[:, None], phi[None, :], thickness, *args)
    modal_source = reconstruct_cosine_field(coefficients, eigenvalues, r, phi)
    side = steady_temperature(radius, phi[None, :], side_z[:, None], *args)
    side_dr, _ = steady_temperature_gradient(radius, phi[None, :], side_z[:, None], *args)
    return {
        "front_modal_W_per_m2": float(np.max(np.abs(h_front*front - conductivity*front_dz - modal_source))),
        "rear_W_per_m2": float(np.max(np.abs(h_back*rear + conductivity*rear_dz))),
        "side_W_per_m2": float(np.max(np.abs(h_side*side + conductivity*side_dr))),
    }


def steady_energy_balance(
    coefficients, eigenvalues, radius, conductivity, thickness,
    h_front, h_back, h_side, *, num_r, num_z, num_phi,
):
    """Return input and surface heat-loss powers and absolute imbalance (W).

    Surface losses use Gauss-Legendre integration in r/z and uniform periodic
    angular integration. Input power is the analytic disk integral of the
    retained source: only m=0 survives angular integration. The integrals are
    evaluated separately, rather than enforcing total_loss=input in code.

    Numerical quadrature counts are explicit caller settings. No relative
    error is returned, since input may be zero for zero or purely m>0 sources.
    """
    if num_phi <= 2 * max(coefficients, default=0):
        raise ValueError("num_phi must be greater than twice the maximum azimuthal order")
    nodes_r, weights_r = leggauss(num_r)
    r = (nodes_r + 1) * radius / 2
    weights_r = weights_r * radius / 2
    nodes_z, weights_z = leggauss(num_z)
    z = (nodes_z + 1) * thickness / 2
    weights_z = weights_z * thickness / 2
    phi = np.linspace(0.0, 2*np.pi, num_phi, endpoint=False)
    args = (coefficients, eigenvalues, conductivity, thickness, h_front, h_back)
    front = steady_temperature(r[:, None], phi[None, :], 0.0, *args)
    rear = steady_temperature(r[:, None], phi[None, :], thickness, *args)
    side = steady_temperature(radius, phi[None, :], z[:, None], *args)
    front_loss = h_front * 2*np.pi * np.sum(weights_r * r * np.mean(front, axis=1))
    rear_loss = h_back * 2*np.pi * np.sum(weights_r * r * np.mean(rear, axis=1))
    side_loss = h_side * 2*np.pi * radius * np.sum(weights_z * np.mean(side, axis=1))
    modal_input = sum(
        c * 2*np.pi * radius * jv(1, lam * radius) / lam
        for c, lam in zip(coefficients.get(0, []), eigenvalues.get(0, []))
    )
    total_loss = front_loss + rear_loss + side_loss
    return {
        "front_loss_W": float(front_loss),
        "rear_loss_W": float(rear_loss),
        "side_loss_W": float(side_loss),
        "total_loss_W": float(total_loss),
        "modal_input_W": float(modal_input),
        "balance_absolute_error_W": float(abs(total_loss - modal_input)),
    }

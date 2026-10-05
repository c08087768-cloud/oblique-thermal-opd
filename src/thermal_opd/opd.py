"""Substrate thermo-optic OPD along paths normal to the plane faces."""

import numpy as np
from scipy.special import jv

from thermal_opd.steady import axial_integral


def steady_opd(r, phi, coefficients, eigenvalues, conductivity, thickness, h_front, h_back, dn_dT):
    """Return dn_dT * integral_0^L DeltaT(r,phi,z) dz in metres.

    r and phi are broadcastable coordinates in metres and radians. Modal
    inputs follow steady_temperature. Integration is analytical for each
    axial response; no sampled 3D temperature array is required.

    This normal-path definition includes substrate refractive-index change
    relative to the unheated state, without removing piston. It does not
    describe obliquely refracted rays, deformation, or stress-optic effects.
    """
    r = np.asarray(r)
    phi = np.asarray(phi)
    opd = np.zeros(np.broadcast_shapes(r.shape, phi.shape))
    for m, coeffs_m in coefficients.items():
        angular = np.cos(m * phi)
        for n, source_coefficient in enumerate(coeffs_m):
            lam = eigenvalues[m][n]
            integrated_temperature = axial_integral(
                lam, source_coefficient, conductivity, thickness, h_front, h_back
            )
            opd += dn_dT * jv(m, lam * r) * angular * integrated_temperature
    return opd

"""Shared Fourier-Bessel spatial sums for steady and transient responses."""
import numpy as np
from scipy.special import jv, jvp


def modal_temperature(r, phi, z, coefficients, eigenvalues, response):
    """Shared spatial sum; response(m,n,z) supplies steady or transient amplitudes."""
    r, phi, z = (np.asarray(a) for a in (r, phi, z))
    temperature = np.zeros(np.broadcast_shapes(r.shape,phi.shape,z.shape))
    for m, values in coefficients.items():
        for n in range(len(values)):
            temperature += jv(m, eigenvalues[m][n]*r)*np.cos(m*phi)*response(m, n, z)
    return temperature


def modal_temperature_gradient(r, phi, z, coefficients, eigenvalues, response):
    """Shared analytic spatial derivatives; response returns (amplitude, dz)."""
    r,phi,z = (np.asarray(a) for a in (r,phi,z))
    shape = np.broadcast_shapes(r.shape,phi.shape,z.shape)
    radial_gradient = np.zeros(shape)
    axial_gradient = np.zeros(shape)
    for m, coeffs_m in coefficients.items():
        angular = np.cos(m * phi)
        for n in range(len(coeffs_m)):
            lam = eigenvalues[m][n]
            axial, derivative = response(m,n,z)
            radial_gradient += lam * jvp(m, lam * r) * angular * axial
            axial_gradient += jv(m, lam * r) * angular * derivative
    return radial_gradient, axial_gradient

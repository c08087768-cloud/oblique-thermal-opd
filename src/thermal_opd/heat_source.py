import numpy as np


def oblique_gaussian_heat_flux(x, y, w0, theta, absorbed_power):
    """Return absorbed surface heat flux (W/m^2) for incidence in the y-z plane.

    ``x`` and ``y`` are surface coordinates in metres and support NumPy
    broadcasting. ``w0`` is the 1/e^2 intensity radius (m) in a plane normal
    to propagation; ``theta`` is the incidence angle in radians, with
    abs(theta) < pi/2. The projected footprint's major axis is along y.

    ``absorbed_power`` (W) is the Gaussian integral over the full plane.
    The flux is not renormalized to the finite substrate disk.
    """

    cos_theta = np.cos(theta)

    exponent = -2.0 * (x**2 + (y * cos_theta) ** 2) / w0**2
    q0 = 2.0 * absorbed_power * cos_theta / (np.pi * w0**2)

    return q0 * np.exp(exponent)


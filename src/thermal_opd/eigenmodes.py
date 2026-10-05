import numpy as np
from scipy.optimize import brentq
from scipy.special import jv, jvp


"""
建立 λ 的扫描区间
        ↓
计算 f(λ)
        ↓
找 f(λ_i) f(λ_{i+1}) < 0
        ↓
说明中间穿过 0
        ↓
brentq 精确求这个根
        ↓
收集前 N 个 λ_mn
"""

def radial_eigenvalue_equation(lam, m, radius, conductivity, h):
    """Return k*lam*J_m'(lam*R) + h*J_m(lam*R) for sidewall heat loss.

    ``lam`` is a radial wavenumber (1/m), ``conductivity`` is k (W/(m K)),
    and ``h`` is the sidewall heat-transfer coefficient (W/(m^2 K)).
    """

    x = lam * radius

    return (
        conductivity * lam * jvp(m, x)
        + h * jv(m, x)
    )


def radial_eigenvalues(m, radius, conductivity, h, num_modes):
    """Return the first positive radial wavenumbers (1/m) in ascending order.

    The current search is used for nonnegative integer m, positive radius,
    conductivity and h, and a positive number of modes. For h=0 the m=0
    zero mode requires separate treatment and is not included here.
    """

    def equation(lam):
        return radial_eigenvalue_equation(
            lam, m, radius, conductivity, h
        )

    roots = []

    # Work with the dimensionless variable x = lambda * R.
    x_max = (num_modes + m / 2 + 5) * np.pi
    x_grid = np.linspace(1e-8, x_max, 10000)

    values = [equation(x / radius) for x in x_grid]

    for i in range(len(x_grid) - 1):
        if values[i] * values[i + 1] < 0:
            lam_left = x_grid[i] / radius
            lam_right = x_grid[i + 1] / radius

            root = brentq(equation, lam_left, lam_right)
            roots.append(root)

            if len(roots) == num_modes:
                break

    if len(roots) < num_modes:
        raise RuntimeError(
            f"Found only {len(roots)} roots, expected {num_modes}."
        )

    return np.array(roots)


def radial_mode_norm(lam, m, radius):
    """Return integral_0^R J_m(lam*r)^2 r dr analytically (m^2).

    The Bessel integral identity holds for positive lam, including the
    Robin eigenvalues. At those roots it is equivalent to the norm in
    supplemental Eq. (S12), without requiring h as another input.
    """
    x = lam * radius
    return 0.5 * radius**2 * (
        jv(m, x)**2 - jv(m - 1, x) * jv(m + 1, x)
    )

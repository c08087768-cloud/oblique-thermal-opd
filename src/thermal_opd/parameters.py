"""Physical inputs for the manuscript, expressed in SI units.

Numerical settings (modal truncation and projection grids) belong in the
calling script. The incidence angle is also supplied by the calling script.
"""

from dataclasses import dataclass

from scipy.constants import Stefan_Boltzmann


@dataclass(frozen=True)
class PhysicalParameters:
    """Nominal geometry, material, environment, and optical inputs.

    ``absorbed_power`` is assigned to the Gaussian on the full plane;
    power intercepted by the finite substrate must be evaluated separately.
    ``dn_dT`` retains the manuscript's value tabulated at 643.8 nm; it is
    not a measurement at the specified 1064 nm optical wavelength.
    """

    radius: float = 12.7e-3              # m
    thickness: float = 5.0e-3            # m
    conductivity: float = 1.38           # W/(m K)
    h_conv: float = 10.0                 # W/(m^2 K)
    emissivity: float = 0.85
    ambient_temperature: float = 300.0   # K
    absorbed_power: float = 110e-3       # W
    beam_radius_ratio: float = 0.30
    dn_dT: float = 0.96e-5               # 1/K, tabulated at 643.8 nm
    wavelength: float = 1064e-9          # m

    @property
    def beam_radius(self) -> float:
        """Return the 1/e^2 intensity radius normal to beam propagation (m)."""
        return self.beam_radius_ratio * self.radius

    @property
    def h_eff(self) -> float:
        """Return convection plus linearized radiation at ambient (W/(m^2 K))."""
        return (
            self.h_conv + 4 * self.emissivity * Stefan_Boltzmann
            * self.ambient_temperature**3
        )

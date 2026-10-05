"""Aperture selection and optical-path-difference metrics."""

import numpy as np


def circular_aperture_mask(x, y, aperture_radius):
    """Return a boundary-inclusive circular mask for broadcastable x,y (m)."""
    if not np.isfinite(aperture_radius) or aperture_radius <= 0:
        raise ValueError("aperture_radius must be finite and positive")
    x, y = np.broadcast_arrays(x, y)
    return x**2 + y**2 <= aperture_radius**2


def opd_peak_to_valley(opd, aperture_mask):
    """Return max(OPD)-min(OPD) on the supplied aperture, in the input units.

    The mask must have the same shape as opd. All selected values must be
    finite; values outside the aperture are excluded, including NaNs.
    """
    opd = np.asarray(opd)
    aperture_mask = np.asarray(aperture_mask, dtype=bool)
    if opd.shape != aperture_mask.shape:
        raise ValueError("opd and aperture_mask must have the same shape")
    selected = opd[aperture_mask]
    if selected.size == 0:
        raise ValueError("aperture is empty")
    if not np.all(np.isfinite(selected)):
        raise ValueError("OPD values inside aperture must be finite")
    return float(np.max(selected) - np.min(selected))

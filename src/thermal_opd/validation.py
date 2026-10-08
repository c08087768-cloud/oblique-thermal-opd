"""Modal boundary checks and independent steady surface heat-loss integrals.

These checks establish consistency with the retained source modes. They do
not establish convergence to the prescribed Gaussian or agreement with FEM.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv

from thermal_opd.projection import reconstruct_cosine_field
from thermal_opd.steady import steady_temperature, steady_temperature_gradient


def transient_boundary_fields(model, times, r, phi, side_z):
    """Sample signed Robin residuals against g(t)*q_modal, in W/m².

    Boundary temperatures and analytic normal-coordinate derivatives are saved
    separately so validators can reconstruct residuals. This checks the retained
    model, not source projection accuracy or transient PDE/J convergence.
    At t=0 both T and the history's load are zero, including a step convention.
    """
    p = model.parameters
    times,r,phi,side_z = (np.asarray(a,dtype=float) for a in (times,r,phi,side_z))
    if any(a.ndim!=1 or not a.size or not np.all(np.isfinite(a)) for a in (times,r,phi,side_z)):
        raise ValueError('Require finite nonempty 1D boundary coordinates/times')
    if times[0]<0 or np.any(np.diff(times)<=0) or r[0]<0 or r[-1]>p.radius or np.any(np.diff(r)<=0):
        raise ValueError('Invalid boundary times/radii')
    if phi[0]!=0 or phi[-1]>=2*np.pi or np.any(np.diff(phi)<=0):
        raise ValueError('Invalid boundary azimuths')
    if side_z[0]!=0 or side_z[-1]!=p.thickness or np.any(np.diff(side_z)<=0):
        raise ValueError('Side boundary must cover z=0 to L')
    source = reconstruct_cosine_field(model.coefficients,model.eigenvalues,r,phi)
    front,rear,side,front_dz,rear_dz,side_dr = [],[],[],[],[],[]
    for t in times:
        front.append(model.temperature(r[:,None],phi[None,:],0.,t))
        rear.append(model.temperature(r[:,None],phi[None,:],p.thickness,t))
        side.append(model.temperature(p.radius,phi[None,:],side_z[:,None],t))
        front_dz.append(model.temperature_gradient(r[:,None],phi[None,:],0.,t)[1])
        rear_dz.append(model.temperature_gradient(r[:,None],phi[None,:],p.thickness,t)[1])
        side_dr.append(model.temperature_gradient(p.radius,phi[None,:],side_z[:,None],t)[0])
    f = dict(times_s=times,r_m=r,phi_rad=phi,side_z_m=side_z,modal_source_W_per_m2=source,
        front_temperature_K=np.array(front),rear_temperature_K=np.array(rear),side_temperature_K=np.array(side),
        front_dT_dz_K_per_m=np.array(front_dz),rear_dT_dz_K_per_m=np.array(rear_dz),side_dT_dr_K_per_m=np.array(side_dr),
        load_factors=np.array([model.history.value(t) for t in times]))
    f['front_load_W_per_m2'] = f['load_factors'][:,None,None]*source
    f.update(transient_residual_fields(f,p))
    return f


def transient_residual_fields(f,p):
    """Reconstruct signed residuals from boundary T, derivatives and load."""
    return dict(
        front_modal_residual_W_per_m2=p.h_eff*f['front_temperature_K']-p.conductivity*f['front_dT_dz_K_per_m']-f['front_load_W_per_m2'],
        rear_residual_W_per_m2=p.h_eff*f['rear_temperature_K']+p.conductivity*f['rear_dT_dz_K_per_m'],
        side_residual_W_per_m2=p.h_eff*f['side_temperature_K']+p.conductivity*f['side_dT_dr_K_per_m'])


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

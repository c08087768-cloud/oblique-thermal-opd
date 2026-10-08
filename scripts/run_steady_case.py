"""Run and save one steady-temperature / refracted-ray thermo-optic OPD case."""

from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from thermal_opd.convergence import continuous_aperture_pv

from thermal_opd.eigenmodes import radial_eigenvalues
from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.metrics import circular_aperture_mask, opd_peak_to_valley
from thermal_opd.opd import refracted_angle, transmitted_ray_mask, steady_refracted_opd
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.projection import fourier_bessel_coefficient
from thermal_opd.steady import steady_temperature, steady_temperature_gradient
from thermal_opd.validation import steady_boundary_residuals, steady_energy_balance
from thermal_opd.plotting import save_temperature_opd_figure


parameters = PhysicalParameters()
incidence_angle_deg = 60.0
aperture_radius_ratio = 0.90

# Jointly verified at 0/30/60 deg against enlarged finite references, 0.1% target.
# Scope: nominal parameters and w0/R=0.30. Keep 201x201 for readable field plots
# (121x121 already met the sampling criterion). Revalidate changed physics.
max_azimuthal_order = 6
num_radial_modes = 8
num_projection_r = 400
num_projection_phi = 32
num_temperature_r = 161
num_temperature_phi = 180
num_opd_xy = 201
num_opd_ray_points = 16
num_pv_seed_radius = 9
num_pv_seed_phi = 48
opd_point_chunk_size = 1024
num_boundary_z = 65
num_energy_r = 80
num_energy_z = 64
num_energy_phi = 180

repository_root = Path(__file__).resolve().parents[1]


def main():
    p = parameters
    if num_projection_phi <= 2 * max_azimuthal_order:
        raise ValueError("Angular projection requires more than 2*M samples")
    if not 0 < aperture_radius_ratio <= 1:
        raise ValueError("The evaluation aperture must lie inside the substrate")
    started = datetime.now(timezone.utc)
    run_id = f"steady_refracted_theta{incidence_angle_deg:g}_M{max_azimuthal_order}_N{num_radial_modes}_{started:%Y%m%dT%H%M%S%fZ}"
    result_directory = repository_root / "results" / run_id
    figure_directory = repository_root / "outputs/figures" / run_id
    print(f"Run: {run_id}", flush=True)
    print(f"R={p.radius*1e3:g} mm, L={p.thickness*1e3:g} mm, w0={p.beam_radius*1e3:g} mm", flush=True)
    print(f"Power on full plane={p.absorbed_power*1e3:g} mW, h_eff={p.h_eff:.9f} W/(m2 K)", flush=True)

    # 1. Heat-flux coefficients. No renormalization of intercepted disk power.
    projection_r = np.linspace(0.0, p.radius, num_projection_r)
    projection_phi = np.linspace(0.0, 2*np.pi, num_projection_phi, endpoint=False)
    theta = np.deg2rad(incidence_angle_deg)
    internal_angle = refracted_angle(theta, p.incident_refractive_index, p.substrate_refractive_index)
    print(f"Refraction angle={np.rad2deg(internal_angle):.6f} deg; rear y shift={p.thickness*np.tan(internal_angle)*1e3:.6f} mm", flush=True)
    source = oblique_gaussian_heat_flux(
        projection_r[:, None]*np.cos(projection_phi),
        projection_r[:, None]*np.sin(projection_phi),
        p.beam_radius, theta, p.absorbed_power,
    )
    print("Computing radial roots and source projections", flush=True)
    roots = {
        m: radial_eigenvalues(m, p.radius, p.conductivity, p.h_eff, num_radial_modes)
        for m in range(0, max_azimuthal_order + 1, 2)
    }
    coefficients = {
        m: np.array([fourier_bessel_coefficient(source, projection_r, projection_phi, m, lam)
                     for lam in values])
        for m, values in roots.items()
    }
    thermal_args = (coefficients, roots, p.conductivity, p.thickness, p.h_eff, p.h_eff)

    # 2. Temperature rises, computed directly from retained coefficients.
    r = np.linspace(0.0, p.radius, num_temperature_r)
    phi = np.linspace(0.0, 2*np.pi, num_temperature_phi, endpoint=False)
    z = np.array([0.0, p.thickness / 2, p.thickness])
    temperatures = np.stack([
        steady_temperature(r[:, None], phi[None, :], plane, *thermal_args) for plane in z
    ])

    # 3. Refracted rays labelled by entrance coordinates; raw OPD keeps piston.
    aperture_radius = aperture_radius_ratio * p.radius
    axis = np.linspace(-aperture_radius, aperture_radius, num_opd_xy)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    entrance_mask = circular_aperture_mask(x, y, aperture_radius)
    mask = entrance_mask & transmitted_ray_mask(x, y, p.radius, p.thickness, internal_angle)
    opd = np.full(x.shape, np.nan)
    print("Computing temperature and aperture OPD", flush=True)
    opd[mask] = steady_refracted_opd(
        x[mask], y[mask], *thermal_args, p.dn_dT, radius=p.radius,
        refraction_angle=internal_angle, num_ray_points=num_opd_ray_points,
        point_chunk_size=opd_point_chunk_size,
    )
    pv = opd_peak_to_valley(opd, mask)
    def optical_at(entry_x, entry_y):
        return steady_refracted_opd(entry_x, entry_y, *thermal_args, p.dn_dT, radius=p.radius,
                                    refraction_angle=internal_angle, num_ray_points=num_opd_ray_points,
                                    point_chunk_size=opd_point_chunk_size)
    continuous_extrema = continuous_aperture_pv(
        optical_at, radius=p.radius, thickness=p.thickness, refraction_angle=internal_angle,
        aperture_radius=aperture_radius, num_seed_radius=num_pv_seed_radius, num_seed_phi=num_pv_seed_phi,
    )

    # 4. Shared modal boundary / energy checks, independent of OPD evaluation.
    print("Checking Robin boundaries and integrated heat loss", flush=True)
    residuals = steady_boundary_residuals(
        r, phi, np.linspace(0.0, p.thickness, num_boundary_z),
        coefficients, roots, p.radius, p.conductivity, p.thickness,
        p.h_eff, p.h_eff, p.h_eff,
    )
    balance = steady_energy_balance(
        coefficients, roots, p.radius, p.conductivity, p.thickness,
        p.h_eff, p.h_eff, p.h_eff,
        num_r=num_energy_r, num_z=num_energy_z, num_phi=num_energy_phi,
    )
    # Compare the modal front flux with the prescribed source as a distinct check.
    # The difference includes finite truncation/projection, not just a BC residual.
    _, front_dz = steady_temperature_gradient(r[:, None], phi[None, :], 0.0, *thermal_args)
    specified = oblique_gaussian_heat_flux(
        r[:, None]*np.cos(phi), r[:, None]*np.sin(phi),
        p.beam_radius, theta, p.absorbed_power,
    )
    peak_flux = float(specified.max())
    front_gaussian_difference = float(np.max(np.abs(p.h_eff*temperatures[0] - p.conductivity*front_dz - specified)))
    # Independent source integral: Gauss radial nodes, periodic angular samples.
    nodes, weights = leggauss(num_energy_r)
    energy_r = (nodes + 1)*p.radius/2
    energy_phi = np.linspace(0.0, 2*np.pi, num_energy_phi, endpoint=False)
    gaussian = oblique_gaussian_heat_flux(
        energy_r[:, None]*np.cos(energy_phi), energy_r[:, None]*np.sin(energy_phi),
        p.beam_radius, theta, p.absorbed_power,
    )
    disk_input = float(2*np.pi*np.sum(weights*p.radius/2*energy_r*np.mean(gaussian, axis=1)))
    energy_relative_error = (balance["balance_absolute_error_W"] / abs(balance["modal_input_W"])
                             if balance["modal_input_W"] != 0 else None)
    normalized_residuals = {key: value / peak_flux if peak_flux != 0 else None
                            for key, value in residuals.items()}

    # 5. Save the settings, coefficients and fields needed to reproduce this run.
    numerical_settings = {
        "M": max_azimuthal_order, "N": num_radial_modes,
        "projection_r": num_projection_r, "projection_phi": num_projection_phi,
        "temperature_r": num_temperature_r, "temperature_phi": num_temperature_phi,
        "opd_xy": num_opd_xy, "opd_ray_points": num_opd_ray_points,
        "pv_seed_radius": num_pv_seed_radius, "pv_seed_phi": num_pv_seed_phi,
        "opd_point_chunk_size": opd_point_chunk_size, "boundary_z": num_boundary_z,
        "energy_r": num_energy_r, "energy_z": num_energy_z, "energy_phi": num_energy_phi,
    }
    source_files = [Path(__file__).resolve(), *sorted((repository_root / "src/thermal_opd").glob("*.py"))]
    summary = {
        "kind": "complete_representative_steady_refracted_ray_case",
        "run_id": run_id, "started_utc": started.isoformat(),
        "physical_parameters_SI": asdict(p), "h_eff_W_per_m2_K": p.h_eff,
        "incidence_angle_deg": incidence_angle_deg, "numerical_settings": numerical_settings,
        "coordinates": "front z=0; rear z=L; z_code=z_paper+L/2; heating major axis y",
        "units": {"coordinates": "m", "temperature_rise": "K", "opd": "m", "flux": "W/m2", "power": "W"},
        "opd_definition": "dn_dT / cos(theta_t) * integral_0^L DeltaT(x0, y0+z*tan(theta_t), z) dz; entrance-plane coordinates; raw piston retained",
        "ray_geometry": {
            "refraction_angle_rad": internal_angle,
            "refraction_angle_deg": float(np.rad2deg(internal_angle)),
            "rear_y_shift_m": float(p.thickness*np.tan(internal_angle)),
            "path_length_m": float(p.thickness/np.cos(internal_angle)),
            "reference": "straight ray in unheated uniform substrate; first-order thermo-optic term",
        },
        "aperture_radius_ratio": aperture_radius_ratio, "aperture_radius_m": aperture_radius,
        "aperture_sample_count": int(mask.sum()),
        "entrance_aperture_sample_count": int(entrance_mask.sum()),
        "excluded_side_exit_sample_count": int(np.count_nonzero(entrance_mask & ~mask)),
        "center_temperature_rises_K": temperatures[:, 0, 0].tolist(),
        "sampled_peak_temperature_rises_K": temperatures.max(axis=(1, 2)).tolist(),
        "opd_pv_m": pv,
        "opd_pv_convention": "sampled PV on saved 201x201 grid; continuous numerical extrema reported separately",
        "opd_continuous_numerical_extrema_m": continuous_extrema,
        "truncation_selection_scope": "0/30/60 deg; nominal material and w0/R=0.30; 0.1% field/PV numerical-change target; finite references",
        "boundary_absolute_residuals_W_per_m2": residuals,
        "boundary_residuals_over_peak_gaussian_flux": normalized_residuals,
        "front_flux_vs_gaussian_max_difference_W_per_m2": front_gaussian_difference,
        "front_flux_vs_gaussian_difference_over_peak_flux": front_gaussian_difference / peak_flux if peak_flux != 0 else None,
        "energy_balance": balance, "modal_energy_balance_relative_error": energy_relative_error,
        "prescribed_disk_input_W": disk_input,
        "modal_input_relative_difference_vs_gaussian": abs(balance["modal_input_W"] - disk_input) / abs(disk_input) if disk_input != 0 else None,
        "eigenvalues_per_m": {str(m): values.tolist() for m, values in roots.items()},
        "source_coefficients_W_per_m2": {str(m): values.tolist() for m, values in coefficients.items()},
        "source_sha256": {file.relative_to(repository_root).as_posix(): hashlib.sha256(file.read_bytes()).hexdigest() for file in source_files},
        "limits": "Single angle; finite numerical settings; straight refracted reference rays; only rear-face transmission in the entrance aperture; thermo-optic contribution; no thermal ray bending or FEM comparison",
        "artifacts": {"fields": "fields.npz", "figure": (figure_directory / "steady_case.png").relative_to(repository_root).as_posix()},
    }
    result_directory.mkdir(parents=True, exist_ok=False)
    (result_directory / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    np.savez_compressed(
        result_directory / "fields.npz", temperature_r_m=r, temperature_phi_rad=phi,
        temperature_z_m=z, temperature_rise_K=temperatures,
        opd_x_m=axis, opd_y_m=axis, entrance_aperture_mask=entrance_mask,
        aperture_mask=mask, opd_m=opd,
    )

    figure_directory.mkdir(parents=True, exist_ok=False)
    save_temperature_opd_figure(figure_directory/'steady_case.png',r,phi,temperatures,x,y,mask,opd,
        title=f'Steady case: {incidence_angle_deg:g} deg, M={max_azimuthal_order}, N={num_radial_modes}',
        pv_m=pv,continuous_pv_m=continuous_extrema['pv'])
    print(f"Center temperature rises [K]: {summary['center_temperature_rises_K']}", flush=True)
    print(f"Refracted-ray aperture OPD-PV: {pv*1e9:.6f} nm", flush=True)
    print(f"Continuous numerical OPD-PV: {continuous_extrema['pv']*1e9:.6f} nm", flush=True)
    print(f"Rear-transmitted samples: {mask.sum()} / {entrance_mask.sum()}", flush=True)
    print(f"Disk / modal input: {disk_input*1e3:.6f} / {balance['modal_input_W']*1e3:.6f} mW", flush=True)
    print(f"Modal energy balance relative error: {energy_relative_error}", flush=True)
    print(f"Boundary residuals / peak flux: {normalized_residuals}", flush=True)
    print(f"Results: {result_directory}\nFigure: {figure_directory / 'steady_case.png'}", flush=True)


if __name__ == "__main__":
    main()

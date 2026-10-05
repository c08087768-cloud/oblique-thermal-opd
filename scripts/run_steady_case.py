"""Run and save one complete steady-temperature / normal-path OPD case."""

from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss

from thermal_opd.eigenmodes import radial_eigenvalues
from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.metrics import circular_aperture_mask, opd_peak_to_valley
from thermal_opd.opd import steady_opd
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.projection import fourier_bessel_coefficient
from thermal_opd.steady import steady_temperature, steady_temperature_gradient
from thermal_opd.validation import steady_boundary_residuals, steady_energy_balance


parameters = PhysicalParameters()
incidence_angle_deg = 60.0
aperture_radius_ratio = 0.90

# Explicit settings for this representative case, following the first sweep.
# They are not a universal convergence guarantee for arbitrary physical inputs.
max_azimuthal_order = 10
num_radial_modes = 20
num_projection_r = 1200
num_projection_phi = 256
num_temperature_r = 161
num_temperature_phi = 180
num_opd_xy = 201
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
    run_id = f"steady_theta{incidence_angle_deg:g}_M{max_azimuthal_order}_N{num_radial_modes}_{started:%Y%m%dT%H%M%S%fZ}"
    result_directory = repository_root / "results" / run_id
    figure_directory = repository_root / "outputs/figures" / run_id
    print(f"Run: {run_id}", flush=True)
    print(f"R={p.radius*1e3:g} mm, L={p.thickness*1e3:g} mm, w0={p.beam_radius*1e3:g} mm", flush=True)
    print(f"Power on full plane={p.absorbed_power*1e3:g} mW, h_eff={p.h_eff:.9f} W/(m2 K)", flush=True)

    # 1. Heat-flux coefficients. No renormalization of intercepted disk power.
    projection_r = np.linspace(0.0, p.radius, num_projection_r)
    projection_phi = np.linspace(0.0, 2*np.pi, num_projection_phi, endpoint=False)
    theta = np.deg2rad(incidence_angle_deg)
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

    # 3. Normal-path thermo-optic OPD. Raw stored values retain their piston.
    aperture_radius = aperture_radius_ratio * p.radius
    axis = np.linspace(-aperture_radius, aperture_radius, num_opd_xy)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    mask = circular_aperture_mask(x, y, aperture_radius)
    opd = np.full(x.shape, np.nan)
    print("Computing temperature and aperture OPD", flush=True)
    opd[mask] = steady_opd(np.hypot(x[mask], y[mask]), np.arctan2(y[mask], x[mask]), *thermal_args, p.dn_dT)
    pv = opd_peak_to_valley(opd, mask)

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
        "opd_xy": num_opd_xy, "boundary_z": num_boundary_z,
        "energy_r": num_energy_r, "energy_z": num_energy_z, "energy_phi": num_energy_phi,
    }
    source_files = [Path(__file__).resolve(), *sorted((repository_root / "src/thermal_opd").glob("*.py"))]
    summary = {
        "kind": "complete_representative_steady_case",
        "run_id": run_id, "started_utc": started.isoformat(),
        "physical_parameters_SI": asdict(p), "h_eff_W_per_m2_K": p.h_eff,
        "incidence_angle_deg": incidence_angle_deg, "numerical_settings": numerical_settings,
        "coordinates": "front z=0; rear z=L; z_code=z_paper+L/2; heating major axis y",
        "units": {"coordinates": "m", "temperature_rise": "K", "opd": "m", "flux": "W/m2", "power": "W"},
        "opd_definition": "dn_dT * integral_0^L temperature_rise dz; normal path; raw values relative to unheated substrate",
        "aperture_radius_ratio": aperture_radius_ratio, "aperture_radius_m": aperture_radius,
        "aperture_sample_count": int(mask.sum()),
        "center_temperature_rises_K": temperatures[:, 0, 0].tolist(),
        "sampled_peak_temperature_rises_K": temperatures.max(axis=(1, 2)).tolist(),
        "opd_pv_m": pv,
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
        "limits": "Single representative angle; finite numerical settings; normal-path thermo-optic contribution; no FEM comparison or universal accuracy claim",
        "artifacts": {"fields": "fields.npz", "figure": (figure_directory / "steady_case.png").relative_to(repository_root).as_posix()},
    }
    result_directory.mkdir(parents=True, exist_ok=False)
    (result_directory / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    np.savez_compressed(
        result_directory / "fields.npz", temperature_r_m=r, temperature_phi_rad=phi,
        temperature_z_m=z, temperature_rise_K=temperatures,
        opd_x_m=axis, opd_y_m=axis, aperture_mask=mask, opd_m=opd,
    )

    # One figure: three temperature planes and the OPD on its smaller aperture.
    plot_phi = np.append(phi, 2*np.pi)
    temperature_x = r[:, None]*np.cos(plot_phi)*1e3
    temperature_y = r[:, None]*np.sin(plot_phi)*1e3
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), layout="constrained")
    levels = np.linspace(temperatures.min(), temperatures.max(), 61)
    for ax, field, title in zip(axes.flat, temperatures, ("Front: z=0", "Mid-plane: z=L/2", "Rear: z=L")):
        contour = ax.contourf(temperature_x, temperature_y, np.column_stack((field, field[:, 0])), levels=levels, cmap="inferno")
        ax.set_title(title)
        fig.colorbar(contour, ax=ax, label="Temperature rise [K]", shrink=0.85)
    display_offset = float(np.min(opd[mask]))
    contour = axes.flat[-1].contourf(x*1e3, y*1e3, (opd-display_offset)*1e9, levels=61, cmap="viridis")
    axes.flat[-1].set_title(f"Normal-path OPD: {aperture_radius_ratio:.2f}R aperture, PV={pv*1e9:.4f} nm")
    fig.colorbar(contour, ax=axes.flat[-1], label="OPD minus aperture minimum [nm]", shrink=0.85)
    for ax in axes.flat:
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
    fig.suptitle(f"Steady case: {incidence_angle_deg:g} deg, M={max_azimuthal_order}, N={num_radial_modes}")
    figure_directory.mkdir(parents=True, exist_ok=False)
    fig.savefig(figure_directory / "steady_case.png", dpi=180)
    plt.close(fig)
    print(f"Center temperature rises [K]: {summary['center_temperature_rises_K']}", flush=True)
    print(f"Normal-path aperture OPD-PV: {pv*1e9:.6f} nm", flush=True)
    print(f"Disk / modal input: {disk_input*1e3:.6f} / {balance['modal_input_W']*1e3:.6f} mW", flush=True)
    print(f"Modal energy balance relative error: {energy_relative_error}", flush=True)
    print(f"Boundary residuals / peak flux: {normalized_residuals}", flush=True)
    print(f"Results: {result_directory}\nFigure: {figure_directory / 'steady_case.png'}", flush=True)


if __name__ == "__main__":
    main()

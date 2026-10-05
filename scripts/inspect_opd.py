"""Inspect normal-path thermo-optic OPD with the manuscript's physical inputs."""

from dataclasses import asdict
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np

from thermal_opd.eigenmodes import radial_eigenvalues
from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.metrics import circular_aperture_mask, opd_peak_to_valley
from thermal_opd.opd import steady_opd
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.projection import fourier_bessel_coefficient


parameters = PhysicalParameters()
incidence_angles_deg = (0.0, 30.0, 60.0)
aperture_radius_ratio = 0.90

# Numerical settings for this inspection, defined in the calling code.
max_azimuthal_order = 8
num_radial_modes = 12
num_projection_r = 400
num_projection_phi = 256
num_opd_xy = 201  # Odd count includes the center and axis/aperture intersections.
num_profile_points = 401

# A small sensitivity comparison, not the full convergence study.
finer_num_projection_r = 800
higher_max_azimuthal_order = 10
higher_num_radial_modes = 16

repository_root = Path(__file__).resolve().parents[1]
figure_path = repository_root / "outputs/figures/thermo_optic_opd.png"
result_directory = repository_root / "results/thermo_optic_opd"


def source_coefficients(angle_deg, eigenvalues, num_r, num_phi):
    p = parameters
    r = np.linspace(0.0, p.radius, num_r)
    phi = np.linspace(0.0, 2 * np.pi, num_phi, endpoint=False)
    q = oblique_gaussian_heat_flux(
        r[:, None] * np.cos(phi), r[:, None] * np.sin(phi),
        p.beam_radius, np.deg2rad(angle_deg), p.absorbed_power,
    )
    return {
        m: np.array([fourier_bessel_coefficient(q, r, phi, m, lam) for lam in roots])
        for m, roots in eigenvalues.items()
    }


def opd_at(r, phi, coefficients, eigenvalues):
    p = parameters
    return steady_opd(
        r, phi, coefficients, eigenvalues,
        p.conductivity, p.thickness, p.h_eff, p.h_eff, p.dn_dT,
    )


def main():
    p = parameters
    aperture_radius = aperture_radius_ratio * p.radius
    axis = np.linspace(-aperture_radius, aperture_radius, num_opd_xy)
    X, Y = np.meshgrid(axis, axis, indexing="xy")
    mask = circular_aperture_mask(X, Y, aperture_radius)
    # Evaluate only inside the aperture, which is entirely inside the substrate.
    r = np.hypot(X[mask], Y[mask])
    phi = np.arctan2(Y[mask], X[mask])
    eigenvalues = {
        m: radial_eigenvalues(m, p.radius, p.conductivity, p.h_eff, num_radial_modes)
        for m in range(0, max_azimuthal_order + 1, 2)
    }
    records = []
    opd_maps = []
    case_coefficients = []
    print("Normal-path substrate thermo-optic OPD; values computed in metres")
    print(f"Aperture: {aperture_radius_ratio:.2f} R = {aperture_radius * 1e3:.3f} mm")
    print(f"M={max_azimuthal_order}, N={num_radial_modes}, projection={num_projection_r}x{num_projection_phi}")
    print(f"Cartesian sampling: {num_opd_xy}x{num_opd_xy}, {mask.sum()} aperture points")
    print("angle [deg]    center OPD [nm]    aperture OPD-PV [nm]")
    for angle_deg in incidence_angles_deg:
        coefficients = source_coefficients(angle_deg, eigenvalues, num_projection_r, num_projection_phi)
        field = np.full(X.shape, np.nan)
        field[mask] = opd_at(r, phi, coefficients, eigenvalues)
        pv = opd_peak_to_valley(field, mask)
        center = float(opd_at(0.0, 0.0, coefficients, eigenvalues))
        records.append({
            "incidence_angle_deg": angle_deg,
            "center_opd_m": center,
            "opd_pv_m": pv,
            "source_coefficients_W_per_m2": {str(m): c.tolist() for m, c in coefficients.items()},
        })
        opd_maps.append(field)
        case_coefficients.append(coefficients)
        print(f"{angle_deg:11.1f}    {center*1e9:15.6f}    {pv*1e9:20.6f}")

    # Compare the last case on the same Cartesian aperture coordinates.
    selected_angle = incidence_angles_deg[-1]
    selected_field = opd_maps[-1]
    selected_coefficients = case_coefficients[-1]
    base_pv = records[-1]["opd_pv_m"]
    fine_coefficients = source_coefficients(
        selected_angle, eigenvalues, finer_num_projection_r, num_projection_phi
    )
    fine_values = opd_at(r, phi, fine_coefficients, eigenvalues)
    fine_pv = float(np.ptp(fine_values))
    higher_roots = {
        m: radial_eigenvalues(m, p.radius, p.conductivity, p.h_eff, higher_num_radial_modes)
        for m in range(0, higher_max_azimuthal_order + 1, 2)
    }
    higher_coefficients = source_coefficients(
        selected_angle, higher_roots, finer_num_projection_r, num_projection_phi
    )
    higher_values = opd_at(r, phi, higher_coefficients, higher_roots)
    higher_pv = float(np.ptp(higher_values))
    sensitivity = {
        "incidence_angle_deg": selected_angle,
        "evaluation_xy": num_opd_xy,
        "base_pv_m": base_pv,
        "fine_projection": {
            "M": max_azimuthal_order, "N": num_radial_modes,
            "projection_r": finer_num_projection_r, "projection_phi": num_projection_phi,
            "pv_m": fine_pv, "relative_pv_difference": abs(fine_pv - base_pv) / abs(fine_pv),
        },
        "higher_modes_and_fine_projection": {
            "M": higher_max_azimuthal_order, "N": higher_num_radial_modes,
            "projection_r": finer_num_projection_r, "projection_phi": num_projection_phi,
            "pv_m": higher_pv, "relative_pv_difference": abs(higher_pv - base_pv) / abs(higher_pv),
        },
    }
    print(f"Fine-projection OPD-PV: {fine_pv*1e9:.6f} nm; relative change={sensitivity['fine_projection']['relative_pv_difference']:.6e}")
    print(f"Higher-mode OPD-PV: {higher_pv*1e9:.6f} nm; relative change={sensitivity['higher_modes_and_fine_projection']['relative_pv_difference']:.6e}")

    # Remove one constant only for display; PV and stored fields use raw OPD.
    display_reference = float(np.min(selected_field[mask]))
    profile_axis = np.linspace(-aperture_radius, aperture_radius, num_profile_points)
    x_profile = opd_at(np.abs(profile_axis), np.where(profile_axis < 0, np.pi, 0.0), selected_coefficients, eigenvalues)
    y_profile = opd_at(np.abs(profile_axis), np.where(profile_axis < 0, -np.pi / 2, np.pi / 2), selected_coefficients, eigenvalues)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    contour = axes[0].contourf(
        X*1e3, Y*1e3, (selected_field-display_reference)*1e9,
        levels=np.linspace(0, base_pv*1e9, 61), cmap="viridis",
    )
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x [mm]")
    axes[0].set_ylabel("y [mm]")
    axes[0].set_title("OPD over the 0.90R aperture")
    fig.colorbar(contour, ax=axes[0], label="OPD minus aperture minimum [nm]", shrink=0.85)
    axes[1].plot(profile_axis*1e3, (x_profile-display_reference)*1e9, label="x axis (minor heating axis)")
    axes[1].plot(profile_axis*1e3, (y_profile-display_reference)*1e9, label="y axis (major heating axis)")
    axes[1].set_xlabel("Position along axis [mm]")
    axes[1].set_ylabel("OPD minus aperture minimum [nm]")
    axes[1].set_title("Normal-path OPD profiles")
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    fig.suptitle(f"Thermo-optic OPD: {selected_angle:g} deg, PV={base_pv*1e9:.3f} nm")
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    summary = {
        "kind": "normal_path_thermo_optic_opd_inspection",
        "physical_parameters_SI": asdict(p),
        "h_eff_W_per_m2_K": p.h_eff,
        "opd_definition": "dn_dT * integral_0^L temperature_rise dz; normal path; unit m",
        "aperture_radius_ratio": aperture_radius_ratio,
        "aperture_radius_m": aperture_radius,
        "numerical_settings": {
            "M": max_azimuthal_order, "N": num_radial_modes,
            "projection_r": num_projection_r, "projection_phi": num_projection_phi,
            "opd_xy": num_opd_xy, "profile_points": num_profile_points,
        },
        "eigenvalues_per_m": {str(m): values.tolist() for m, values in eigenvalues.items()},
        "cases": records,
        "sensitivity": sensitivity,
        "display_reference_m": display_reference,
    }
    result_directory.mkdir(parents=True, exist_ok=True)
    (result_directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(
        result_directory / "opd_fields.npz", x_m=axis, y_m=axis,
        aperture_mask=mask, incidence_angles_deg=incidence_angles_deg,
        opd_m=np.stack(opd_maps), profile_axis_m=profile_axis,
        selected_angle_deg=selected_angle, x_axis_opd_m=x_profile, y_axis_opd_m=y_profile,
    )
    print(f"Figure: {figure_path}")
    print(f"Results: {result_directory}")


if __name__ == "__main__":
    main()

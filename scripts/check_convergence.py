"""Vary one numerical setting at a time for temperature and normal-path OPD.

The finest case in each sweep is a numerical reference, not an exact solution.
This script does not compare against FEM or impose an accuracy threshold.
"""

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np

from thermal_opd.eigenmodes import radial_eigenvalues
from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.metrics import circular_aperture_mask
from thermal_opd.opd import steady_opd
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.projection import fourier_bessel_coefficient
from thermal_opd.steady import steady_temperature


parameters = PhysicalParameters()
incidence_angle_deg = 60.0
aperture_radius_ratio = 0.90

# Numerical settings belong to this experiment; all internal units are SI.
azimuthal_orders = (0, 2, 4, 6, 8, 10)
radial_mode_counts = (6, 8, 12, 16, 20)
projection_r_counts = (100, 200, 400, 800, 1200)
projection_phi_counts = (24, 32, 64, 128, 256)
opd_xy_counts = (80, 120, 160, 201, 301)
num_temperature_r = 161
num_temperature_phi = 180
common_num_opd_xy = 201

repository_root = Path(__file__).resolve().parents[1]
result_directory = repository_root / "results/convergence"
figure_path = repository_root / "outputs/figures/convergence.png"


@dataclass(frozen=True)
class Settings:
    M: int
    N: int
    projection_r: int
    projection_phi: int
    opd_xy: int


def aperture_coordinates(count):
    radius = aperture_radius_ratio * parameters.radius
    axis = np.linspace(-radius, radius, count)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    mask = circular_aperture_mask(x, y, radius)
    return np.hypot(x[mask], y[mask]), np.arctan2(y[mask], x[mask])


def relative_disk_l2(fields, reference, r):
    """Per-plane relative L2 change with disk area weight r dr dphi."""
    numerator = np.trapezoid(np.mean((fields - reference)**2, axis=2) * r, x=r, axis=1)
    denominator = np.trapezoid(np.mean(reference**2, axis=2) * r, x=r, axis=1)
    return np.sqrt(numerator / denominator)


def main():
    p = parameters
    baseline = Settings(
        max(azimuthal_orders), max(radial_mode_counts),
        max(projection_r_counts), max(projection_phi_counts), common_num_opd_xy,
    )
    sweeps = {
        "azimuthal_order": ("M", azimuthal_orders),
        "radial_modes": ("N", radial_mode_counts),
        "projection_radial": ("projection_r", projection_r_counts),
        "projection_angular": ("projection_phi", projection_phi_counts),
        "opd_sampling": ("opd_xy", opd_xy_counts),
    }
    if min(projection_phi_counts) <= 2 * baseline.M:
        raise ValueError("Every angular projection grid must contain more than 2*M samples")
    print(f"Temperature / normal-path OPD convergence at {incidence_angle_deg:g} deg", flush=True)
    print(f"Fixed settings before each single-parameter sweep: {asdict(baseline)}", flush=True)
    roots = {
        m: radial_eigenvalues(m, p.radius, p.conductivity, p.h_eff, baseline.N)
        for m in range(0, baseline.M + 1, 2)
    }
    r = np.linspace(0.0, p.radius, num_temperature_r)
    phi = np.linspace(0.0, 2 * np.pi, num_temperature_phi, endpoint=False)
    z_planes = (0.0, p.thickness / 2, p.thickness)
    thermal_args = (p.conductivity, p.thickness, p.h_eff, p.h_eff)
    opd_args = (*thermal_args, p.dn_dT)
    common_r, common_phi = aperture_coordinates(common_num_opd_xy)

    # OPD is linear in source coefficients. Cache unit responses on the common
    # evaluation grid, reusing the production solver rather than a new formula.
    print("Preparing reusable OPD modal responses on the common aperture grid", flush=True)
    opd_basis = {
        m: np.stack([
            steady_opd(common_r, common_phi, {m: np.array([1.0])}, {m: np.array([lam])}, *opd_args)
            for lam in values
        ])
        for m, values in roots.items()
    }
    projection_cache = {}
    solution_cache = {}
    evaluation_cache = {}

    def solve(settings):
        key = (settings.M, settings.N, settings.projection_r, settings.projection_phi)
        if key in solution_cache:
            return solution_cache[key]
        grid_key = (settings.projection_r, settings.projection_phi)
        if grid_key not in projection_cache:
            projection_r = np.linspace(0.0, p.radius, settings.projection_r)
            projection_phi = np.linspace(0.0, 2 * np.pi, settings.projection_phi, endpoint=False)
            q = oblique_gaussian_heat_flux(
                projection_r[:, None] * np.cos(projection_phi),
                projection_r[:, None] * np.sin(projection_phi),
                p.beam_radius, np.deg2rad(incidence_angle_deg), p.absorbed_power,
            )
            projection_cache[grid_key] = {
                m: np.array([fourier_bessel_coefficient(q, projection_r, projection_phi, m, lam)
                             for lam in values])
                for m, values in roots.items()
            }
        selected_roots = {m: values[:settings.N] for m, values in roots.items() if m <= settings.M}
        coefficients = {m: projection_cache[grid_key][m][:settings.N] for m in selected_roots}
        fields = np.stack([
            steady_temperature(r[:, None], phi[None, :], z, coefficients, selected_roots, *thermal_args)
            for z in z_planes
        ])
        solution_cache[key] = coefficients, selected_roots, fields
        return solution_cache[key]

    def evaluate(settings):
        if settings in evaluation_cache:
            return evaluation_cache[settings]
        coefficients, selected_roots, fields = solve(settings)
        if settings.opd_xy == common_num_opd_xy:
            opd = np.zeros_like(common_r)
            for m, values in coefficients.items():
                opd += np.sum(values[:, None] * opd_basis[m][:settings.N], axis=0)
        else:
            evaluation_r, evaluation_phi = aperture_coordinates(settings.opd_xy)
            opd = steady_opd(evaluation_r, evaluation_phi, coefficients, selected_roots, *opd_args)
        record = {
            **asdict(settings),
            "opd_pv_m": float(np.ptp(opd)),
            "aperture_sample_count": int(opd.size),
            "center_temperature_rises_K": fields[:, 0, 0].tolist(),
            "sampled_peak_temperature_rises_K": fields.max(axis=(1, 2)).tolist(),
        }
        evaluation_cache[settings] = record, fields
        return record, fields

    # Check cache summation against an ordinary production call at several points.
    check_coefficients, check_roots, _ = solve(baseline)
    indices = np.linspace(0, common_r.size - 1, 9, dtype=int)
    expected = steady_opd(common_r[indices], common_phi[indices], check_coefficients, check_roots, *opd_args)
    cached = sum(np.sum(values[:, None] * opd_basis[m][:baseline.N, indices], axis=0)
                 for m, values in check_coefficients.items())
    np.testing.assert_allclose(cached, expected, rtol=1e-12, atol=1e-20)
    cache_max_difference = float(np.max(np.abs(cached - expected)))
    groups = {}
    csv_records = []
    for name, (setting_name, counts) in sweeps.items():
        reference_settings = replace(baseline, **{setting_name: counts[-1]})
        reference, reference_fields = evaluate(reference_settings)
        rows = []
        previous = None
        print(f"\n{name}: value, OPD-PV [nm], change vs finest, front-temperature L2 change", flush=True)
        for count in counts:
            record, fields = evaluate(replace(baseline, **{setting_name: count}))
            row = dict(record)
            row["varied_value"] = count
            row["relative_pv_change_vs_reference"] = abs(record["opd_pv_m"] - reference["opd_pv_m"]) / abs(reference["opd_pv_m"])
            row["temperature_relative_disk_l2_vs_reference"] = relative_disk_l2(fields, reference_fields, r).tolist()
            row["relative_pv_change_vs_previous"] = (
                None if previous is None else abs(record["opd_pv_m"] - previous) / abs(record["opd_pv_m"])
            )
            previous = record["opd_pv_m"]
            rows.append(row)
            print(f"{count:4d}  {row['opd_pv_m']*1e9:13.6f}  {row['relative_pv_change_vs_reference']:.6e}  {row['temperature_relative_disk_l2_vs_reference'][0]:.6e}", flush=True)
            flat_row = {k: v for k, v in row.items() if not isinstance(v, list)}
            flat_row["sweep"] = name
            for label, i in zip(("front", "mid", "rear"), range(3)):
                flat_row[f"{label}_center_temperature_K"] = row["center_temperature_rises_K"][i]
                flat_row[f"{label}_peak_temperature_K"] = row["sampled_peak_temperature_rises_K"][i]
                flat_row[f"{label}_relative_disk_l2_vs_reference"] = row["temperature_relative_disk_l2_vs_reference"][i]
            csv_records.append(flat_row)
        groups[name] = {"varied_setting": setting_name, "reference_settings": asdict(reference_settings), "cases": rows}

    summary = {
        "kind": "single_parameter_temperature_and_normal_path_opd_sweeps",
        "physical_parameters_SI": asdict(p),
        "h_eff_W_per_m2_K": p.h_eff,
        "incidence_angle_deg": incidence_angle_deg,
        "aperture_radius_ratio": aperture_radius_ratio,
        "opd_definition": "dn_dT * integral_0^L temperature_rise dz; normal path; unit m",
        "temperature_evaluation": {"r": num_temperature_r, "phi": num_temperature_phi, "z_m": z_planes},
        "reference_convention": "last (finest) case in each sweep; not an exact or FEM reference",
        "limits": "one angle; fixed temperature evaluation grid; sampling of PV extrema can mask changes in the rest of the OPD field; no acceptance threshold",
        "opd_cache_max_absolute_difference_m": cache_max_difference,
        "sweeps": groups,
    }
    result_directory.mkdir(parents=True, exist_ok=True)
    (result_directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (result_directory / "comparison.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_records[0]))
        writer.writeheader()
        writer.writerows(csv_records)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), layout="constrained")
    for ax, (name, group) in zip(axes.flat, groups.items()):
        rows = group["cases"]
        ax.plot([row["varied_value"] for row in rows], [row["opd_pv_m"]*1e9 for row in rows], "o-")
        ax.set_xlabel(group["varied_setting"])
        ax.set_ylabel("OPD-PV [nm]")
        ax.set_title(name.replace("_", " "))
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.grid(alpha=0.3)
    ax = axes.flat[-1]
    for name in ("azimuthal_order", "radial_modes"):
        rows = groups[name]["cases"][:-1]  # Reference change is zero, undefined on log axis.
        ax.semilogy([row["varied_value"] for row in rows],
                    [row["temperature_relative_disk_l2_vs_reference"][0] for row in rows],
                    "o-", label=groups[name]["varied_setting"])
    ax.set_xlabel("M or N")
    ax.set_ylabel("Front-temperature relative disk L2 change")
    ax.set_title("Change against each finest case")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.suptitle(f"Single-parameter sweeps at {incidence_angle_deg:g} deg; numerical references only")
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    print(f"\nFigure: {figure_path}\nResults: {result_directory}", flush=True)


if __name__ == "__main__":
    main()

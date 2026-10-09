"""Persist and independently reload transient Robin boundary diagnostics.

Sidecar files allow existing runs to gain boundary checks without overwriting
their fields, original parameters or original calculation provenance.
"""
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from thermal_opd.run_io import data_directory

from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.transient import PowerHistory, TransientModel
from thermal_opd.validation import transient_boundary_fields, transient_residual_fields


BOUNDARY_FILES = ('boundary_checks.npz','boundary_summary.json','boundary_residuals.csv')
RESIDUAL_KEYS = ('front_modal_residual_W_per_m2','rear_residual_W_per_m2','side_residual_W_per_m2')
CSV_COLUMNS = ('time_s','front_modal_max_abs_W_per_m2','rear_max_abs_W_per_m2','side_max_abs_W_per_m2',
               'front_modal_over_full_power_peak_flux','rear_over_full_power_peak_flux','side_over_full_power_peak_flux')


def model_from_saved_transient(summary, fields):
    """Restore exactly the saved physical parameters, history, J, C and roots."""
    kind = summary['kind']
    if kind not in ('transient_fine_fields','transient_modal_fem_study'):
        raise ValueError('Boundary check needs a transient field/FEM run, not a selection run')
    history = summary['precision_scope']['history'] if kind=='transient_fine_fields' else summary['history']
    J = summary['selected_axial_modes'] if kind=='transient_fine_fields' else summary['axial_counts'][-1]
    coefficients = {int(m):fields[f'source_coefficients_m{m}_W_per_m2'] for m in fields['orders']}
    roots = {int(m):fields[f'radial_roots_m{m}_per_m'] for m in fields['orders']}
    return TransientModel(PhysicalParameters(**summary['physical_parameters_SI']),coefficients,roots,
                          num_axial_modes=J,history=PowerHistory(**history))


def boundary_metrics(fields):
    return np.stack([np.max(np.abs(fields[key]),axis=(1,2)) for key in RESIDUAL_KEYS],axis=1)


def write_transient_boundary_check(directory,model,times,r,phi,side_z,*,angle_deg,
                                   relative_tolerance,absolute_tolerance,provenance=None):
    """Write actual boundary values and CSV/JSON metrics; return PASS/FAIL report.

    Relative errors use the full-power prescribed Gaussian peak, fixed across
    time. The absolute tolerance handles zero source without dividing by zero.
    This does not test pointwise FEM fluxes or source-projection convergence.
    """
    directory = Path(directory)
    if not np.isfinite(angle_deg) or abs(angle_deg)>=90:
        raise ValueError('Boundary incident angle must be finite and between -90 and 90 deg')
    if not all(np.isfinite(t) and t>=0 for t in (relative_tolerance,absolute_tolerance)):
        raise ValueError('Boundary tolerances must be finite and nonnegative')
    f = transient_boundary_fields(model,times,r,phi,side_z)
    p = model.parameters
    peak = float(oblique_gaussian_heat_flux(0.,0.,p.beam_radius,np.deg2rad(angle_deg),p.absorbed_power))
    scale = abs(peak)
    absolute = boundary_metrics(f)
    relative = absolute/scale if scale>0 else np.zeros_like(absolute)
    limit = absolute_tolerance+relative_tolerance*scale
    report = dict(version=1,status='PASS' if np.max(absolute)<=limit else 'FAIL',
        checked_at_utc=datetime.now(timezone.utc).isoformat(),physical_parameters_SI=asdict(p),
        history=asdict(model.history),num_axial_modes=len(model.axial.beta),angle_deg=float(angle_deg),
        source_target='g(t)*q_modal; projection error against the Gaussian is separate',
        derivative_method='analytic Fourier-Bessel and axial derivatives',
        equations=['h*T-k*dT/dz-g(t)*q_modal','h*T+k*dT/dz','h*T+k*dT/dr'],
        boundary_order=['front_modal','rear','side'],times_s=f['times_s'].tolist(),
        grid=dict(r=len(r),phi=len(phi),side_z=len(side_z)),
        full_power_peak_gaussian_flux_W_per_m2=scale,
        relative_tolerance=float(relative_tolerance),absolute_tolerance_W_per_m2=float(absolute_tolerance),
        allowed_absolute_residual_W_per_m2=float(limit),
        max_absolute_residuals_W_per_m2=absolute.tolist(),
        max_residuals_over_full_power_peak_flux=relative.tolist(),
        zero_source_note='At zero prescribed power, normalized values are 0; the absolute residual gate still applies.',
        scope='Sampled semi-analytical Robin consistency only; not PDE/J convergence or FEM pointwise flux validation.',
        provenance=provenance)
    np.savez_compressed(data_directory(directory)/BOUNDARY_FILES[0],**f)
    with (data_directory(directory)/BOUNDARY_FILES[2]).open('w',newline='',encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(CSV_COLUMNS)
        writer.writerows(np.column_stack((f['times_s'],absolute,relative)))
    (directory/BOUNDARY_FILES[1]).write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report


def validate_transient_boundary_check(directory,summary):
    """Recompute stored residuals, reconstruct the model and check snapshot T."""
    directory = Path(directory)
    locations=(data_directory(directory)/BOUNDARY_FILES[0],directory/BOUNDARY_FILES[1],data_directory(directory)/BOUNDARY_FILES[2])
    present = [path.is_file() for path in locations]
    def require(condition,message):
        if not condition:
            raise ValueError(message)
    if not any(present) and not summary.get('boundary_validation_required',False):
        return dict(status='NOT_RUN',note='Legacy run has no saved transient boundary check.')
    require(all(present),'Missing transient boundary check files')
    report = json.loads((directory/BOUNDARY_FILES[1]).read_text(encoding='utf-8'))
    require(report['version']==1,'Unsupported transient boundary check version')
    with np.load(data_directory(directory)/'fields.npz',allow_pickle=False) as saved, np.load(locations[0],allow_pickle=False) as f:
        model = model_from_saved_transient(summary,saved)
        p = model.parameters
        require(report['physical_parameters_SI']==asdict(p) and report['history']==asdict(model.history)
                and report['num_axial_modes']==len(model.axial.beta) and report['angle_deg']==summary['angle_deg'],
                'Transient boundary model differs from saved run')
        times = saved['times_s']
        require(np.array_equal(times,summary['times_s']) and np.array_equal(f['times_s'],times)
                and np.array_equal(report['times_s'],times),'Transient boundary times differ')
        fine = summary['kind']=='transient_fine_fields'
        r,phi = (saved[key] for key in (('temperature_r_m','temperature_phi_rad') if fine else ('r_m','phi_rad')))
        require(np.array_equal(f['r_m'],r) and np.array_equal(f['phi_rad'],phi),'Transient boundary grid differs from saved fields')
        require(report['grid']==dict(r=r.size,phi=phi.size,side_z=f['side_z_m'].size),'Transient boundary grid metadata differs')
        expected = transient_boundary_fields(model,times,r,phi,f['side_z_m'])
        require(set(f.files)==set(expected),'Transient boundary array keys differ')
        for key,value in expected.items():
            require(f[key].shape==value.shape and np.all(np.isfinite(f[key]))
                    and np.allclose(f[key],value,rtol=1e-11,atol=1e-9),f'Transient boundary array differs: {key}')
        T = saved['temperature_K'] if fine else saved['temperature_sa_K']
        for key,plane in (('front_temperature_K',0),('rear_temperature_K',2)):
            require(np.allclose(f[key],T[:,plane],rtol=1e-11,atol=1e-12),'Transient boundary temperature differs from actual saved snapshot')
        residuals = transient_residual_fields(f,p)
        for key,value in residuals.items():
            require(np.allclose(f[key],value,rtol=1e-11,atol=1e-12),'Saved transient boundary residual differs from T/gradient/load')
        absolute = boundary_metrics(residuals)
        peak = abs(float(oblique_gaussian_heat_flux(0.,0.,p.beam_radius,np.deg2rad(summary['angle_deg']),p.absorbed_power)))
        relative = absolute/peak if peak>0 else np.zeros_like(absolute)
        rel_tol,abs_tol = report['relative_tolerance'],report['absolute_tolerance_W_per_m2']
        require(all(np.isfinite(t) and t>=0 for t in (rel_tol,abs_tol)),'Invalid transient boundary tolerance')
        limit = abs_tol+rel_tol*peak
        require(np.isclose(report['full_power_peak_gaussian_flux_W_per_m2'],peak,rtol=1e-12,atol=0.)
                and np.isclose(report['allowed_absolute_residual_W_per_m2'],limit,rtol=1e-12,atol=0.),
                'Transient boundary normalization differs')
        for key,value in (('max_absolute_residuals_W_per_m2',absolute),('max_residuals_over_full_power_peak_flux',relative)):
            require(np.asarray(report[key]).shape==value.shape and np.allclose(report[key],value,rtol=1e-11,atol=1e-14),
                    'Transient boundary summary metrics differ')
        with locations[2].open(newline='',encoding='utf-8') as stream:
            rows = list(csv.reader(stream))
        require(tuple(rows[0])==CSV_COLUMNS,'Transient boundary CSV columns differ')
        table = np.asarray(rows[1:],dtype=float)
        expected_table = np.column_stack((times,absolute,relative))
        require(table.shape==expected_table.shape and np.allclose(table,expected_table,rtol=1e-11,atol=1e-14),
                'Transient boundary CSV differs')
        require(report['status']=='PASS' and np.max(absolute)<=limit,'Transient boundary residual exceeds tolerance')
    return dict(status='PASS',max_absolute_residuals_W_per_m2=np.max(absolute,axis=0).tolist(),
                boundary_order=report['boundary_order'],allowed_absolute_residual_W_per_m2=limit,
                report_file=BOUNDARY_FILES[1],table_file=str(locations[2].relative_to(directory)))

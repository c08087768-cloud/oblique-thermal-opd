"""Steady physical constraints, independent checks and saved experiment review."""
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv
from thermal_opd.projection import reconstruct_cosine_field
from thermal_opd.steady import steady_temperature,steady_temperature_gradient
from thermal_opd.run_io import data_directory,figure_directory,read_json,read_numeric_data,validation_report
from thermal_opd.convergence import field_errors,within_tolerance
from thermal_opd.validation_common import (require,validate_field_errors,mesh_boundary_edges,points_inside_polygon,
    finite_difference_laplacian,saved_validation_checks,saved_validation_result)

STEADY_RUN_KINDS=frozenset(('steady_joint_truncation_selection','steady_prismatic_fem_comparison',
                          'complete_representative_steady_refracted_ray_case'))

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


def steady_pde_residual(temperature,p,r,phi,z,*,space_step):
    """Independent steady cylindrical conduction residual, W/m³."""
    laplacian=finite_difference_laplacian(temperature,p,r,phi,z,space_step=space_step)
    return np.zeros_like(laplacian)-p.conductivity*laplacian


def validate_steady_selection(path, summary):
    with read_numeric_data(data_directory(path)/'fields.npz') as fields:
        require(all(np.all(np.isfinite(fields[key])) for key in fields.files), 'Nonfinite selection field')
        weights = fields['temperature_area_weights']
        require(len(summary['cases']) == len(summary['angles_deg']), 'Missing selection angle')
        for case in summary['cases']:
            label = f"theta{case['angle_deg']:g}"
            require(case['reference_pass'] and within_tolerance(case['reference_stability_errors'], summary['reference_stability_tolerance']), 'Reference stability failed')
            selected = next(c for c in case['candidates'] if c['settings']==summary['selected_settings'])
            require(selected['accepted'] and within_tolerance(selected['errors'], summary['tolerance']), 'Selected candidate failed')
            require(case['denser_common_grid_pass'] and within_tolerance(case['denser_common_grid_errors'], summary['tolerance']), 'Dense verification failed')
            errors = field_errors(fields[label+'_selected_temperature_K'], fields[label+'_reference_temperature_K'], weights,
                                  fields[label+'_selected_opd_m'], fields[label+'_reference_opd_m'])
            validate_field_errors(errors, selected['errors'], skip_pv=True)
            fine = field_errors(fields[label+'_selected_temperature_K'], fields[label+'_reference_temperature_K'], weights,
                                fields[label+'_verification_selected_opd_m'], fields[label+'_verification_reference_opd_m'])
            validate_field_errors(fine, case['denser_common_grid_errors'], skip_pv=True)
        require(summary['selected_figure_xy'] is not None, 'No accepted figure grid')
        for group in summary['sampling_checks']:
            row = next(r for r in group['cases'] if r['xy_count']==summary['selected_figure_xy'])
            require(row['relative_pv_sampling_error'] <= summary['tolerance'], 'Figure sampling failed')


def validate_steady_fem(path, summary):
    p, h = summary['physical_parameters_SI'], summary['h_eff_W_per_m2_K']
    require(summary['sanity_pass'], 'FEM sanity failed')
    require(len(summary['cases']) == len(summary['angles_deg']), 'Missing FEM angle')
    for case in summary['cases']:
        label = f"theta{case['angle_deg']:g}"
        with read_numeric_data(data_directory(path)/f'{label}_fields.npz') as fields:
            tm, om = fields['temperature_mask'], fields['opd_common_mask']
            require(tm.dtype == bool and om.dtype == bool and tm.ndim==2 and om.ndim==2, 'Invalid FEM mask format')
            require(int(tm.sum())==case['temperature_common_sample_count'], 'Temperature mask count mismatch')
            require(int(om.sum())==case['opd_common_sample_count'], 'OPD mask count mismatch')
            require(np.all(~om | fields['opd_physical_mask']), 'Common aperture exceeds physical mask')
            ts, os = fields['temperature_sa_K'][..., tm], fields['opd_sa_m'][om]
            tf, of = fields['temperature_fem_K'][..., tm], fields['opd_fem_m'][..., om]
            require(np.all(np.isnan(fields['temperature_fem_K'][..., ~tm])), 'Unexpected exterior FEM temperatures')
            require(np.all(np.isnan(fields['opd_fem_m'][..., ~om])), 'Unexpected exterior FEM OPD')
            require(np.allclose(np.ptp(os), case['sa_opd_pv_common_aperture_m'], rtol=1e-10, atol=1e-20), 'SA PV mismatch')
            require(tf.shape[0]==len(case['meshes']), 'Missing FEM mesh field')
            for j, row in enumerate(case['meshes']):
                actual = field_errors(tf[j], ts, np.ones(ts.shape[1]), of[j], os)
                validate_field_errors(actual, row['errors_vs_sa'])
                changes = field_errors(tf[j], tf[-1], np.ones(ts.shape[1]), of[j], of[-1])
                validate_field_errors(changes, row['changes_vs_finest_fem'])
                require(np.allclose(np.ptp(of[j]), row['fem_opd_pv_common_aperture_m'], rtol=1e-10, atol=1e-20), 'FEM PV mismatch')
                require(row['relative_linear_residual'] <= summary['sanity_tolerance'], 'FEM linear residual failed')
                with read_numeric_data(data_directory(path)/f'{label}_mesh{j}.npz') as mesh:
                    xy, tri, z, t = (mesh[k] for k in ('disk_nodes_xy_m', 'disk_triangles', 'axial_nodes_z_m', 'temperature_K'))
                    require(t.shape==(z.size, xy.shape[1]) and np.all(np.isfinite(t)), 'Bad FEM nodal temperature')
                    require(np.all(np.diff(z)>0) and z[0]==0 and z[-1]==p['thickness'], 'Bad FEM axial coordinates')
                    require(tri.min()>=0 and tri.max()<xy.shape[1], 'Bad FEM topology')
                    edges = mesh_boundary_edges(tri)
                    v1, v2 = xy[:, tri[1]]-xy[:, tri[0]], xy[:, tri[2]]-xy[:, tri[0]]
                    areas = np.abs(v1[0]*v2[1]-v1[1]*v2[0])/2
                    lengths = np.linalg.norm(xy[:, edges[:, 1]]-xy[:, edges[:, 0]], axis=0)
                    losses = (h*np.sum(areas*t[0, tri].mean(axis=0)),
                              h*np.sum(areas*t[-1, tri].mean(axis=0)),
                              h*np.trapezoid(np.sum(t[:, edges].mean(axis=2)*lengths, axis=1), x=z))
                    for name, value in zip(('front_loss_W', 'rear_loss_W', 'side_loss_W'), losses):
                        require(np.isclose(value, row[name], rtol=1e-10, atol=1e-14), f'Nodal heat loss mismatch: {name}')
                    require(abs(sum(losses)-row['input_W'])/abs(row['input_W']) <= summary['sanity_tolerance'], 'Reloaded FEM energy failed')
                    tx, ty = fields['temperature_x_m'][tm], fields['temperature_y_m'][tm]
                    ox, oy = fields['opd_x_m'][om], fields['opd_y_m'][om]
                    shift = p['thickness']*np.tan(np.deg2rad(case['refraction_angle_deg']))
                    require(np.all(points_inside_polygon(tx, ty, xy, edges)), 'Temperature point outside FEM polygon')
                    require(np.all(points_inside_polygon(ox, oy, xy, edges) & points_inside_polygon(ox, oy+shift, xy, edges)), 'Ray outside FEM polygon')
            q_check = case['load_quadrature_check']
            actual_q = field_errors(fields['temperature_higher_quadrature_K'][..., tm], tf[-1], np.ones(ts.shape[1]),
                                    fields['opd_higher_quadrature_m'][om], of[-1])
            validate_field_errors(actual_q, q_check['errors'])
            require(q_check['passed'] and within_tolerance(q_check['errors'], summary['mesh_change_target']/10), 'FEM load quadrature failed')
            expected_change = within_tolerance(case['meshes'][-2]['changes_vs_finest_fem'], summary['mesh_change_target'])
            require(case['finest_mesh_change_within_target']==expected_change, 'Incorrect FEM convergence flag')
        require((figure_directory(path)/f'{label}_comparison.png').is_file(), 'Missing FEM comparison figure')
    require((data_directory(path)/'comparison.csv').is_file(), 'Missing FEM comparison table')
    established = summary['mode']=='full' and all(c['finest_mesh_change_within_target'] for c in summary['cases'])
    require(summary['mesh_convergence_established']==established, 'Incorrect full mesh-convergence flag')


def validate_steady_fields(path,summary):
    """Check representative steady runs, including their legacy flat layout."""
    from thermal_opd.model import PreparedModel,ModalResolution
    from thermal_opd.parameters import PhysicalParameters
    from thermal_opd.opd import transmitted_ray_mask
    p=PhysicalParameters(**summary['physical_parameters_SI'])
    s=summary['numerical_settings']
    c={int(m):np.asarray(v) for m,v in summary['source_coefficients_W_per_m2'].items()}
    roots={int(m):np.asarray(v) for m,v in summary['eigenvalues_per_m'].items()}
    prepared=PreparedModel(p,summary['incidence_angle_deg'],ModalResolution(s['M'],s['N'],s['projection_r'],s['projection_phi'],s['opd_ray_points']),c,roots)
    solution=prepared.solve_steady()
    with read_numeric_data(data_directory(path)/'fields.npz') as f:
        r,phi,z=f['temperature_r_m'],f['temperature_phi_rad'],f['temperature_z_m']
        require(np.array_equal(z,[0.,p.thickness/2,p.thickness]),'Bad steady axial planes')
        expected=np.stack([solution.temperature(r[:,None],phi[None,:],zz) for zz in z])
        require(f['temperature_rise_K'].shape==expected.shape and np.all(np.isfinite(f['temperature_rise_K']))
                and np.allclose(f['temperature_rise_K'],expected,rtol=1e-11,atol=1e-13),'Steady temperature differs from saved modal model')
        x,y=f['opd_x_m'],f['opd_y_m']
        if x.ndim==y.ndim==1:
            x,y=np.meshgrid(x,y)
        mask=f['aperture_mask']
        require(mask.dtype==bool and mask.shape==x.shape==y.shape==f['opd_m'].shape,'Bad steady optical shape')
        require(np.all(np.isnan(f['opd_m'][~mask])),'Unexpected steady exterior OPD')
        physical=transmitted_ray_mask(x,y,p.radius,p.thickness,prepared.refraction_angle)
        entrance=x*x+y*y<=(p.radius*summary['aperture_radius_ratio'])**2
        require(np.array_equal(mask,physical&entrance),'Steady aperture mask differs from geometry')
        optical=solution.opd(x[mask],y[mask],point_chunk_size=s['opd_point_chunk_size'])
        require(np.allclose(f['opd_m'][mask],optical,rtol=1e-11,atol=1e-20),'Steady OPD differs from saved model')
        require(np.isclose(np.ptp(optical),summary['opd_pv_m'],rtol=1e-11,atol=1e-20),'Steady PV differs')
        bc=steady_boundary_residuals(r,phi,np.linspace(0,p.thickness,s['boundary_z']),c,roots,p.radius,p.conductivity,p.thickness,p.h_eff,p.h_eff,p.h_eff)
        require(max(bc.values())<1e-7,'Steady boundary residual above tolerance')
        for key,value in bc.items():
            require(np.isclose(value,summary['boundary_absolute_residuals_W_per_m2'][key],atol=1e-10,rtol=1e-8),'Steady boundary metrics differ')
        energy=steady_energy_balance(c,roots,p.radius,p.conductivity,p.thickness,p.h_eff,p.h_eff,p.h_eff,
                                     num_r=s['energy_r'],num_z=s['energy_z'],num_phi=s['energy_phi'])
        for key,value in energy.items():
            require(np.isclose(value,summary['energy_balance'][key],rtol=1e-10,atol=1e-12),'Steady energy metrics differ')
        require(energy['balance_absolute_error_W']<=1e-9+abs(energy['modal_input_W'])*1e-7,'Steady energy balance above tolerance')


def steady_saved_validation_report(directory,summary):
    """Record only checks actually completed by steady saved-result review."""
    checks=saved_validation_checks(directory,summary)
    kind=summary.get('kind')
    if kind=='steady_joint_truncation_selection':
        checks['joint_modal_truncation']['status']='passed'
    if kind=='complete_representative_steady_refracted_ray_case':
        checks['thermal_boundaries'].update(status='passed',evidence='summary.json',
            metrics=summary.get('boundary_absolute_residuals_W_per_m2') or {},threshold=1e-7)
        checks['steady_energy_balance']['status']='passed'
    if kind=='steady_prismatic_fem_comparison':
        checks['fem_energy_balance']['status']='passed'
        checks['fem_weak_equation_residual']['status']='passed'
    if summary.get('mesh_convergence_established',False):
        checks['fem_space_convergence']['status']='passed'
    return validation_report(list(checks.values()),config={'physical_parameters_SI':summary.get('physical_parameters_SI'),
        'kind':kind,'scope':summary.get('precision_scope')},run_directory=directory)


def validate_saved_steady_run(directory,summary=None):
    """Reload steady fields and recompute numerical gates; accept original flat runs."""
    path=Path(directory)
    summary=read_json(path/'summary.json') if summary is None else dict(summary)
    kind=summary.get('kind')
    if kind=='complete_representative_steady_refracted_ray_case' and 'status' not in summary:
        summary['status']='COMPLETE'
    require(summary.get('status')=='COMPLETE','Run is not complete')
    if kind=='steady_joint_truncation_selection':
        validate_steady_selection(path,summary)
    elif kind=='steady_prismatic_fem_comparison':
        validate_steady_fem(path,summary)
    elif kind=='complete_representative_steady_refracted_ray_case':
        validate_steady_fields(path,summary)
    else:
        raise ValueError(f'Unsupported run kind: {kind}')
    return saved_validation_result(summary,steady_saved_validation_report(path,summary))

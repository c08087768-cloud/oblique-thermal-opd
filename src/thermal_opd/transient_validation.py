"""Transient constraints, independent references and saved experiment review.

Boundary sidecars retain their original schema and flat-run reading support.
"""
import csv
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from thermal_opd.run_io import data_directory

from thermal_opd.heat_source import oblique_gaussian_heat_flux
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.model import PreparedModel,ModalResolution
from thermal_opd.transient import PowerHistory
from scipy.linalg import eigh_tridiagonal
from thermal_opd.projection import reconstruct_cosine_field
from thermal_opd.run_io import figure_directory,resolve_run_path,read_json,read_numeric_data,write_json,validation_report
from thermal_opd.convergence import field_errors,within_tolerance
from thermal_opd.validation_common import (require,validate_field_errors,mesh_boundary_edges,points_inside_polygon,
    finite_difference_laplacian,saved_validation_checks,saved_validation_result)

TRANSIENT_RUN_KINDS=frozenset(('transient_modal_fem_study','transient_joint_truncation_selection','transient_fine_fields'))


BOUNDARY_FILES = ('boundary_checks.npz','boundary_summary.json','boundary_residuals.csv')
RESIDUAL_KEYS = ('front_modal_residual_W_per_m2','rear_residual_W_per_m2','side_residual_W_per_m2')
CSV_COLUMNS = ('time_s','front_modal_max_abs_W_per_m2','rear_max_abs_W_per_m2','side_max_abs_W_per_m2',
               'front_modal_over_full_power_peak_flux','rear_over_full_power_peak_flux','side_over_full_power_peak_flux')


def solution_from_saved_transient(summary, fields):
    """Restore exactly the saved physical parameters, history, J, C and roots."""
    kind = summary['kind']
    if kind not in ('transient_fine_fields','transient_modal_fem_study'):
        raise ValueError('Boundary check needs a transient field/FEM run, not a selection run')
    history = summary['precision_scope']['history'] if kind=='transient_fine_fields' else summary['history']
    J = summary['selected_axial_modes'] if kind=='transient_fine_fields' else summary['axial_counts'][-1]
    coefficients = {int(m):fields[f'source_coefficients_m{m}_W_per_m2'] for m in fields['orders']}
    roots = {int(m):fields[f'radial_roots_m{m}_per_m'] for m in fields['orders']}
    prepared=PreparedModel(PhysicalParameters(**summary['physical_parameters_SI']),summary['angle_deg'],
                           ModalResolution(**summary['modal_settings']),coefficients,roots)
    return prepared.solve_transient(num_axial_modes=J,history=PowerHistory(**history))


def boundary_metrics(fields):
    return np.stack([np.max(np.abs(fields[key]),axis=(1,2)) for key in RESIDUAL_KEYS],axis=1)


def write_transient_boundary_check(directory,model,times,r,phi,side_z,*,
                                   relative_tolerance,absolute_tolerance,provenance=None):
    """Write actual boundary values and CSV/JSON metrics; return PASS/FAIL report.

    Relative errors use the full-power prescribed Gaussian peak, fixed across
    time. The absolute tolerance handles zero source without dividing by zero.
    This does not test pointwise FEM fluxes or source-projection convergence.
    """
    directory = Path(directory)
    angle_deg=model.prepared.angle_deg
    if not all(np.isfinite(t) and t>=0 for t in (relative_tolerance,absolute_tolerance)):
        raise ValueError('Boundary tolerances must be finite and nonnegative')
    f = transient_boundary_fields(model,times,r,phi,side_z)
    p = model.prepared.parameters
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
    write_json(directory/BOUNDARY_FILES[1],report)
    return report


def validate_transient_boundary_check(directory,summary):
    """Recompute stored residuals, reconstruct the model and check snapshot T."""
    directory = Path(directory)
    locations=(data_directory(directory)/BOUNDARY_FILES[0],directory/BOUNDARY_FILES[1],data_directory(directory)/BOUNDARY_FILES[2])
    present = [path.is_file() for path in locations]
    if not any(present) and not summary.get('boundary_validation_required',False):
        return dict(status='NOT_RUN',note='Legacy run has no saved transient boundary check.')
    require(all(present),'Missing transient boundary check files')
    report = read_json(directory/BOUNDARY_FILES[1])
    require(report['version']==1,'Unsupported transient boundary check version')
    with read_numeric_data(data_directory(directory)/'fields.npz') as saved, read_numeric_data(locations[0]) as f:
        model = solution_from_saved_transient(summary,saved)
        p = model.prepared.parameters
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


def transient_boundary_fields(model, times, r, phi, side_z):
    """Sample signed Robin residuals against g(t)*q_modal, in W/m².

    Boundary temperatures and analytic normal-coordinate derivatives are saved
    separately so validators can reconstruct residuals. This checks the retained
    model, not source projection accuracy or transient PDE/J convergence.
    At t=0 both T and the history's load are zero, including a step convention.
    """
    p = model.prepared.parameters
    times,r,phi,side_z = (np.asarray(a,dtype=float) for a in (times,r,phi,side_z))
    if any(a.ndim!=1 or not a.size or not np.all(np.isfinite(a)) for a in (times,r,phi,side_z)):
        raise ValueError('Require finite nonempty 1D boundary coordinates/times')
    if times[0]<0 or np.any(np.diff(times)<=0) or r[0]<0 or r[-1]>p.radius or np.any(np.diff(r)<=0):
        raise ValueError('Invalid boundary times/radii')
    if phi[0]!=0 or phi[-1]>=2*np.pi or np.any(np.diff(phi)<=0):
        raise ValueError('Invalid boundary azimuths')
    if side_z[0]!=0 or side_z[-1]!=p.thickness or np.any(np.diff(side_z)<=0):
        raise ValueError('Side boundary must cover z=0 to L')
    source = reconstruct_cosine_field(model.prepared.coefficients,model.prepared.roots,r,phi)
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


def transient_pde_residual(temperature,p,r,phi,z,time,*,space_step,time_step):
    """Independent rho cp T_t-k Laplacian(T) residual, W/m³."""
    r,phi,z=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(phi,dtype=float),np.asarray(z,dtype=float))
    h,dt=space_step,time_step
    if not all(np.isfinite(v) and v>0 for v in (h,dt)) or time<=dt:
        raise ValueError('Require positive finite steps and t>dt')
    laplacian=finite_difference_laplacian(lambda r,phi,z:temperature(r,phi,z,time),p,r,phi,z,space_step=h)
    temporal=(temperature(r,phi,z,time+dt)-temperature(r,phi,z,time-dt))/(2*dt)
    return p.density*p.specific_heat*temporal-p.conductivity*laplacian


def transient_uniform_flux_1d_reference(thickness,conductivity,volumetric_capacity,flux,h_front,h_back,times,*,num_intervals):
    """Independent 1D finite-volume diffusion with exact-in-time step response.

    For a cylinder with insulated side and spatially uniform front flux, T
    depends only on z. Half control volumes at the faces impose Robin losses.
    Symmetric tridiagonal diagonalization avoids the prism FEM assembly and
    the semi-analytical Fourier-Bessel/axial characteristic equations.
    This reference still requires its own spatial refinement.
    """
    from operator import index
    n=index(num_intervals)
    times=np.asarray(times,dtype=float)
    if n<2 or times.ndim!=1 or not np.all(np.isfinite(times)) or np.any(times<0):
        raise ValueError('Invalid 1D reference grid or times')
    if not all(np.isfinite(v) and v>0 for v in (thickness,conductivity,volumetric_capacity)) or min(h_front,h_back)<0 or h_front+h_back<=0:
        raise ValueError('Require positive material/geometry and dissipative Robin faces')
    z=np.linspace(0,thickness,n+1)
    dz=thickness/n
    mass=np.full(n+1,volumetric_capacity*dz);mass[[0,-1]]*=.5
    diagonal=np.full(n+1,2*conductivity/dz)
    diagonal[0]=conductivity/dz+h_front
    diagonal[-1]=conductivity/dz+h_back
    off=np.full(n,-conductivity/dz)
    rates,vectors=eigh_tridiagonal(diagonal/mass,off/np.sqrt(mass[:-1]*mass[1:]))
    forcing=flux*vectors[0]/np.sqrt(mass[0])
    response=-np.expm1(-times[:,None]*rates[None,:])/rates[None,:]*forcing[None,:]
    temperature=(response@vectors.T)/np.sqrt(mass)[None,:]
    return z,temperature


def validate_transient_fem(path, summary):
    """Recompute mode sums, field errors and nodal storage/loss from saved files."""
    p,h = summary['physical_parameters_SI'],summary['h_eff_W_per_m2_K']
    tolerance = summary['sanity_tolerance']
    require(summary['sanity_pass'],'Transient sanity failed')
    require(not summary['transient_truncation_established'] and not summary['fem_convergence_established'],
            'Diagnostic study cannot claim final convergence')
    with read_numeric_data(data_directory(path)/'fields.npz') as f:
        require(all(np.all(np.isfinite(f[key])) for key in f.files),'Nonfinite transient field')
        times,weights = f['times_s'],f['area_weights_m2']
        require(np.array_equal(times,summary['times_s']) and times[0]==0 and np.all(np.diff(times)>0),'Bad transient times')
        temperature,opd = f['temperature_sa_K'],f['opd_sa_m']
        require(np.all(temperature[0]==0) and np.all(opd[0]==0),'Nonzero initial SA field')
        for total,orders in (('temperature_sa_K','order_temperature_K'),('opd_sa_m','order_opd_m'),
                             ('steady_temperature_K','steady_order_temperature_K'),('steady_opd_m','steady_order_opd_m'),
                             ('source_projected_W_per_m2','source_order_W_per_m2')):
            require(np.allclose(f[total],f[orders].sum(axis=0),rtol=1e-11,atol=1e-18),'Thermal mode sum mismatch')
        energies = np.sum(f['order_temperature_K']**2*weights[None,None,None],axis=(-2,-1))
        shares = np.divide(energies,energies.sum(axis=0),out=np.zeros_like(energies),where=energies.sum(axis=0)!=0)
        require(np.allclose(shares,f['temperature_plane_l2_squared_shares'],rtol=1e-11,atol=1e-14),'Modal norm shares differ')
        require(np.allclose(np.ptp(f['order_opd_m'],axis=-1),f['order_sampled_opd_pv_m'],rtol=1e-11,atol=1e-20),'Modal optical PV differs')
        long_errors = field_errors(f['long_time_temperature_K'].reshape(3,-1),f['steady_temperature_K'].reshape(3,-1),
                                  weights.ravel(),f['long_time_opd_m'],f['steady_opd_m'])
        validate_field_errors(long_errors,summary['long_time_errors_vs_steady'])
        for j,t in enumerate(times):
            if t==0:
                continue
            for i in range(len(summary['axial_counts'])):
                err = field_errors(f['axial_check_temperature_K'][i,j].reshape(3,-1),temperature[j].reshape(3,-1),weights.ravel(),
                                   f['axial_check_opd_m'][i,j],opd[j])
                validate_field_errors(err,summary['axial_errors_vs_largest_J'][i][j])
            err = field_errors(temperature[j].reshape(3,-1),f['reference_temperature_K'][j].reshape(3,-1),weights.ravel(),
                               opd[j],f['reference_opd_m'][j])
            validate_field_errors(err,summary['transverse_errors_vs_finite_reference'][j])
            err = field_errors(temperature[j].reshape(3,-1),temperature[j].reshape(3,-1),weights.ravel(),
                               opd[j],f['higher_ray_quadrature_opd_m'][j])
            validate_field_errors(err,summary['ray_quadrature_errors'][j])
        mask,common,tm = f['opd_physical_mask'],f['opd_fem_common_mask'],f['temperature_fem_common_mask']
        require(mask.dtype==common.dtype==tm.dtype==bool,'Bad transient masks')
        require(np.all(~common|mask),'FEM optical aperture exceeds physical aperture')
        x,y = f['opd_x_m'],f['opd_y_m']
        angle = np.deg2rad(summary['refraction_angle_deg'])
        shift = p['thickness']*np.tan(angle)
        require(np.all(x[mask]**2+y[mask]**2<=(p['radius']*summary['aperture_radius_ratio'])**2*(1+1e-14)),
                'Ray entrance outside aperture')
        require(np.all(x[mask]**2+(y[mask]+shift)**2<=p['radius']**2*(1+1e-14)),'Ray exits cylinder')
        require(np.array_equal(f['opd_sa_fem_common_m'],opd[:,common[mask]]),'Common SA OPD differs')
        require(f['temperature_fem_K'].shape[0]==len(summary['fem_runs']),'Missing transient FEM run')
        for i,row in enumerate(summary['fem_runs']):
            require(row['max_relative_linear_residual']<tolerance and row['max_energy_balance_error_relative_to_full_input']<tolerance,
                    'Transient FEM solver/energy failed')
            tf,of = f['temperature_fem_K'][i],f['opd_fem_m'][i]
            difference = tf-temperature[...,tm]
            require(np.allclose(np.max(np.abs(difference),axis=-1),row['temperature_absolute_max_difference_K'],rtol=1e-10,atol=1e-13),
                    'Absolute temperature errors differ')
            rms = np.sqrt(np.sum(difference**2*weights[tm],axis=-1)/weights[tm].sum())
            require(np.allclose(rms,row['temperature_area_rms_difference_K'],rtol=1e-10,atol=1e-13),'Temperature RMS errors differ')
            require(np.all(tf[0]==0) and np.all(of[0]==0),'Nonzero initial FEM field')
            for j,t in enumerate(times):
                if t>0:
                    err = field_errors(tf[j],temperature[j][...,tm],weights[tm],of[j],f['opd_sa_fem_common_m'][j])
                    validate_field_errors(err,row['errors_vs_sa'][j])
            with read_numeric_data(data_directory(path)/f'fem_run{i}.npz') as mesh:
                xy,tri,z,T = (mesh[key] for key in ('disk_nodes_xy_m','disk_triangles','axial_nodes_z_m','temperature_K'))
                require(T.shape==(times.size,z.size,xy.shape[1]) and np.all(np.isfinite(T)),'Bad transient FEM nodal fields')
                require(np.array_equal(mesh['times_s'],times) and np.all(T[0]==0),'Bad FEM initial time')
                require(z[0]==0 and z[-1]==p['thickness'] and np.all(np.diff(z)>0),'Bad FEM layers')
                require(tri.min()>=0 and tri.max()<xy.shape[1],'Bad FEM mesh topology')
                v1,v2 = xy[:,tri[1]]-xy[:,tri[0]],xy[:,tri[2]]-xy[:,tri[0]]
                areas = np.abs(v1[0]*v2[1]-v1[1]*v2[0])/2
                edges = mesh_boundary_edges(tri)
                lengths = np.linalg.norm(xy[:,edges[:,1]]-xy[:,edges[:,0]],axis=0)
                stored = p['density']*p['specific_heat']*np.trapezoid(np.sum(T[:,:,tri].mean(axis=-2)*areas,axis=-1),x=z,axis=-1)
                require(np.all(points_inside_polygon(x[common],y[common],xy,edges)&points_inside_polygon(x[common],y[common]+shift,xy,edges)),
                        'Transient ray outside FEM polygon')
                steps = {key:mesh['steps'][:,k] for k,key in enumerate(mesh['step_keys'])}
                require(all(np.all(np.isfinite(v)) for v in steps.values()),'Nonfinite transient FEM step')
                require(np.all(steps['dt_s']>0) and steps['time_s'][-1]==times[-1],'Invalid time integration range')
                require(np.allclose(np.diff(np.r_[0.,steps['time_s']]),steps['dt_s'],rtol=1e-10,atol=1e-12),'Incorrect step times')
                actual_storage = np.diff(np.r_[0.,steps['stored_energy_J']])/steps['dt_s']
                require(np.allclose(actual_storage,steps['storage_rate_W'],rtol=1e-8,atol=1e-12),'Stored-energy rate mismatch')
                balance = abs(actual_storage+steps['heat_loss_W']-steps['input_W'])
                scale = max(abs(row['full_input_W']),np.finfo(float).tiny)
                require(np.max(balance)/scale<tolerance,'Saved transient energy failed')
                from thermal_opd.transient import PowerHistory
                history = PowerHistory(**summary['history'])
                expected_input = row['full_input_W']*np.array([history.value(t) for t in steps['time_s']])
                require(np.allclose(steps['input_W'],expected_input,rtol=1e-12,atol=1e-14),'Time-dependent load differs')
                for j,t in enumerate(times[1:],1):
                    k = np.flatnonzero(np.isclose(steps['time_s'],t,rtol=0,atol=1e-12))
                    require(k.size==1,'FEM missed observation time')
                    require(np.isclose(stored[j],steps['stored_energy_J'][k[0]],rtol=1e-10,atol=1e-12),'Nodal stored energy mismatch')
                    loss = h*(np.sum(areas*T[j,0,tri].mean(axis=0))+np.sum(areas*T[j,-1,tri].mean(axis=0))
                              +np.trapezoid(np.sum(T[j,:,edges].mean(axis=1)*lengths[:,None],axis=0),x=z))
                    require(np.isclose(loss,steps['heat_loss_W'][k[0]],rtol=1e-9,atol=1e-12),'Nodal heat loss mismatch')
        require((figure_directory(path)/'transient_modes.png').is_file() and (data_directory(path)/'modal_contributions.csv').is_file(),'Missing modal study artifacts')


def validate_transient_selection(path,summary):
    require(len(summary['cases'])==len(summary['angles_deg']),'Missing transient selection angle')
    with read_numeric_data(data_directory(path)/'fields.npz') as f:
        require(all(np.all(np.isfinite(f[key])) for key in f.files),'Nonfinite transient selection field')
        weights = f['temperature_area_weights']
        for case in summary['cases']:
            label = f'theta{case["angle_deg"]:g}'
            selected = next(row for row in case['candidates'] if row['settings']==summary['selected_settings'])
            require(case['reference_pass'] and selected['accepted'] and case['denser_common_grid_pass'] and case['sampling_pass'],
                    'Transient selection gate failed')
            for j in range(len(summary['times_s'])):
                for left,right,pv_left,pv_right,stored,limit in (
                    ('reference','enlarged','reference','enlarged',case['reference_errors_by_time'],summary['reference_stability_tolerance']),
                    ('selected','reference','selected','reference',selected['errors_by_time'],summary['tolerance'])):
                    err = field_errors(f[label+'_'+left+'_temperature_K'][j],f[label+'_'+right+'_temperature_K'][j],weights,
                                       f[label+'_'+left+'_opd_m'][j],f[label+'_'+right+'_opd_m'][j])
                    a,b = f[label+'_'+pv_left+'_continuous_pv_m'][j],f[label+'_'+pv_right+'_continuous_pv_m'][j]
                    require(a>0 and b>0,'Invalid continuous transient PV')
                    err['opd_relative_pv'] = abs(a-b)/b
                    validate_field_errors(err,stored[j])
                    require(within_tolerance(err,limit),'Transient reference/candidate above target')
                err = field_errors(f[label+'_selected_temperature_K'][j],f[label+'_reference_temperature_K'][j],weights,
                                   f[label+'_fine_selected_opd_m'][j],f[label+'_fine_reference_opd_m'][j])
                err['opd_relative_pv'] = abs(f[label+'_selected_continuous_pv_m'][j]-f[label+'_reference_continuous_pv_m'][j])/f[label+'_reference_continuous_pv_m'][j]
                validate_field_errors(err,case['denser_errors_by_time'][j])
                require(within_tolerance(err,summary['tolerance']),'Dense transient comparison above target')
                pv = f[label+'_selected_continuous_pv_m'][j]
                sampling = abs(np.ptp(f[label+'_figure_opd_m'][j])-pv)/pv
                require(np.isclose(sampling,case['sampling_errors_by_time'][j],rtol=1e-9,atol=1e-12),'Transient sampling metric differs')
                require(sampling<=summary['tolerance'],'Transient figure sampling above target')


def validate_transient_fields(path,summary):
    p = summary['physical_parameters_SI']
    if summary['transient_truncation_established']:
        selected_path = resolve_run_path(path,summary['selection_directory'])
        validate_saved_transient_run(selected_path)
        selected = read_json(selected_path/'summary.json')
        require(selected['kind']=='transient_joint_truncation_selection' and selected['physical_parameters_SI']==p,
                'Fine transient selection has different physics')
        require(summary['angle_deg'] in selected['angles_deg'] and summary['modal_settings']==selected['selected_settings']['modal']
                and summary['selected_axial_modes']==selected['selected_settings']['J'],'Fine transient outside selected scope')
        require(summary['precision_scope']['times_s']==selected['times_s'],'Selected times differ')
    require(not summary['fem_convergence_established'],'Fine SA fields do not establish FEM convergence')
    with read_numeric_data(data_directory(path)/'fields.npz') as f:
        times,r,phi,z = (f[key] for key in ('times_s','temperature_r_m','temperature_phi_rad','temperature_z_m'))
        T,O,mask = f['temperature_K'],f['opd_m'],f['aperture_mask']
        s = summary['numerical_settings']
        require(np.array_equal(times,summary['times_s']) and times[0]==0 and np.all(np.diff(times)>0),'Invalid fine transient times')
        require(T.shape==(times.size,3,s['temperature_r'],s['temperature_phi']) and np.all(np.isfinite(T)),'Bad fine temperature shape/values')
        require(np.all(np.isfinite(r)) and np.all(np.diff(r)>0) and r[0]==0 and r[-1]<=p['radius'],'Bad fine radii')
        require(np.all(np.isfinite(phi)) and np.all(np.diff(phi)>0) and phi[0]==0 and phi[-1]<2*np.pi,'Bad fine angles')
        require(np.array_equal(z,[0.,p['thickness']/2,p['thickness']]),'Bad fine axial planes')
        require(mask.dtype==bool and mask.shape==(s['opd_xy'],s['opd_xy']) and O.shape==(times.size,)+mask.shape,'Bad fine OPD shape')
        require(np.all(np.isnan(O[:,~mask])),'Unexpected exterior fine OPD values')
        require(np.all(np.isfinite(O[:,mask])) and np.all(T[0]==0) and np.all(O[0,mask]==0),'Bad fine initial/optical values')
        x,y = f['opd_x_m'],f['opd_y_m']
        shift = p['thickness']*np.tan(np.deg2rad(summary['refraction_angle_deg']))
        require(np.all(x[mask]**2+y[mask]**2<=(p['radius']*summary['aperture_radius_ratio'])**2*(1+1e-14)), 'Fine entry outside aperture')
        require(np.all(x[mask]**2+(y[mask]+shift)**2<=p['radius']**2*(1+1e-14)), 'Fine ray leaves cylinder')
        require(np.allclose(T,f['order_temperature_K'].sum(axis=0),rtol=1e-11,atol=1e-14),'Fine thermal mode sum differs')
        require(np.allclose(O[:,mask],f['order_opd_m'][...,mask].sum(axis=0),rtol=1e-11,atol=1e-20),'Fine optical mode sum differs')
        sampled = np.ptp(O[:,mask],axis=1)
        require(np.allclose(sampled,f['sampled_pv_m'],rtol=1e-11,atol=1e-20),'Fine sampled PV differs')
        require(len(summary['snapshots'])==times.size,'Missing fine snapshot')
        for j,row in enumerate(summary['snapshots']):
            require(row['time_s']==times[j] and np.isclose(sampled[j],row['sampled_pv_m'],rtol=1e-11,atol=1e-20),'Saved snapshot PV differs')
            if summary['transient_truncation_established'] and times[j] in summary['precision_scope']['times_s']:
                require(row['continuous_pv_m']>0 and abs(sampled[j]-row['continuous_pv_m'])/row['continuous_pv_m']<=summary['precision_target'],
                        'Fine snapshot PV sampling above target')
        if summary['transient_truncation_established']:
            model = solution_from_saved_transient(summary,f)
            ii,jj = np.array([0,r.size//2,r.size-1]),np.array([0,phi.size//4,phi.size//2])
            for j in (1,times.size-2,times.size-1):
                for plane,depth in enumerate(z):
                    expected = model.temperature(r[ii],phi[jj],depth,times[j])
                    require(np.allclose(T[j,plane,ii,jj],expected,rtol=1e-11,atol=1e-13),'Fine field does not match saved model')
            w = np.broadcast_to(np.maximum(r[:,None],p['radius']*1e-9),T.shape[-2:]).ravel()
            errors = field_errors(T[-1].reshape(3,-1),f['steady_temperature_K'].reshape(3,-1),w,O[-1,mask],f['steady_opd_m'][mask])
            validate_field_errors(errors,summary['long_time_errors_vs_same_modes_steady'])
            require(within_tolerance(errors,summary['precision_target']),'Fine long-time limit above target')
            with read_numeric_data(data_directory(path)/'history.npz') as history:
                require(all(np.all(np.isfinite(history[key])) for key in history.files),'Nonfinite fine history')
                ht,ho,hT = history['times_s'],history['opd_m'],history['temperature_K']
                require(ht[0]==0 and np.all(np.diff(ht)>0) and ht[-1]==times[-1],'Invalid fine history times')
                require(np.allclose(np.ptp(ho,axis=1),history['sampled_pv_m'],rtol=1e-11,atol=1e-20),'Fine history PV differs')
                for j,time in enumerate(times):
                    k = np.flatnonzero(ht==time)
                    require(k.size==1,'Fine history missed snapshot')
                    require(np.allclose(hT[k[0],:,0],T[j,:,0,0],rtol=1e-11,atol=1e-13),'Fine center history differs')
                    require(np.allclose(ho[k[0]],O[j,mask],rtol=1e-11,atol=1e-20),'Fine history optical field differs')
    figure_path = resolve_run_path(path,summary['figure_directory']) if summary.get('figure_directory') else figure_directory(path)
    for filename in summary['figures']:
        require((figure_path/filename).is_file(),'Missing fine transient figure')


def transient_saved_validation_report(directory,summary,boundary_check=None):
    """Record initial, time, long-time and FEM evidence in the original report schema."""
    checks=saved_validation_checks(directory,summary)
    kind=summary.get('kind')
    if kind=='transient_joint_truncation_selection' or summary.get('transient_truncation_established',False):
        checks['joint_modal_truncation']['status']='passed'
    if kind in ('transient_fine_fields','transient_modal_fem_study'):
        checks['initial_condition'].update(status='passed',metrics={'max_initial_temperature_K':0.,'max_initial_opd_m':0.},threshold=0.)
    long_errors=summary.get('long_time_errors_vs_same_modes_steady',summary.get('long_time_errors_vs_steady'))
    if long_errors is not None:
        checks['long_time_steady_limit'].update(metrics=long_errors,
            status='passed' if within_tolerance(long_errors,checks['long_time_steady_limit']['threshold']) else 'failed')
    checks['thermal_boundaries'].update(metrics=boundary_check or {},
        threshold=(boundary_check or {}).get('allowed_absolute_residual_W_per_m2'))
    if boundary_check is not None and boundary_check['status']=='PASS':
        checks['thermal_boundaries']['status']='passed'
    if kind=='transient_modal_fem_study':
        checks['fem_energy_balance']['status']='passed'
        checks['fem_weak_equation_residual']['status']='passed'
    if summary.get('mesh_convergence_established',False):
        checks['fem_space_convergence']['status']='passed'
    return validation_report(list(checks.values()),config={'physical_parameters_SI':summary.get('physical_parameters_SI'),
        'kind':kind,'scope':summary.get('precision_scope')},run_directory=directory)


def validate_saved_transient_run(directory,summary=None):
    """Reload transient fields, recompute gates and independently inspect boundaries."""
    path=Path(directory)
    summary=read_json(path/'summary.json') if summary is None else summary
    kind=summary.get('kind')
    require(summary.get('status')=='COMPLETE','Run is not complete')
    if kind=='transient_modal_fem_study':
        validate_transient_fem(path,summary)
    elif kind=='transient_joint_truncation_selection':
        validate_transient_selection(path,summary)
    elif kind=='transient_fine_fields':
        validate_transient_fields(path,summary)
    else:
        raise ValueError(f'Unsupported run kind: {kind}')
    boundary_check=None
    if kind in ('transient_modal_fem_study','transient_fine_fields'):
        boundary_check=validate_transient_boundary_check(path,summary)
    return saved_validation_result(summary,transient_saved_validation_report(path,summary,boundary_check),boundary_check)

"""Reload result files and verify fields/metrics rather than trusting status text."""

import json
from pathlib import Path
import numpy as np

from thermal_opd.convergence import field_errors, within_tolerance


def require(condition, message):
    if not condition:
        raise ValueError(message)


def matching_errors(actual, stored, *, skip_pv=False):
    for key, value in actual.items():
        if skip_pv and key == 'opd_relative_pv':
            continue
        require(np.allclose(value, stored[key], rtol=1e-8, atol=1e-12), f'Saved field metric differs: {key}')


def validate_selection(path, summary):
    with np.load(path/'fields.npz', allow_pickle=False) as fields:
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
            matching_errors(errors, selected['errors'], skip_pv=True)
            fine = field_errors(fields[label+'_selected_temperature_K'], fields[label+'_reference_temperature_K'], weights,
                                fields[label+'_verification_selected_opd_m'], fields[label+'_verification_reference_opd_m'])
            matching_errors(fine, case['denser_common_grid_errors'], skip_pv=True)
        require(summary['selected_figure_xy'] is not None, 'No accepted figure grid')
        for group in summary['sampling_checks']:
            row = next(r for r in group['cases'] if r['xy_count']==summary['selected_figure_xy'])
            require(row['relative_pv_sampling_error'] <= summary['tolerance'], 'Figure sampling failed')


def boundary_edges(triangles):
    edges = np.concatenate((triangles[[0, 1]].T, triangles[[1, 2]].T, triangles[[2, 0]].T))
    edges, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    return edges[counts==1]


def inside_polygon(x, y, p, edges):
    inside = np.ones(np.broadcast_shapes(np.shape(x), np.shape(y)), dtype=bool)
    scale = np.max(np.abs(p))
    for edge in edges:
        a, b = p[:, edge].T
        normal = np.array([b[1]-a[1], a[0]-b[0]])
        normal *= 1 if np.dot(normal, a+b) >= 0 else -1
        inside &= normal[0]*(x-a[0])+normal[1]*(y-a[1]) <= 64*np.finfo(float).eps*scale**2
    return inside


def validate_fem(path, summary):
    p, h = summary['physical_parameters_SI'], summary['h_eff_W_per_m2_K']
    require(summary['sanity_pass'], 'FEM sanity failed')
    require(len(summary['cases']) == len(summary['angles_deg']), 'Missing FEM angle')
    for case in summary['cases']:
        label = f"theta{case['angle_deg']:g}"
        with np.load(path/f'{label}_fields.npz', allow_pickle=False) as fields:
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
                matching_errors(actual, row['errors_vs_sa'])
                changes = field_errors(tf[j], tf[-1], np.ones(ts.shape[1]), of[j], of[-1])
                matching_errors(changes, row['changes_vs_finest_fem'])
                require(np.allclose(np.ptp(of[j]), row['fem_opd_pv_common_aperture_m'], rtol=1e-10, atol=1e-20), 'FEM PV mismatch')
                require(row['relative_linear_residual'] <= summary['sanity_tolerance'], 'FEM linear residual failed')
                with np.load(path/f'{label}_mesh{j}.npz', allow_pickle=False) as mesh:
                    xy, tri, z, t = (mesh[k] for k in ('disk_nodes_xy_m', 'disk_triangles', 'axial_nodes_z_m', 'temperature_K'))
                    require(t.shape==(z.size, xy.shape[1]) and np.all(np.isfinite(t)), 'Bad FEM nodal temperature')
                    require(np.all(np.diff(z)>0) and z[0]==0 and z[-1]==p['thickness'], 'Bad FEM axial coordinates')
                    require(tri.min()>=0 and tri.max()<xy.shape[1], 'Bad FEM topology')
                    edges = boundary_edges(tri)
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
                    require(np.all(inside_polygon(tx, ty, xy, edges)), 'Temperature point outside FEM polygon')
                    require(np.all(inside_polygon(ox, oy, xy, edges) & inside_polygon(ox, oy+shift, xy, edges)), 'Ray outside FEM polygon')
            q_check = case['load_quadrature_check']
            actual_q = field_errors(fields['temperature_higher_quadrature_K'][..., tm], tf[-1], np.ones(ts.shape[1]),
                                    fields['opd_higher_quadrature_m'][om], of[-1])
            matching_errors(actual_q, q_check['errors'])
            require(q_check['passed'] and within_tolerance(q_check['errors'], summary['mesh_change_target']/10), 'FEM load quadrature failed')
            expected_change = within_tolerance(case['meshes'][-2]['changes_vs_finest_fem'], summary['mesh_change_target'])
            require(case['finest_mesh_change_within_target']==expected_change, 'Incorrect FEM convergence flag')
        require((path/f'{label}_comparison.png').is_file(), 'Missing FEM comparison figure')
    require((path/'comparison.csv').is_file(), 'Missing FEM comparison table')
    established = summary['mode']=='full' and all(c['finest_mesh_change_within_target'] for c in summary['cases'])
    require(summary['mesh_convergence_established']==established, 'Incorrect full mesh-convergence flag')


def validate_transient(path, summary):
    """Recompute mode sums, field errors and nodal storage/loss from saved files."""
    p,h = summary['physical_parameters_SI'],summary['h_eff_W_per_m2_K']
    tolerance = summary['sanity_tolerance']
    require(summary['sanity_pass'],'Transient sanity failed')
    require(not summary['transient_truncation_established'] and not summary['fem_convergence_established'],
            'Diagnostic study cannot claim final convergence')
    with np.load(path/'fields.npz',allow_pickle=False) as f:
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
        matching_errors(long_errors,summary['long_time_errors_vs_steady'])
        for j,t in enumerate(times):
            if t==0:
                continue
            for i in range(len(summary['axial_counts'])):
                err = field_errors(f['axial_check_temperature_K'][i,j].reshape(3,-1),temperature[j].reshape(3,-1),weights.ravel(),
                                   f['axial_check_opd_m'][i,j],opd[j])
                matching_errors(err,summary['axial_errors_vs_largest_J'][i][j])
            err = field_errors(temperature[j].reshape(3,-1),f['reference_temperature_K'][j].reshape(3,-1),weights.ravel(),
                               opd[j],f['reference_opd_m'][j])
            matching_errors(err,summary['transverse_errors_vs_finite_reference'][j])
            err = field_errors(temperature[j].reshape(3,-1),temperature[j].reshape(3,-1),weights.ravel(),
                               opd[j],f['higher_ray_quadrature_opd_m'][j])
            matching_errors(err,summary['ray_quadrature_errors'][j])
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
                    matching_errors(err,row['errors_vs_sa'][j])
            with np.load(path/f'fem_run{i}.npz',allow_pickle=False) as mesh:
                xy,tri,z,T = (mesh[key] for key in ('disk_nodes_xy_m','disk_triangles','axial_nodes_z_m','temperature_K'))
                require(T.shape==(times.size,z.size,xy.shape[1]) and np.all(np.isfinite(T)),'Bad transient FEM nodal fields')
                require(np.array_equal(mesh['times_s'],times) and np.all(T[0]==0),'Bad FEM initial time')
                require(z[0]==0 and z[-1]==p['thickness'] and np.all(np.diff(z)>0),'Bad FEM layers')
                require(tri.min()>=0 and tri.max()<xy.shape[1],'Bad FEM mesh topology')
                v1,v2 = xy[:,tri[1]]-xy[:,tri[0]],xy[:,tri[2]]-xy[:,tri[0]]
                areas = np.abs(v1[0]*v2[1]-v1[1]*v2[0])/2
                edges = boundary_edges(tri)
                lengths = np.linalg.norm(xy[:,edges[:,1]]-xy[:,edges[:,0]],axis=0)
                stored = p['density']*p['specific_heat']*np.trapezoid(np.sum(T[:,:,tri].mean(axis=-2)*areas,axis=-1),x=z,axis=-1)
                require(np.all(inside_polygon(x[common],y[common],xy,edges)&inside_polygon(x[common],y[common]+shift,xy,edges)),
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
        require((path/'transient_modes.png').is_file() and (path/'modal_contributions.csv').is_file(),'Missing modal study artifacts')


def validate_transient_selection(path,summary):
    require(len(summary['cases'])==len(summary['angles_deg']),'Missing transient selection angle')
    with np.load(path/'fields.npz',allow_pickle=False) as f:
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
                    matching_errors(err,stored[j])
                    require(within_tolerance(err,limit),'Transient reference/candidate above target')
                err = field_errors(f[label+'_selected_temperature_K'][j],f[label+'_reference_temperature_K'][j],weights,
                                   f[label+'_fine_selected_opd_m'][j],f[label+'_fine_reference_opd_m'][j])
                err['opd_relative_pv'] = abs(f[label+'_selected_continuous_pv_m'][j]-f[label+'_reference_continuous_pv_m'][j])/f[label+'_reference_continuous_pv_m'][j]
                matching_errors(err,case['denser_errors_by_time'][j])
                require(within_tolerance(err,summary['tolerance']),'Dense transient comparison above target')
                pv = f[label+'_selected_continuous_pv_m'][j]
                sampling = abs(np.ptp(f[label+'_figure_opd_m'][j])-pv)/pv
                require(np.isclose(sampling,case['sampling_errors_by_time'][j],rtol=1e-9,atol=1e-12),'Transient sampling metric differs')
                require(sampling<=summary['tolerance'],'Transient figure sampling above target')


def validate_fine_transient(path,summary):
    p = summary['physical_parameters_SI']
    if summary['transient_truncation_established']:
        selected_path = Path(summary['selection_directory'])
        validate_saved_run(selected_path)
        selected = json.loads((selected_path/'summary.json').read_text(encoding='utf-8'))
        require(selected['kind']=='transient_joint_truncation_selection' and selected['physical_parameters_SI']==p,
                'Fine transient selection has different physics')
        require(summary['angle_deg'] in selected['angles_deg'] and summary['modal_settings']==selected['selected_settings']['modal']
                and summary['selected_axial_modes']==selected['selected_settings']['J'],'Fine transient outside selected scope')
        require(summary['precision_scope']['times_s']==selected['times_s'],'Selected times differ')
    require(not summary['fem_convergence_established'],'Fine SA fields do not establish FEM convergence')
    with np.load(path/'fields.npz',allow_pickle=False) as f:
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
            from thermal_opd.parameters import PhysicalParameters
            from thermal_opd.transient import TransientModel
            coefficients = {int(m):f[f'source_coefficients_m{m}_W_per_m2'] for m in f['orders']}
            roots = {int(m):f[f'radial_roots_m{m}_per_m'] for m in f['orders']}
            model = TransientModel(PhysicalParameters(**p),coefficients,roots,num_axial_modes=summary['selected_axial_modes'])
            ii,jj = np.array([0,r.size//2,r.size-1]),np.array([0,phi.size//4,phi.size//2])
            for j in (1,times.size-2,times.size-1):
                for plane,depth in enumerate(z):
                    expected = model.temperature(r[ii],phi[jj],depth,times[j])
                    require(np.allclose(T[j,plane,ii,jj],expected,rtol=1e-11,atol=1e-13),'Fine field does not match saved model')
            w = np.broadcast_to(np.maximum(r[:,None],p['radius']*1e-9),T.shape[-2:]).ravel()
            errors = field_errors(T[-1].reshape(3,-1),f['steady_temperature_K'].reshape(3,-1),w,O[-1,mask],f['steady_opd_m'][mask])
            matching_errors(errors,summary['long_time_errors_vs_same_modes_steady'])
            require(within_tolerance(errors,summary['precision_target']),'Fine long-time limit above target')
            with np.load(path/'history.npz',allow_pickle=False) as history:
                require(all(np.all(np.isfinite(history[key])) for key in history.files),'Nonfinite fine history')
                ht,ho,hT = history['times_s'],history['opd_m'],history['temperature_K']
                require(ht[0]==0 and np.all(np.diff(ht)>0) and ht[-1]==times[-1],'Invalid fine history times')
                require(np.allclose(np.ptp(ho,axis=1),history['sampled_pv_m'],rtol=1e-11,atol=1e-20),'Fine history PV differs')
                for j,time in enumerate(times):
                    k = np.flatnonzero(ht==time)
                    require(k.size==1,'Fine history missed snapshot')
                    require(np.allclose(hT[k[0],:,0],T[j,:,0,0],rtol=1e-11,atol=1e-13),'Fine center history differs')
                    require(np.allclose(ho[k[0]],O[j,mask],rtol=1e-11,atol=1e-20),'Fine history optical field differs')
    figure_path = Path(summary.get('figure_directory',path))
    for filename in summary['figures']:
        require((figure_path/filename).is_file(),'Missing fine transient figure')


def validate_saved_run(directory):
    path = Path(directory)
    summary = json.loads((path/'summary.json').read_text(encoding='utf-8'))
    require(summary.get('status')=='COMPLETE', 'Run is not complete')
    kind = summary.get('kind')
    if kind == 'steady_joint_truncation_selection':
        validate_selection(path, summary)
    elif kind == 'steady_prismatic_fem_comparison':
        validate_fem(path, summary)
    elif kind == 'transient_modal_fem_study':
        validate_transient(path, summary)
    elif kind == 'transient_joint_truncation_selection':
        validate_transient_selection(path,summary)
    elif kind == 'transient_fine_fields':
        validate_fine_transient(path,summary)
    else:
        raise ValueError(f'Unsupported run kind: {kind}')
    boundary_check = None
    if kind in ('transient_modal_fem_study','transient_fine_fields'):
        from thermal_opd.transient_validation import validate_transient_boundary_check
        boundary_check = validate_transient_boundary_check(path,summary)
    return {'status': 'PASS', 'kind': kind, 'run_directory': str(path.resolve()),
            'mesh_convergence_established': summary.get('mesh_convergence_established'),
            'transient_truncation_established': summary.get('transient_truncation_established'),
            'fem_convergence_established': summary.get('fem_convergence_established'),
            'transient_boundary_check': boundary_check,
            'note': 'PASS verifies saved run consistency and stated numerical gates; smoke PASS is not final convergence.'}

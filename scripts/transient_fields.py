"""Fine transient fields, called by run_transient_case --mode fields.

Numerical settings stay in this script. Shared physics and snapshot plotting
are imported from thermal_opd; this does not run a large FEM batch.
"""
import csv
from dataclasses import asdict
import json
from pathlib import Path
from time import perf_counter
import matplotlib.pyplot as plt
import numpy as np
from thermal_opd.convergence import continuous_aperture_pv,field_errors
from thermal_opd.model import ModalResolution,solve_modal_source
from thermal_opd.opd import refracted_angle,transmitted_ray_mask,steady_refracted_opd,prepare_modal_ray_series
from thermal_opd.plotting import save_temperature_opd_figure
from thermal_opd.run_io import create_run,provenance,write_summary
from thermal_opd.run_validation import validate_saved_run
from thermal_opd.steady import steady_temperature,axial_response
from thermal_opd.transient import TransientModel
from thermal_opd.transient_validation import write_transient_boundary_check

temperature_r_points,temperature_phi_points,opd_xy_points = 161,180,201
aperture_radius_ratio = .9
point_chunk_size = 512
history_positive_points = 121
long_time_decay_exponent = 12.
pv_seed_radius,pv_seed_phi = 9,48
boundary_z_points = 65
boundary_relative_tolerance,boundary_absolute_tolerance_W_per_m2 = 1e-10,1e-9


def load_selection(directory,p,angle,history):
    directory = Path(directory).resolve()
    validate_saved_run(directory)
    s = json.loads((directory/'summary.json').read_text(encoding='utf-8'))
    if s['kind']!='transient_joint_truncation_selection' or s['physical_parameters_SI']!=asdict(p):
        raise ValueError('Need completed transient selection with identical physical parameters')
    if angle not in s['angles_deg'] or history!='step' or s['history']['kind']!='step':
        raise ValueError('Selected precision scope does not cover this angle/loading history')
    return s


def run_fine_fields(args,p,root):
    if args.selection is None:
        raise ValueError('--mode fields requires --selection pointing to a completed transient selection')
    selection = load_selection(args.selection,p,args.angle,args.history)
    modal = ModalResolution(**selection['selected_settings']['modal'])
    J = selection['selected_settings']['J']
    if opd_xy_points<selection['selected_figure_xy'] or selection['aperture_radius_ratio']!=aperture_radius_ratio:
        raise ValueError('Fine figure settings do not satisfy the selection sampling/aperture')
    directory = create_run(root,f'transient_fields_theta{args.angle:g}',args.output_dir)
    figures = Path(root)/'outputs/figures'/directory.name
    figures.mkdir(parents=True,exist_ok=False)
    summary = dict(kind='transient_fine_fields',status='RUNNING',angle_deg=args.angle,physical_parameters_SI=asdict(p),
        aperture_radius_ratio=aperture_radius_ratio,selection_directory=str(Path(args.selection).resolve()),
        transient_truncation_established=True,fem_convergence_established=False,boundary_validation_required=True,
        modal_settings=asdict(modal),selected_axial_modes=J,precision_target=selection['tolerance'],
        precision_scope=dict(angles_deg=selection['angles_deg'],times_s=selection['times_s'],history=selection['history']),
        numerical_settings=dict(temperature_r=temperature_r_points,temperature_phi=temperature_phi_points,opd_xy=opd_xy_points),
        figure_directory=str(figures),provenance=provenance(root,__file__),snapshots=[],figures=[])
    write_summary(directory,summary)
    print(f'Fine fields: {directory}\nFigures: {figures}',flush=True)
    started = perf_counter()
    try:
        internal = refracted_angle(np.deg2rad(args.angle),p.incident_refractive_index,p.substrate_refractive_index)
        summary['refraction_angle_deg'] = float(np.rad2deg(internal))
        c,roots = solve_modal_source(p,args.angle,modal)
        model = TransientModel(p,c,roots,num_axial_modes=J)
        long_time = max(selection['times_s'][-1],long_time_decay_exponent/model.slowest_decay_rate)
        times = np.array([0.,*selection['times_s'],long_time])
        summary['times_s'],summary['long_time_s'] = times.tolist(),float(long_time)
        r = np.linspace(0.,p.radius,temperature_r_points)
        phi = np.arange(temperature_phi_points)*2*np.pi/temperature_phi_points
        z = np.array([0.,p.thickness/2,p.thickness])
        axis = np.linspace(-aperture_radius_ratio*p.radius,aperture_radius_ratio*p.radius,opd_xy_points)
        x,y = np.meshgrid(axis,axis)
        mask = (x*x+y*y<=(aperture_radius_ratio*p.radius)**2)&transmitted_ray_mask(x,y,p.radius,p.thickness,internal)
        print('Computing dense snapshots on the selected modes',flush=True)
        T = np.stack([np.stack([model.temperature(r[:,None],phi[None,:],plane,t) for plane in z]) for t in times])
        O = np.full((times.size,)+x.shape,np.nan)
        O[:,mask] = model.refracted_opd_series(x[mask],y[mask],times,radius=p.radius,refraction_angle=internal,
                                              num_ray_points=modal.ray_points,point_chunk_size=point_chunk_size)
        sampled_pv = np.ptp(O[:,mask],axis=1)
        saved = dict(times_s=times,temperature_r_m=r,temperature_phi_rad=phi,temperature_z_m=z,temperature_K=T,
                     opd_x_m=x,opd_y_m=y,aperture_mask=mask,opd_m=O,sampled_pv_m=sampled_pv,orders=np.array(list(c)))
        for m in c:
            saved[f'source_coefficients_m{m}_W_per_m2'],saved[f'radial_roots_m{m}_per_m'] = c[m],roots[m]
        # Source coefficients and axial/ray modules are exactly those selected.
        case = next(case for case in selection['cases'] if case['angle_deg']==args.angle)
        row = next(row for row in case['candidates'] if row['settings']==selection['selected_settings'])
        continuous = [None,*row['extrema']]
        evaluator = model.prepare_refracted_opd(long_time,radius=p.radius,refraction_angle=internal,
                                              num_ray_points=modal.ray_points,point_chunk_size=point_chunk_size)
        continuous.append(continuous_aperture_pv(evaluator,radius=p.radius,thickness=p.thickness,refraction_angle=internal,
            aperture_radius=aperture_radius_ratio*p.radius,num_seed_radius=pv_seed_radius,num_seed_phi=pv_seed_phi))
        for j,t in enumerate(times):
            if t==0:
                summary['snapshots'].append(dict(time_s=0.,sampled_pv_m=0.,continuous_pv_m=0.,figure=None))
                continue
            filename = f'transient_t{t:.6g}s.png'
            summary['snapshots'].append(dict(time_s=float(t),sampled_pv_m=float(sampled_pv[j]),
                continuous_pv_m=continuous[j]['pv'],extrema=continuous[j],figure=filename,
                relative_pv_sampling_error=float(abs(sampled_pv[j]-continuous[j]['pv'])/continuous[j]['pv'])))
            save_temperature_opd_figure(figures/filename,r,phi,T[j],x,y,mask,O[j],
                title=f'Transient: {args.angle:g} deg, t={t:.6g} s, M={modal.M}, N={modal.N}, J={J}',
                pv_m=sampled_pv[j],continuous_pv_m=continuous[j]['pv'])
            summary['figures'].append(filename)
        thermal_args = (c,roots,p.conductivity,p.thickness,p.h_eff,p.h_eff)
        saved['steady_temperature_K'] = np.stack([steady_temperature(r[:,None],phi[None,:],plane,*thermal_args) for plane in z])
        saved['steady_opd_m'] = np.full(x.shape,np.nan)
        saved['steady_opd_m'][mask] = steady_refracted_opd(x[mask],y[mask],*thermal_args,p.dn_dT,radius=p.radius,
            refraction_angle=internal,num_ray_points=modal.ray_points,point_chunk_size=point_chunk_size)
        steady_integrator = prepare_modal_ray_series(c,roots,
            lambda m,n,depth:axial_response(depth,roots[m][n],c[m][n],p.conductivity,p.thickness,p.h_eff,p.h_eff)[0][None,:],
            1,p.thickness,p.dn_dT,radius=p.radius,refraction_angle=internal,num_ray_points=modal.ray_points,
            point_chunk_size=point_chunk_size)
        steady_extrema = continuous_aperture_pv(lambda xx,yy:steady_integrator(xx,yy)[0],radius=p.radius,
            thickness=p.thickness,refraction_angle=internal,aperture_radius=aperture_radius_ratio*p.radius,
            num_seed_radius=pv_seed_radius,num_seed_phi=pv_seed_phi)
        summary['steady_extrema'] = steady_extrema
        save_temperature_opd_figure(figures/'steady_limit.png',r,phi,saved['steady_temperature_K'],x,y,mask,saved['steady_opd_m'],
            title=f'Same-mode steady limit: {args.angle:g} deg, M={modal.M}, N={modal.N}',
            pv_m=np.ptp(saved['steady_opd_m'][mask]),continuous_pv_m=steady_extrema['pv'])
        summary['figures'].append('steady_limit.png')
        long_errors = field_errors(T[-1].reshape(3,-1),saved['steady_temperature_K'].reshape(3,-1),
            np.broadcast_to(np.maximum(r[:,None],p.radius*1e-9),T.shape[-2:]).ravel(),O[-1,mask],saved['steady_opd_m'][mask])
        summary['long_time_errors_vs_same_modes_steady'] = long_errors
        # Same fine coordinates for the retained thermal order contributions.
        order_T,order_O = [],[]
        for m in c:
            one = TransientModel(p,{m:c[m]},{m:roots[m]},num_axial_modes=J)
            order_T.append(np.stack([np.stack([one.temperature(r[:,None],phi[None,:],plane,t) for plane in z]) for t in times]))
            optical = np.full_like(O,np.nan)
            optical[:,mask] = one.refracted_opd_series(x[mask],y[mask],times,radius=p.radius,refraction_angle=internal,
                num_ray_points=modal.ray_points,point_chunk_size=point_chunk_size)
            order_O.append(optical)
        saved['order_temperature_K'],saved['order_opd_m'] = np.array(order_T),np.array(order_O)
        wr = r*(r[1]-r[0]);wr[-1] *= .5
        energies = np.sum(saved['order_temperature_K']**2*wr[None,None,None,:,None],axis=(-2,-1))
        shares = np.divide(energies,energies.sum(axis=0),out=np.zeros_like(energies),where=energies.sum(axis=0)!=0)
        saved['temperature_plane_l2_squared_shares'] = shares
        saved['order_sampled_opd_pv_m'] = np.ptp(saved['order_opd_m'][...,mask],axis=-1)
        with (directory/'modal_contributions.csv').open('w',newline='',encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['time_s','thermal_order_m','front_temperature_L2_squared_share','single_order_sampled_opd_PV_nm','total_sampled_PV_nm'])
            for j,t in enumerate(times):
                for i,m in enumerate(c):
                    writer.writerow([t,m,shares[i,j,0],saved['order_sampled_opd_pv_m'][i,j]*1e9,sampled_pv[j]*1e9])
        fig,axes = plt.subplots(1,2,figsize=(12,4),layout='constrained')
        for i,m in enumerate(c):
            axes[0].semilogx(times[1:],shares[i,1:,0],marker='o',label=f'm={m}')
            axes[1].semilogx(times[1:],saved['order_sampled_opd_pv_m'][i,1:]*1e9,marker='o',label=f'thermal m={m}')
        axes[1].semilogx(times[1:],sampled_pv[1:]*1e9,'k--',label='total')
        axes[0].set_ylabel('Front temperature L2 squared share')
        axes[1].set_ylabel('Single-order sampled OPD PV [nm]')
        for ax in axes:
            ax.set_xlabel('Time [s]');ax.grid(alpha=.25);ax.legend()
        fig.savefig(figures/'transient_modes.png',dpi=180);plt.close(fig)
        summary['figures'].append('transient_modes.png')
        np.savez_compressed(directory/'fields.npz',**saved)
        print('Checking front/rear/side Robin residuals at every snapshot',flush=True)
        boundary = write_transient_boundary_check(directory,model,times,r,phi,np.linspace(0,p.thickness,boundary_z_points),
            angle_deg=args.angle,relative_tolerance=boundary_relative_tolerance,
            absolute_tolerance=boundary_absolute_tolerance_W_per_m2,provenance=summary['provenance'])
        if boundary['status']!='PASS':
            raise RuntimeError('Transient Robin boundary residual failed')
        print('Computing temperature histories and dense-aperture PV history',flush=True)
        ht = np.unique(np.r_[0.,times,np.geomspace(times[1],long_time,history_positive_points)])
        hx,hy = np.array([0.,p.radius/2,0.]),np.array([0.,0.,p.radius/2])
        history_T = np.stack([np.stack([model.temperature(np.hypot(hx,hy),np.arctan2(hy,hx),plane,t) for plane in z]) for t in ht])
        history_O = model.refracted_opd_series(x[mask],y[mask],ht,radius=p.radius,refraction_angle=internal,
            num_ray_points=modal.ray_points,point_chunk_size=point_chunk_size)
        history_pv = np.ptp(history_O,axis=1)
        history_steady_T = np.stack([steady_temperature(np.hypot(hx,hy),np.arctan2(hy,hx),plane,*thermal_args) for plane in z])
        np.savez_compressed(directory/'history.npz',times_s=ht,point_x_m=hx,point_y_m=hy,z_m=z,
            temperature_K=history_T,opd_m=history_O,sampled_pv_m=history_pv,steady_temperature_K=history_steady_T,
            steady_sampled_pv_m=np.ptp(saved['steady_opd_m'][mask]))
        with (directory/'snapshots.csv').open('w',newline='',encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['time_s','center_front_K','center_middle_K','center_rear_K','sampled_PV_nm','continuous_numerical_PV_nm'])
            for j,t in enumerate(times):
                writer.writerow([t,*T[j,:,0,0],sampled_pv[j]*1e9,summary['snapshots'][j]['continuous_pv_m']*1e9])
        with (directory/'history.csv').open('w',newline='',encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['time_s','center_front_K','center_middle_K','center_rear_K','x_half_radius_front_K','y_half_radius_front_K','sampled_PV_nm'])
            for j,t in enumerate(ht):
                writer.writerow([t,*history_T[j,:,0],history_T[j,0,1],history_T[j,0,2],history_pv[j]*1e9])
        save_evolution(figures,saved,times.size-1)
        save_early_depths(figures,saved)
        save_history(figures,ht,history_T,history_pv,history_steady_T,np.ptp(saved['steady_opd_m'][mask]),summary['snapshots'])
        summary['figures'] += ['transient_evolution.png','transient_early_depths.png','transient_history.png']
        summary.update(status='COMPLETE',elapsed_seconds=perf_counter()-started,
            precision_note='0.1% finite-reference selection at recorded angles/key times; zero exact, long-time steady limit checked separately. Intermediate history samples are not individually certified.',
            fem_note='This is the fine semi-analytical field run. FEM space/time convergence remains separate.')
        write_summary(directory,summary)
        print(f'COMPLETE: {directory}',flush=True)
    except BaseException as error:
        summary.update(status='ABORTED' if isinstance(error,KeyboardInterrupt) else 'FAILED',error=str(error))
        write_summary(directory,summary)
        raise


def save_evolution(directory,f,last):
    indices = range(1,last)
    fig,axes = plt.subplots(2,len(indices),figsize=(19,7),layout='constrained')
    r,phi = f['temperature_r_m'],np.r_[f['temperature_phi_rad'],2*np.pi]
    tx,ty = r[:,None]*np.cos(phi)*1e3,r[:,None]*np.sin(phi)*1e3
    temperature_high = np.max(f['temperature_K'][list(indices),0])
    opd_high = max(np.ptp(f['opd_m'][j][f['aperture_mask']])*1e9 for j in indices)
    for column,j in enumerate(indices):
        T = f['temperature_K'][j,0]
        contour = axes[0,column].contourf(tx,ty,np.column_stack((T,T[:,0])),levels=np.linspace(0,temperature_high,61),cmap='inferno',extend='min')
        O = f['opd_m'][j]
        optical = (O-np.min(O[f['aperture_mask']]))*1e9
        ray = axes[1,column].contourf(f['opd_x_m']*1e3,f['opd_y_m']*1e3,optical,levels=np.linspace(0,opd_high,61),cmap='viridis')
        axes[0,column].set_title(f't={f["times_s"][j]:g} s')
        for ax in axes[:,column]:
            ax.set_aspect('equal');ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]')
    fig.colorbar(contour,ax=axes[0,:],label='Front temperature rise [K]')
    fig.colorbar(ray,ax=axes[1,:],label='OPD minus aperture minimum [nm]')
    fig.suptitle('Transient evolution: fixed scales across time within each row')
    fig.savefig(directory/'transient_evolution.png',dpi=180);plt.close(fig)


def save_history(directory,t,T,pv,steady_T,steady_pv,snapshots):
    fig,axes = plt.subplots(2,2,figsize=(12,8),layout='constrained')
    early = t<=75.
    colors = ('tab:red','tab:orange','tab:blue')
    for i,(name,color) in enumerate(zip(('front','middle','rear'),colors)):
        axes[0,0].plot(t[early],T[early,i,0],color=color,label=f'center {name}')
        axes[0,1].semilogx(t[1:],T[1:,i,0],color=color,label=f'center {name}')
        axes[0,1].axhline(steady_T[i,0],color=color,linestyle='--',alpha=.6)
    axes[1,0].plot(t[early],pv[early]*1e9,label='201x201 sampled PV')
    rows = [row for row in snapshots if 0<row['time_s']<=75]
    axes[1,0].scatter([row['time_s'] for row in rows],[row['continuous_pv_m']*1e9 for row in rows],
                      color='black',s=20,label='continuous numerical PV')
    axes[1,1].semilogx(t[1:],pv[1:]*1e9,label='201x201 sampled PV')
    axes[1,1].axhline(steady_pv*1e9,color='black',linestyle='--',label='same modes steady')
    for i,ax in enumerate(axes.flat):
        ax.set_xlabel('Time [s]');ax.set_ylabel('Temperature rise [K]' if i<2 else 'OPD PV [nm]')
        ax.grid(alpha=.25);ax.legend()
    axes[0,0].set_title('Early center temperature history')
    axes[0,1].set_title('Approach to steady temperature')
    axes[1,0].set_title('Early refracted-ray OPD history')
    axes[1,1].set_title('Approach to steady OPD')
    fig.savefig(directory/'transient_history.png',dpi=180);plt.close(fig)


def save_early_depths(directory,f):
    """Separate scales expose the extremely small early rear temperature."""
    r,phi = f['temperature_r_m'],np.r_[f['temperature_phi_rad'],2*np.pi]
    x,y = r[:,None]*np.cos(phi)*1e3,r[:,None]*np.sin(phi)*1e3
    fig,axes = plt.subplots(2,3,figsize=(13,8),layout='constrained')
    for row,j in enumerate((1,2)):
        for plane,name in enumerate(('front','middle','rear')):
            T = f['temperature_K'][j,plane]
            contour = axes[row,plane].contourf(x,y,np.column_stack((T,T[:,0])),levels=61,cmap='inferno')
            axes[row,plane].set_title(f't={f["times_s"][j]:g} s, {name}')
            axes[row,plane].set_aspect('equal');axes[row,plane].set_xlabel('x [mm]');axes[row,plane].set_ylabel('y [mm]')
            fig.colorbar(contour,ax=axes[row,plane],label='Temperature rise [K]',shrink=.85)
    fig.suptitle('Early depth response: each panel uses its own temperature scale')
    fig.savefig(directory/'transient_early_depths.png',dpi=180);plt.close(fig)

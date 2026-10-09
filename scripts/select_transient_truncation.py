"""Joint transient truncation at the steady study's 0.1% numerical target.

Step loading, nominal physical inputs, 0/30/60 degrees and specified positive
observation times only. No FEM is run. Numerical continuous PV uses the same
multistart optimizer as the steady selection; it is not a certified bound.
"""
from thermal_opd.run_io import script_defaults
from thermal_opd.model import ThermalModel
from thermal_opd.run_io import data_directory
import argparse
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import numpy as np
from numpy.polynomial.legendre import leggauss
from thermal_opd.convergence import field_errors,within_tolerance,continuous_aperture_pv
from thermal_opd.model import ModalResolution
from thermal_opd.opd import refracted_angle,transmitted_ray_mask
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.run_io import create_run,provenance,write_summary
from thermal_opd.transient import PowerHistory

parameters = PhysicalParameters()
angles_deg = (0.,30.,60.)
times_s = (.625,2.5,20.,50.,75.)
tolerance = 1e-3
reference_tolerance = tolerance/5
reference = ModalResolution(14,32,2400,256,64)
enlarged_reference = ModalResolution(16,40,3200,512,80)
reference_J,enlarged_reference_J = 48,64
candidates = (ModalResolution(6,8,400,32,16),ModalResolution(8,12,800,64,16),
              ModalResolution(8,16,1200,128,16),ModalResolution(10,20,1200,128,32),
              ModalResolution(12,24,1600,128,32))
candidate_J = 12
temperature_r_points,temperature_phi_points = 96,128
common_xy_points,verification_xy_points,figure_xy_points = 61,121,201
aperture_radius_ratio = .9
pv_seed_radius,pv_seed_phi = 9,48
point_chunk_size = 512
repository_root = Path(__file__).resolve().parents[1]


def optical_coordinates(p,angle,count):
    axis = np.linspace(-aperture_radius_ratio*p.radius,aperture_radius_ratio*p.radius,count)
    x,y = np.meshgrid(axis,axis)
    mask = (x*x+y*y<=(aperture_radius_ratio*p.radius)**2)&transmitted_ray_mask(x,y,p.radius,p.thickness,angle)
    return x[mask],y[mask]


def optical(model,xy,angle,ray_points):
    return model.opd_series(*xy, times_s, num_ray_points=ray_points, point_chunk_size=point_chunk_size)


def compare(T,R,weights,O,RO,pv=None,ref_pv=None):
    errors = [field_errors(T[i],R[i],weights,O[i],RO[i]) for i in range(len(times_s))]
    if pv is not None:
        for i,error in enumerate(errors):
            error['opd_relative_pv'] = float(abs(pv[i]-ref_pv[i])/ref_pv[i])
    return errors


def extrema(model,angle,ray_points):
    output = []
    for t in times_s:
        evaluator = model.prepare_opd(t, num_ray_points=ray_points, point_chunk_size=point_chunk_size)
        result = continuous_aperture_pv(evaluator,radius=parameters.radius,thickness=parameters.thickness,
            refraction_angle=angle,aperture_radius=aperture_radius_ratio*parameters.radius,
            num_seed_radius=pv_seed_radius,num_seed_phi=pv_seed_phi)
        if result['successful_optimizations']<1:
            raise RuntimeError('No successful continuous-PV optimization')
        output.append(result)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    args = parser.parse_args()
    directory = create_run(repository_root,'transient_truncation',args.output_dir)
    p = parameters
    summary = dict(kind='transient_joint_truncation_selection',status='RUNNING',physical_parameters_SI=asdict(p),
        h_eff_W_per_m2_K=p.h_eff,history=asdict(PowerHistory()),angles_deg=list(angles_deg),times_s=list(times_s),
        tolerance=tolerance,reference_stability_tolerance=reference_tolerance,
        reference={'modal':asdict(reference),'J':reference_J},
        enlarged_reference={'modal':asdict(enlarged_reference),'J':enlarged_reference_J},
        aperture_radius_ratio=aperture_radius_ratio,common_xy_points=common_xy_points,
        verification_xy_points=verification_xy_points,selected_figure_xy=figure_xy_points,
        provenance=provenance(repository_root,__file__),cases=[])
    write_summary(directory, summary, defaults=script_defaults(globals()))
    print(f'Run: {directory}',flush=True)
    started = perf_counter()
    nodes,w = leggauss(temperature_r_points)
    r = (nodes+1)*p.radius/2
    phi = np.arange(temperature_phi_points)*2*np.pi/temperature_phi_points
    rt,pt = np.broadcast_arrays(r[:,None],phi[None,:])
    weights = np.broadcast_to((w*p.radius/2*r)[:,None],rt.shape).ravel()
    planes = np.array([0.,p.thickness/2,p.thickness])
    saved = dict(temperature_r_m=rt.ravel(),temperature_phi_rad=pt.ravel(),temperature_z_m=planes,
                 temperature_area_weights=weights,times_s=np.array(times_s))
    retained = []
    try:
        for angle in angles_deg:
            internal = refracted_angle(np.deg2rad(angle),p.incident_refractive_index,p.substrate_refractive_index)
            xy = optical_coordinates(p,internal,common_xy_points)
            label = f'theta{angle:g}'
            def evaluate(settings,J):
                tick = perf_counter()
                prepared = ThermalModel(p,angle,settings).prepare()
                model = prepared.solve_transient(num_axial_modes=J)
                T = np.stack([np.stack([model.temperature(r[:,None],phi[None,:],z,t).ravel() for z in planes]) for t in times_s])
                O = optical(model,xy,internal,settings.ray_points)
                print(f'  fields M{settings.M}/N{settings.N}/J{J} ({perf_counter()-tick:.1f}s)',flush=True)
                return model,T,O
            print(f'Angle {angle:g}: two high-order references',flush=True)
            rm,RT,RO = evaluate(reference,reference_J)
            bm,BT,BO = evaluate(enlarged_reference,enlarged_reference_J)
            print(f'Angle {angle:g}: reference continuous PV',flush=True)
            rp,bp = extrema(rm,internal,reference.ray_points),extrema(bm,internal,enlarged_reference.ray_points)
            rpv,bpv = np.array([v['pv'] for v in rp]),np.array([v['pv'] for v in bp])
            stability = compare(RT,BT,weights,RO,BO,rpv,bpv)
            case = dict(angle_deg=angle,refraction_angle_deg=float(np.rad2deg(internal)),
                        reference_errors_by_time=stability,reference_pass=all(within_tolerance(e,reference_tolerance) for e in stability),
                        reference_extrema=rp,enlarged_extrema=bp,candidates=[])
            print(f'  reference pass={case["reference_pass"]}',flush=True)
            entries = []
            for settings in candidates:
                model,T,O = evaluate(settings,candidate_J)
                errors = compare(T,RT,weights,O,RO)
                # Reject inexpensive field failures before expensive extrema.
                ep = extrema(model,internal,settings.ray_points) if all(within_tolerance(e,tolerance) for e in errors) else None
                if ep is not None:
                    pv = np.array([v['pv'] for v in ep])
                    errors = compare(T,RT,weights,O,RO,pv,rpv)
                accepted = ep is not None and all(within_tolerance(e,tolerance) for e in errors)
                row = dict(settings={'modal':asdict(settings),'J':candidate_J},errors_by_time=errors,accepted=accepted,extrema=ep)
                case['candidates'].append(row)
                entries.append((model,T,O,ep))
                print(f'  candidate M{settings.M}/N{settings.N}: accepted={accepted}; maximum change={max(np.max(v) for e in errors for v in e.values()):.4g}',flush=True)
            retained.append((internal,rm,RT,RO,entries))
            for key,value in (('reference_temperature_K',RT),('enlarged_temperature_K',BT),('reference_opd_m',RO),('enlarged_opd_m',BO),
                              ('reference_continuous_pv_m',rpv),('enlarged_continuous_pv_m',bpv)):
                saved[label+'_'+key] = value
            saved[label+'_opd_x_m'],saved[label+'_opd_y_m'] = xy
            summary['cases'].append(case)
            write_summary(directory, summary, defaults=script_defaults(globals()))
        feasible = [i for i in range(len(candidates)) if all(case['reference_pass'] and case['candidates'][i]['accepted'] for case in summary['cases'])]
        if not feasible:
            raise RuntimeError('No candidate passes every angle/time and reference-stability gate')
        index = feasible[0]
        selected = candidates[index]
        summary['selected_settings'] = {'modal':asdict(selected),'J':candidate_J}
        for case,(internal,rm,RT,RO,entries) in zip(summary['cases'],retained):
            label = f'theta{case["angle_deg"]:g}'
            model,T,O,ep = entries[index]
            pv = np.array([v['pv'] for v in ep])
            print(f'{label}: verify 121x121 and 201x201 sampling',flush=True)
            xy_fine = optical_coordinates(p,internal,verification_xy_points)
            fine_r,fine_s = optical(rm,xy_fine,internal,reference.ray_points),optical(model,xy_fine,internal,selected.ray_points)
            fine_errors = compare(T,RT,weights,fine_s,fine_r,pv,saved[label+'_reference_continuous_pv_m'])
            case['denser_errors_by_time'] = fine_errors
            case['denser_common_grid_pass'] = all(within_tolerance(e,tolerance) for e in fine_errors)
            figure = optical(model,optical_coordinates(p,internal,figure_xy_points),internal,selected.ray_points)
            sampling = abs(np.ptp(figure,axis=1)-pv)/pv
            case['sampling_errors_by_time'],case['sampling_pass'] = sampling.tolist(),bool(np.all(sampling<=tolerance))
            for key,value in (('selected_temperature_K',T),('selected_opd_m',O),('selected_continuous_pv_m',pv),
                              ('fine_reference_opd_m',fine_r),('fine_selected_opd_m',fine_s),('figure_opd_m',figure)):
                saved[label+'_'+key] = value
            write_summary(directory, summary, defaults=script_defaults(globals()))
            if not (case['denser_common_grid_pass'] and case['sampling_pass']):
                raise RuntimeError(f'{label}: dense-field or figure-sampling gate failed')
        np.savez_compressed(data_directory(directory)/'fields.npz',**saved)
        summary.update(status='COMPLETE',elapsed_seconds=perf_counter()-started)
        write_summary(directory, summary, defaults=script_defaults(globals()))
        print(f'COMPLETE: selected {summary["selected_settings"]}; {directory}',flush=True)
    except BaseException as error:
        np.savez_compressed(data_directory(directory)/'partial_fields.npz',**saved)
        summary.update(status='ABORTED' if isinstance(error,KeyboardInterrupt) else 'FAILED',error=str(error))
        write_summary(directory, summary, defaults=script_defaults(globals()))
        raise


if __name__=='__main__':
    main()

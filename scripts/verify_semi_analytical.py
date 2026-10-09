"""Opt-in independent semi-analytical checks; numerical settings below.

Reports unexecuted convergence studies explicitly. Never launched by solve().
"""
from thermal_opd.run_io import script_defaults
from thermal_opd.run_io import data_directory
import argparse
from dataclasses import asdict,replace
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from thermal_opd.model import ThermalModel,PreparedModel,ModalResolution
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.transient import PowerHistory
from thermal_opd.run_io import create_run,provenance,write_summary
from thermal_opd.validation import transient_boundary_fields,steady_boundary_residuals,steady_energy_balance
from thermal_opd.verification import verification_check,verification_report,finite_difference_pde_residual

parameters=PhysicalParameters()
resolution=ModalResolution(6,8,400,32,32)
num_axial_modes=12
times_s=(.625,2.5,20.,50.,75.)
early_time_s=.001
early_axial_counts=(12,48,128,256,512)
reference_thickness_m,reference_time_s,reference_flux_W_per_m2=.1,.1,1000.
reference_radial_wavenumber_per_m,reference_axial_modes=1.,1200
reference_tolerance=2e-6
pde_tolerance,energy_tolerance,limit_tolerance=1e-4,1e-5,1e-3
boundary_tolerance_W_per_m2=1e-7
space_step_ratio,time_step_s=2e-4,1e-4
repository_root=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--angle',type=float,default=60.)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    directory=create_run(repository_root,'semi_analytical_verification',args.output_dir)
    p=parameters
    prepared=ThermalModel(p,args.angle,resolution).prepare()
    steady=prepared.solve_steady()
    transient=prepared.solve_transient(num_axial_modes=num_axial_modes,history=PowerHistory())
    config=dict(physical_parameters_SI=asdict(p),angle_deg=args.angle,resolution=asdict(resolution),
                num_axial_modes=num_axial_modes,times_s=list(times_s),early_time_s=early_time_s,
                early_axial_counts=list(early_axial_counts),space_step_ratio=space_step_ratio,time_step_s=time_step_s)
    checks=[];fields={}
    def check(name,value,threshold,scope,field=None):
        checks.append(verification_check(name,'passed' if value<=threshold else 'failed',
            metrics={'max_error':float(value)},threshold=threshold,scope=scope,evidence=field))
    r,phi,z=np.linspace(0,p.radius,25),np.arange(32)*2*np.pi/32,np.linspace(0,p.thickness,33)
    bc=transient_boundary_fields(transient,np.array([0.,*times_s]),r,phi,z)
    for name,key in (('front_boundary','front_modal_residual_W_per_m2'),('rear_boundary','rear_residual_W_per_m2'),('side_boundary','side_residual_W_per_m2')):
        fields[name]=bc[key]
        check(name,np.max(abs(bc[key])),boundary_tolerance_W_per_m2,'Retained heat flux; sampled Robin consistency',name)
    sb=steady_boundary_residuals(r,phi,z,prepared.coefficients,prepared.roots,p.radius,p.conductivity,p.thickness,p.h_eff,p.h_eff,p.h_eff)
    fields['steady_boundaries']=np.array(list(sb.values()))
    check('steady_boundaries',max(sb.values()),boundary_tolerance_W_per_m2,'Steady retained heat flux','steady_boundaries')
    rr,pp,zz=np.meshgrid(np.linspace(.15*p.radius,.8*p.radius,5),np.arange(8)*2*np.pi/8,
                         np.linspace(.15*p.thickness,.85*p.thickness,5),indexing='ij')
    rr,pp,zz=rr.ravel(),pp.ravel(),zz.ravel()
    fields['initial_temperature']=transient.temperature(rr,pp,zz,0.)
    check('initial_condition',np.max(abs(fields['initial_temperature'])),0.,'Interior volume samples; explicit zero initial branch','initial_temperature')
    peak=2*p.absorbed_power*np.cos(np.deg2rad(args.angle))/(np.pi*p.beam_radius**2)
    scale=max(peak/p.thickness,1.)
    h=space_step_ratio*min(p.radius,p.thickness)
    residuals=[]
    for step in (h,h/2):
        residuals.append(np.stack([finite_difference_pde_residual(transient.temperature,p,rr,pp,zz,t,
            space_step=step,time_step=time_step_s)/scale for t in times_s]))
    fields['pde_normalized_residuals']=np.array(residuals)
    check('transient_pde',np.max(abs(fields['pde_normalized_residuals'])),pde_tolerance,
          'Interior finite differences at h and h/2; normalized by full-power peak flux/L','pde_normalized_residuals')
    steady_res=np.stack([finite_difference_pde_residual(lambda r,phi,z,t:steady.temperature(r,phi,z),p,rr,pp,zz,1.,
        space_step=step,time_step=time_step_s)/scale for step in (h,h/2)])
    fields['steady_pde_normalized_residuals']=steady_res
    check('steady_pde',np.max(abs(steady_res)),pde_tolerance,'Independent interior finite differences','steady_pde_normalized_residuals')
    long_time=max(12/transient.slowest_decay_rate,max(times_s))
    st=steady.temperature(rr,pp,zz)
    lt=transient.temperature(rr,pp,zz,long_time)
    fields['long_time_normalized_difference']=(lt-st)/max(np.max(abs(st)),np.finfo(float).tiny)
    check('long_time_steady_limit',np.max(abs(fields['long_time_normalized_difference'])),limit_tolerance,
          f'Step; t={long_time:g}s; independent steady solution interface','long_time_normalized_difference')
    zero_p=replace(p,absorbed_power=0.)
    zero=PreparedModel(zero_p,args.angle,resolution,{m:np.zeros_like(v) for m,v in prepared.coefficients.items()},prepared.roots)
    zero_t=zero.solve_transient(num_axial_modes=num_axial_modes)
    fields['zero_source']=np.r_[zero.solve_steady().temperature(rr,pp,zz),zero_t.temperature(rr,pp,zz,.625),zero_t.opd([0.,.001],[0.,-.001],.625)]
    check('zero_source',np.max(abs(fields['zero_source'])),0.,'Steady/transient and optical zero baseline','zero_source')
    # Independent short-time, deep-slab limit; this tests axial load normalization
    # and rho*cp, rather than providing a general finite-cylinder certificate.
    reference_p=replace(p,thickness=reference_thickness_m,h_conv=0.,emissivity=0.)
    reference_prepared=PreparedModel(reference_p,0.,resolution,
        {0:np.array([reference_flux_W_per_m2])},{0:np.array([reference_radial_wavenumber_per_m])})
    config['known_axial_reference']=dict(physical_parameters_SI=asdict(reference_p),
        source_coefficient_W_per_m2=reference_flux_W_per_m2,radial_wavenumber_per_m=reference_radial_wavenumber_per_m,
        num_axial_modes=reference_axial_modes,time_s=reference_time_s)
    reference_solution=reference_prepared.solve_transient(num_axial_modes=reference_axial_modes)
    actual=reference_solution.temperature(0.,0.,0.,reference_time_s)
    expected=2*reference_flux_W_per_m2/reference_p.conductivity*np.sqrt(reference_p.thermal_diffusivity*reference_time_s/np.pi)
    fields['semi_infinite_reference_values_K']=np.array([actual,expected])
    fields['semi_infinite_reference_relative_error']=np.array([(actual-expected)/expected])
    check('known_short_time_axial_reference',abs(fields['semi_infinite_reference_relative_error'][0]),reference_tolerance,
          'Near-uniform single radial mode at r=0, insulated thick slab, short-time constant-flux analytic limit; not a finite-cylinder boundary benchmark',
          'semi_infinite_reference_relative_error')
    # Integrated energy: finite differences of volume temperature; separate surface integrals.
    a,w=leggauss(48);rg=(a+1)*p.radius/2;zg=(a+1)*p.thickness/2
    wr=w*p.radius/2*rg;wz=w*p.thickness/2;ph=np.arange(64)*2*np.pi/64
    from thermal_opd.projection import reconstruct_cosine_field
    q=reconstruct_cosine_field(prepared.coefficients,prepared.roots,rg,ph)
    input_power=2*np.pi*np.dot(wr,q.mean(axis=1))
    balance=[]
    for t in times_s:
        Tplus=transient.temperature(rg[:,None,None],ph[None,:,None],zg[None,None,:],t+time_step_s)
        Tminus=transient.temperature(rg[:,None,None],ph[None,:,None],zg[None,None,:],t-time_step_s)
        storage=p.density*p.specific_heat*2*np.pi*np.einsum('rz,r,z->',((Tplus-Tminus)/(2*time_step_s)).mean(axis=1),wr,wz)
        faces=transient.temperature(rg[:,None],ph[None,:],0.,t)+transient.temperature(rg[:,None],ph[None,:],p.thickness,t)
        side=transient.temperature(p.radius,ph[None,:],zg[:,None],t)
        loss=p.h_eff*2*np.pi*(np.dot(wr,faces.mean(axis=1))+p.radius*np.dot(wz,side.mean(axis=1)))
        balance.append((storage+loss-input_power)/max(abs(input_power),1e-30))
    fields['transient_energy_normalized_residual']=np.array(balance)
    check('transient_energy_balance',max(abs(fields['transient_energy_normalized_residual'])),energy_tolerance,
          'Independent finite-time volume derivative and separate three surface integrals','transient_energy_normalized_residual')
    energy=steady_energy_balance(prepared.coefficients,prepared.roots,p.radius,p.conductivity,p.thickness,p.h_eff,p.h_eff,p.h_eff,num_r=48,num_z=48,num_phi=64)
    fields['steady_energy_normalized_residual']=np.array([energy['balance_absolute_error_W']/max(abs(energy['modal_input_W']),1e-30)])
    check('steady_energy_balance',fields['steady_energy_normalized_residual'][0],energy_tolerance,'Independent surface quadrature','steady_energy_normalized_residual')
    projected=[]
    for factor in (1,2,4):
        refined=ThermalModel(p,args.angle,replace(resolution,projection_r=resolution.projection_r*factor,
                                                projection_phi=resolution.projection_phi*factor)).prepare()
        projected.append(reconstruct_cosine_field(refined.coefficients,refined.roots,r,phi))
    fields['projection_grid_normalized_changes']=(np.array(projected[:-1])-projected[-1])/max(np.max(abs(projected[-1])),1e-30)
    check('projection_quadrature_convergence',np.max(abs(fields['projection_grid_normalized_changes'])),limit_tolerance,
          'Fixed M/N/roots; projection radial/angular grids multiplied by 1,2,4; projected-flux max error over reference peak',
          'projection_grid_normalized_changes')
    early=np.stack([prepared.solve_transient(num_axial_modes=J).temperature(r[:,None],phi[None,:],0.,early_time_s) for J in early_axial_counts])
    fields['early_front_temperature_K']=early
    early_scale=max(np.max(abs(early[-1])),np.finfo(float).tiny)
    errors=np.max(abs(early[:-1]-early[-1]),axis=(1,2))/early_scale
    fields['early_selected_J_error']=np.array([errors[early_axial_counts.index(num_axial_modes)]])
    fields['early_reference_stability']=np.array([errors[-1]])
    check('initial_neighborhood_selected_J',fields['early_selected_J_error'][0],limit_tolerance,
          f'Explicit additional scope t={early_time_s}s, J={num_axial_modes}; failure does not invalidate previously selected t>=0.625s','early_selected_J_error')
    check('initial_neighborhood_reference_stability',errors[-1],limit_tolerance,
          f'J{early_axial_counts[-2]} vs J{early_axial_counts[-1]} at t={early_time_s}s','early_reference_stability')
    for name in ('joint_modal_truncation','fem_independent_verification','cross_method_comparison'):
        checks.append(verification_check(name,'not_run',scope='Separate experiment/evidence required; this run does not certify it'))
    report=verification_report(checks,config=config,run_directory=directory)
    np.savez_compressed(data_directory(directory)/'verification_fields.npz',**fields)
    (directory/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    write_summary(directory, dict(kind='scientific_verification', status='COMPLETE', verification_type='semi_analytical', physical_parameters_SI=asdict(p), config=config, verification_status=report['status'], provenance=provenance(repository_root, __file__)), defaults=script_defaults(globals()))
    print(f"COMPLETE: {directory}; scientific status={report['status']}; counts={report['counts']}")


if __name__=='__main__':
    main()

"""Independent uniform-front-flux benchmark for prism FEM, not SA comparison."""
from thermal_opd.run_io import script_defaults
from thermal_opd.run_io import data_directory
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from thermal_opd.fem import solve_steady_fem,solve_transient_fem
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.transient import PowerHistory
from thermal_opd.run_io import create_run,provenance,write_summary
from thermal_opd.transient_validation import transient_uniform_flux_1d_reference
from thermal_opd.run_io import validation_check,validation_report

parameters=PhysicalParameters()
uniform_flux_W_per_m2=1000.
reference_intervals=(256,512)
axial_layers=(8,16,32)
time_steps_s=(.25,.125,.0625)
times_s=(0.,.625,2.5,20.)
reference_tolerance,local_fem_tolerance=1e-4,2e-2
repository_root=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    directory=create_run(repository_root,'fem_verification',args.output_dir)
    p=parameters;q=uniform_flux_W_per_m2
    config=dict(physical_parameters_SI=asdict(p),front_flux_W_per_m2=q,h_front=p.h_eff,h_back=p.h_eff,h_side=0.,
        reference_intervals=list(reference_intervals),axial_layers=list(axial_layers),time_steps_s=list(time_steps_s),times_s=list(times_s),
        reference_tolerance=reference_tolerance,local_fem_tolerance=local_fem_tolerance)
    fields={};checks=[]
    def check(name,value,threshold,scope):
        fields[name]=np.array([value])
        checks.append(validation_check(name,'passed' if value<=threshold else 'failed',metrics={'max_error':float(value)},threshold=threshold,scope=scope,evidence=name))
    refs=[]
    for n in reference_intervals:
        z,T=transient_uniform_flux_1d_reference(p.thickness,p.conductivity,p.density*p.specific_heat,q,p.h_eff,p.h_eff,times_s,num_intervals=n)
        refs.append((z,T));fields[f'reference_{n}_z_m']=z;fields[f'reference_{n}_T_K']=T
    z,T=refs[-1]
    previous=np.stack([np.interp(z,refs[0][0],row) for row in refs[0][1]])
    check('independent_reference_stability',np.max(abs(previous[1:]-T[1:]))/np.max(abs(T[1:])),reference_tolerance,'1D finite-volume space refinement; exact modal time response')
    thermal=(p.radius,p.thickness,p.conductivity,p.h_eff,p.h_eff,0.,lambda x,y:q+0*x)
    steady=solve_steady_fem(*thermal,disk_refinements=1,axial_layers=axial_layers[-1],quadrature_order=4,solver='direct',rtol=1e-10)
    expected=q*(p.conductivity+p.h_eff*(p.thickness-z))/(p.conductivity*(2*p.h_eff)+p.h_eff**2*p.thickness)
    check('steady_known_linear_solution',np.max(abs(steady.temperature_at(0*z,0*z,z)-expected))/np.max(abs(expected)),1e-10,'Exact uniform-flux insulated-side linear axial solution')
    configurations=[(n,time_steps_s[0]) for n in axial_layers]+[(axial_layers[-1],dt) for dt in time_steps_s[1:]]
    errors=[];energy=[];linear=[]
    for i,(n,dt) in enumerate(configurations):
        result=solve_transient_fem(*thermal,density=p.density,specific_heat=p.specific_heat,times=times_s,
            time_step=dt,history=PowerHistory(),disk_refinements=1,axial_layers=n,quadrature_order=4,solver='direct',rtol=1e-10)
        sampled=np.stack([f.temperature_at(0*z,0*z,z) for f in result.frames])
        fields[f'configuration_{i}_T_K']=sampled
        errors.append(float(np.max(np.max(abs(sampled[1:]-T[1:]),axis=1)/np.max(abs(T[1:]),axis=1))))
        energy.append(max(v['balance_absolute_error_W'] for v in result.steps)/result.diagnostics['full_input_W'])
        linear.append(max(v['relative_linear_residual'] for v in result.steps))
        check(f'initial_condition_{i}',np.max(abs(result.frames[0].temperature_K)),0.,'All prism FEM nodes')
    fields['configuration_errors']=np.array(errors)
    fields['configurations']=np.array(configurations)
    check('transient_known_1d_reference',errors[-1],local_fem_tolerance,'Local benchmark only; not final Gaussian mesh/time certificate')
    for name,indices in (('space_refinement',range(len(axial_layers))),('time_refinement',range(len(axial_layers)-1,len(configurations)))):
        values=np.array([errors[i] for i in indices])
        check(name,max(0.,np.max(np.diff(values))),1e-12,'Error should decrease under refinement against independently refined 1D reference')
    check('energy_balance',max(energy),1e-9,'Per-step storage plus surface losses equals applied input')
    check('linear_equation_residual',max(linear),1e-9,'Discrete weak-form equations, including natural boundary loads')
    for name in ('gaussian_case_mesh_convergence','gaussian_case_time_convergence','pointwise_boundary_flux','cross_method_comparison'):
        checks.append(validation_check(name,'not_run',scope='Separate general-case verification; weak-form/energy checks do not prove pointwise flux'))
    report=validation_report(checks,config=config,run_directory=directory)
    np.savez_compressed(data_directory(directory)/'verification_fields.npz',**fields)
    (directory/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    write_summary(directory, dict(kind='scientific_verification', status='COMPLETE', verification_type='fem', config=config, physical_parameters_SI=asdict(p), verification_status=report['status'], provenance=provenance(repository_root, __file__)), defaults=script_defaults(globals()))
    print(f"COMPLETE: {directory}; scientific status={report['status']}; counts={report['counts']}")


if __name__=='__main__':
    main()

"""Scientific evidence records and independent reference/checking paths.

These functions are opt-in; ordinary model preparation never calls them.
No claim of convergence is inferred from successful software execution.
"""
from thermal_opd.run_io import data_directory
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import eigh_tridiagonal


def verification_check(name,status,*,metrics=None,threshold=None,scope=None,evidence=None):
    if status not in ('passed','failed','not_run'):
        raise ValueError('Scientific status must be passed, failed or not_run')
    return dict(name=name,status=status,metrics=metrics or {},threshold=threshold,scope=scope,evidence=evidence)


def verification_report(checks,*,config,run_directory):
    names=[row['name'] for row in checks]
    if len(set(names))!=len(names):
        raise ValueError('Duplicate verification check names')
    statuses=[row['status'] for row in checks]
    status='failed' if 'failed' in statuses else ('not_run' if not checks or 'not_run' in statuses else 'passed')
    return dict(schema_version=1,status=status,checked_at_utc=datetime.now(timezone.utc).isoformat(),
                run_directory=str(run_directory),config=config,checks=checks,
                counts={key:statuses.count(key) for key in ('passed','failed','not_run')})


def finite_difference_pde_residual(temperature,p,r,phi,z,time,*,space_step,time_step):
    """Independent central differences of rho cp T_t-k Laplacian(T), W/m³.

    Interior r>0 samples only; phi step is chosen to give the same arc length
    as the radial/axial displacement. Uses field values, no solver derivatives.
    Caller studies finite-difference step sensitivity separately.
    """
    r,phi,z=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(phi,dtype=float),np.asarray(z,dtype=float))
    h,dt=space_step,time_step
    if not all(np.isfinite(v) and v>0 for v in (h,dt)) or time<=dt:
        raise ValueError('Require positive finite steps and t>dt')
    if np.any(r<=h) or np.any(r+h>=p.radius) or np.any(z<=h) or np.any(z+h>=p.thickness):
        raise ValueError('PDE stencil must stay inside the cylinder away from r=0')
    center=temperature(r,phi,z,time)
    plus,minus=temperature(r+h,phi,z,time),temperature(r-h,phi,z,time)
    radial=(plus-2*center+minus)/h**2+(plus-minus)/(2*h*r)
    dphi=h/r
    angular=(temperature(r,phi+dphi,z,time)-2*center+temperature(r,phi-dphi,z,time))/(dphi**2*r**2)
    axial=(temperature(r,phi,z+h,time)-2*center+temperature(r,phi,z-h,time))/h**2
    temporal=(temperature(r,phi,z,time+dt)-temperature(r,phi,z,time-dt))/(2*dt)
    return p.density*p.specific_heat*temporal-p.conductivity*(radial+angular+axial)


def uniform_flux_1d_reference(thickness,conductivity,volumetric_capacity,flux,h_front,h_back,times,*,num_intervals):
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


def validate_scientific_artifacts(directory,summary):
    """Validate evidence/flags without turning a scientifically failed run green."""
    import json
    from pathlib import Path
    path=Path(directory)
    report=json.loads((path/'verification.json').read_text(encoding='utf-8'))
    with np.load(data_directory(path)/'verification_fields.npz',allow_pickle=False) as f:
        if not all(np.all(np.isfinite(f[k])) for k in f.files):
            raise ValueError('Nonfinite scientific verification evidence')
        for row in report['checks']:
            verification_check(row['name'],row['status'])
            if row['status']=='not_run':
                continue
            key=row['evidence']
            actual=float(np.max(np.abs(f[key])))
            threshold=row['threshold']
            expected='passed' if actual<=threshold else 'failed'
            if not np.isfinite(threshold) or threshold<0 or row['status']!=expected or not np.isclose(actual,row['metrics']['max_error'],rtol=1e-12,atol=1e-15):
                raise ValueError('Scientific verification metric/status differs from saved evidence')
    recomputed=verification_report(report['checks'],config=report['config'],run_directory=path)
    if report['status']!=recomputed['status'] or report['counts']!=recomputed['counts'] or summary['verification_status']!=report['status']:
        raise ValueError('Scientific verification aggregate status differs')
    return report


def verification_from_validated_run(directory,summary,boundary_check=None):
    """Map executed saved-run checks to scoped scientific evidence records.

    Called only after run_validation has recomputed their underlying metrics.
    Unknown/unperformed checks remain not_run even if file integrity passes.
    """
    kind=summary.get('kind')
    checks=[]
    def add(name,executed,scope,evidence=None,metrics=None,threshold=None):
        checks.append(verification_check(name,'passed' if executed else 'not_run',scope=scope,evidence=evidence,
                                         metrics=metrics,threshold=threshold))
    selection=kind in ('steady_joint_truncation_selection','transient_joint_truncation_selection')
    add('joint_modal_truncation',selection or summary.get('transient_truncation_established',False),
        'Only the recorded angles, times, history and finite-reference tolerances',summary.get('selection_directory',str(directory)))
    add('projection_quadrature_convergence',False,'Joint modal/projection refinement does not isolate projection quadrature')
    add('pde_residual',False,'Independent PDE experiment/evidence required')
    transient=kind in ('transient_fine_fields','transient_modal_fem_study')
    add('initial_condition',transient,'Saved zero-time semi-analytical field; FEM nodes additionally checked in FEM studies')
    if transient:
        checks[-1].update(metrics={'max_initial_temperature_K':0.,'max_initial_opd_m':0.},threshold=0.)
    add('initial_neighborhood',False,'Exact t=0 assignment does not establish t→0+ convergence')
    long_errors=summary.get('long_time_errors_vs_same_modes_steady',summary.get('long_time_errors_vs_steady'))
    add('long_time_steady_limit',long_errors is not None,'Recorded long-time comparison with same retained steady modes',
        'summary.json',metrics=long_errors,threshold=summary.get('precision_target',1e-3))
    if long_errors is not None:
        from thermal_opd.convergence import within_tolerance
        checks[-1]['status']='passed' if within_tolerance(long_errors,checks[-1]['threshold']) else 'failed'
    steady=kind=='complete_representative_steady_refracted_ray_case'
    add('thermal_boundaries',steady or (boundary_check is not None and boundary_check['status']=='PASS'),
        'Sampled semi-analytical front/rear/side residuals against retained heat flux',
        'summary.json' if steady else 'boundary_summary.json',
        metrics=summary.get('boundary_absolute_residuals_W_per_m2') if steady else boundary_check,
        threshold=1e-7 if steady else (boundary_check or {}).get('allowed_absolute_residual_W_per_m2'))
    add('steady_energy_balance',steady,'Independent steady surface integrals; only the retained source',
        'summary.json',metrics=summary.get('energy_balance'),threshold='1e-9 W + 1e-7*abs(modal_input_W)')
    fem=kind in ('transient_modal_fem_study','steady_prismatic_fem_comparison')
    add('fem_energy_balance',fem,'Reloaded nodal storage/loss integrals and recorded input')
    add('fem_weak_equation_residual',fem,'Discrete weak-form linear residual; not pointwise normal flux')
    add('fem_independent_reference',False,'Independent FEM reference experiment/evidence required')
    add('fem_space_convergence',bool(summary.get('mesh_convergence_established',False)),
        'Only recorded mesh-refinement target; smoke is not a final certificate')
    add('fem_time_convergence',False,'Separate transient time-refinement acceptance required')
    add('cross_method_comparison',False,'Fields/metrics may be available; agreement is not an independent verification certificate')
    return verification_report(checks,config={'physical_parameters_SI':summary.get('physical_parameters_SI'),
        'kind':kind,'scope':summary.get('precision_scope')},run_directory=directory)

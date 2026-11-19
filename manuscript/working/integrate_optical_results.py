"""Integrate only completed, checked optical-update results into paper assets.

Run after recompute_paper_optics stages. The normal-incidence relative FEM
metrics are inherited by an exact positive scaling (unchanged straight path
and pupil); oblique optical metrics come from new FEM node/ray evaluations.
Historical temperature-only refinement is retained because the thermal PDE,
load, boundaries and thermal parameters are unchanged. It is not relabeled
as a new FEM optical-refinement certificate. Added 2026-11-03.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
from thermal_opd.parameters import PhysicalParameters
from dataclasses import asdict

ROOT=Path(__file__).resolve().parents[2]
PAPER=ROOT/'manuscript'
AUDIT=ROOT/'outputs/experiment_audit/20261017T121521368371Z_f107298c'
INDEX=json.loads((AUDIT/'summary.json').read_text())
STAGES={row['name']:row for row in INDEX['stages']}
INPUTS={}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def historical(stage,name='summary.json'):
    key=STAGES[stage]['run_directory']+'/'+name
    path=AUDIT/'data/snapshot'/key
    assert sha(path)==INDEX['files'][key]
    INPUTS[path.relative_to(ROOT).as_posix()]=sha(path)
    return json.loads(path.read_text())


def peak(errors,key):return max(float(v)for e in errors for v in np.atleast_1d(e[key]))


def update_macro(text,name,transform):
    pattern=r'\\newcommand\{\\'+name+r'\}\{.*?(?=\n\\newcommand|\Z)'
    return re.sub(pattern,lambda m:transform(m[0]),text,flags=re.S)


def rows_replace(text,rows):
    for angle,row in rows.items():
        pattern=r'^\$'+str(angle)+r'\^\\circ\$ & .*?\\\\(?=\\bottomrule|$)'
        text,n=re.subn(pattern,lambda _:row,text,flags=re.M)
        assert n==1,(angle,n)
    return text


def integrate(run):
    p=PhysicalParameters()
    sa={a:json.loads((run/f'sa_theta{a}/summary.json').read_text())for a in [0,30,60]}
    for a in sa:
        assert sa[a]['physical_parameters_SI']==asdict(p)
        assert json.loads((run/f'sa_check_theta{a}.json').read_text())['passed']
    steady_old=historical('steady_fem')
    steady_cases={int(c['angle_deg']):c for c in steady_old['cases']}
    steady={}
    for a in [0,30,60]:
        if a==0:
            e=steady_cases[a]['meshes'][-1]['errors_vs_sa']
            steady[a]=dict(temperature_relative_linf=e['temperature_relative_linf'],optical={
                'pv_relative':[e['opd_relative_pv']],
                'linf_over_reference_pv':[e['opd_piston_removed_linf_over_reference_pv']]})
        else:steady[a]=json.loads((run/f'steady_fem_theta{a}.json').read_text())
    transient={}
    old_zero=historical('transient_fem_theta0')
    transient[0]=dict(times_s=old_zero['times_s'],errors_vs_sa=old_zero['fem_runs'][0]['errors_vs_sa'])
    for a in [30,60]:
        transient[a]=json.loads((run/f'new_opd_fem_theta{a}.json').read_text())
        assert transient[a]['status']=='COMPLETE' and transient[a]['physical_parameters_SI']==asdict(p)
    # Exact zero-angle scaling must also reproduce the newly computed SA PVs.
    old_sa0=historical('transient_fields_theta0')
    ratio=p.dn_dT/old_sa0['physical_parameters_SI']['dn_dT']
    for old in old_sa0['snapshots']:
        t=old['time_s']
        if 0<t<=75:
            new=next(v for v in sa[0]['history']if v['time_s']==t)
            assert np.isclose(new['sampled_pv_nm']*1e-9,old['sampled_pv_m']*ratio,rtol=1e-10,atol=1e-20)
    main_rows={};all_rows=[]
    metrics={}
    for a in [0,30,60]:
        good=[e for t,e in zip(transient[a]['times_s'],transient[a]['errors_vs_sa'])if t>=.0025]
        all_positive=[e for t,e in zip(transient[a]['times_s'],transient[a]['errors_vs_sa'])if t>0]
        s=steady[a]
        row=[100*s['optical']['pv_relative'][0],100*peak(good,'opd_relative_pv'),
             100*peak(good,'opd_piston_removed_linf_over_reference_pv')]
        main_rows[a]=f'${a}^\\circ$ & '+' & '.join(f'{v:.6f}'for v in row)+r'\\'
        metrics[a]=dict(steady_pv_percent=row[0],transient_pv_percent=row[1],transient_linf_percent=row[2])
        all_rows.append(f'Steady & ${a}^\\circ$ & '+ ' & '.join(f'{v:.6f}'for v in [
            100*max(s['temperature_relative_linf']),100*s['optical']['linf_over_reference_pv'][0],row[0]])+r'\\')
        all_rows.append(f'Transient & ${a}^\\circ$ & '+' & '.join(f'{100*peak(all_positive,k):.6f}'for k in [
            'temperature_relative_linf','opd_piston_removed_linf_over_reference_pv','opd_relative_pv'])+r'\\')
    ref=historical('transient_fem_convergence','verification.json')
    checks={c['name']:c for c in ref['checks']}
    temp_rows={}
    for a in [0,30,60]:
        values=[]
        for kind in ['xy','z','time','joint']:
            data=checks[f'theta{a}_{kind}_convergence']['metrics']
            errors=[e for t,e in zip(data['times_s'],data['errors_by_time'])if t>=.0025]
            v=100*max(float(x)for e in errors for key,val in e.items()if key.startswith('temperature_')for x in np.atleast_1d(val))
            assert v <= 100*checks[f'theta{a}_{kind}_convergence']['threshold']
            values.append(v)
        temp_rows[a]=f'${a}^\\circ$ & '+' & '.join(f'{v:.6f}'for v in values)+r'\\'
    tables=(PAPER/'tables.tex').read_text(encoding='utf-8')
    tables=tables.replace(r'$9.6\times10^{-6}$',r'$8.4154\times10^{-6}$').replace('$1.44963$','$1.450076$')
    tables=tables.replace(r'$\gamma_T$ is a retained nominal value tabulated at $643.8$ nm, not a wavelength-specific measurement at $1064$ nm.',
        r'$n_s$ and $\gamma_T$ are evaluated from the same absolute Corning 7980 Sellmeier fit at $1064$ nm and $300$ K and held constant in the linear optical calculation.')
    steady_rows={}
    for a in [0,30,60]:
        with np.load(run/f'sa_theta{a}/fields.npz') as f:temp=f['steady_temperature_K'][:,0,0]
        steady_rows[a]=f'${a}^\\circ$ & '+' & '.join(f'{v:.3f}'for v in [*temp,sa[a]['steady_sampled_pv_nm']])+r'\\'
    tables=update_macro(tables,'SteadyFieldTable',lambda s:rows_replace(s,steady_rows))
    tables=update_macro(tables,'ValidationTable',lambda s:rows_replace(s,main_rows))
    tables=update_macro(tables,'TransientRefinementTable',lambda s:rows_replace(s,temp_rows)
        .replace('across temperature and OPD metrics','across temperature metrics'))
    (PAPER/'tables.tex').write_text(tables,encoding='utf-8')
    # Export just the two Figure2 frames, without modifying the original export.
    export=PAPER/'data/optical_update_20261103'
    export.mkdir(exist_ok=True)
    with np.load(run/'new_opd_fem_theta60.npz') as f:
        ids=[int(np.flatnonzero(f['times_s']==t)[0])for t in [.625,75.]]
        out=export/'transient_60deg_comparison_frames.npz'
        np.savez_compressed(out,opd_x_m=f['opd_x_m'],opd_y_m=f['opd_y_m'],
            opd_physical_mask=f['physical_mask'],opd_fem_common_mask=f['common_mask'],
            times_s=f['times_s'][ids],opd_sa_m=f['opd_sa_physical_m'][ids],opd_fem_m=f['opd_fem_m'][ids])
    export_manifest=dict(export_sha256=sha(out),source_files_unchanged=True,
        physical_parameters_SI=asdict(p),source_file=(run/'new_opd_fem_theta60.npz').relative_to(ROOT).as_posix(),
        source_sha256=sha(run/'new_opd_fem_theta60.npz'),note='New refracted-path evaluations; original exported frames are preserved separately.')
    (export/'transient_60deg_export_manifest.json').write_text(json.dumps(export_manifest,indent=2),encoding='utf-8')
    # Supplement: retain five sections, make the refinement scope explicit.
    sup=(PAPER/'supplement.tex').read_text(encoding='utf-8')
    prefix='The steady temperature reference uses' if 'The steady temperature reference uses' in sup else 'The steady reference uses'
    a=sup.index(prefix);b=sup.index('\nThe positive observation times',a)
    sup=sup[:a]+r'''The steady temperature reference uses $M=14$, $N=32$, a $2400\times256$ source grid, and 64 ray points; the enlarged temperature reference uses $M=16$, $N=40$, a $3200\times512$ source grid, and 80 ray points.
The transient temperature references use the same two transverse configurations with $J=512$ and $768$, respectively.
The selected steady configuration is $M=6$, $N=8$, $400\times32$, and 16 ray points; the selected transient configuration is $M=10$, $N=20$, $J=256$, $1200\times128$, and 32 ray points.
The finite-reference target is $0.1\%$, with additional transient temperature-gradient and near-surface depth checks.
For optical evaluation, the steady reference uses $M=10$, $N=20$, $1200\times128$, and 32 ray points.
The transient optical reference uses $M=14$, $N=28$, $J=384$, $1800\times256$, and 64 ray points; a separate check doubles the selected transient ray quadrature from 32 to 64 points.
All optical evaluations use the matched $1064$ nm, $300$ K inputs in the main-text parameter table.
''' +sup[b:]
    sup=sup.replace('Figure~\\ref{fig:supp-fem-refinement} shows the recorded refinement changes at each positive observation time.',
        'Figure~\\ref{fig:supp-fem-refinement} shows the temperature-refinement changes at each positive observation time.')
    sup=sup.replace('Transient FEM refinement changes for','Transient FEM temperature-refinement changes for')
    sup=sup.replace('takes the maximum across the recorded temperature and OPD metrics at a given time.',
        'takes the maximum across the recorded temperature metrics at a given time.')
    sup=sup.replace('Optical differences use the common-pupil PV of the corresponding refined FEM as their denominator, with piston removed for optical field differences.\n','')
    sup=sup.replace('These curves compare FEM discretizations; the FEM--semi-analytical differences are reported separately in Table~\\ref{tab:supp-comparison}.',
        'These curves assess the thermal FEM discretization; optical method-to-method differences on the refracted paths are reported separately in Table~\\ref{tab:supp-comparison}.')
    sup=sup.replace('The FEM refinement comparisons use a $201\\times201$ Cartesian thermal grid and a $101\\times101$ optical grid on their respective common masks.',
        'The thermal FEM refinement comparisons use a $201\\times201$ Cartesian thermal grid on the common mesh domain.')
    continuous={a:{k:v['pv']*1e9 for k,v in sa[a]['continuous_extrema'].items()}for a in sa}
    for marker,key in [('For $0^\\circ$, $30^\\circ$, and $60^\\circ$, the saved continuous steady PV values are','steady')]:
        pattern=re.escape(marker)+r'.*? nm, respectively\.'
        line=marker+' '+', '.join(f'${continuous[a][key]:.3f}$'for a in [0,30])+f', and ${continuous[60][key]:.3f}$ nm, respectively.'
        sup,n=re.subn(pattern,lambda _:line,sup);assert n==1
    pattern=r'For the transient fields, continuous-pupil OPD-PV estimates at .*? nm\.'
    line='For the transient fields, continuous-pupil OPD-PV estimates at $0.625$ s are '+', '.join(f'${continuous[a]["0.625"]:.3f}$'for a in [0,30])+f', and ${continuous[60]["0.625"]:.3f}$ nm for $0^\\circ$, $30^\\circ$, and $60^\\circ$; at $75$ s they are '+', '.join(f'${continuous[a]["75"]:.3f}$'for a in [0,30])+f', and ${continuous[60]["75"]:.3f}$ nm.'
    sup,n=re.subn(pattern,lambda _:line,sup);assert n==1
    a=sup.index('Steady & $0^\\circ$');b=sup.index('\\bottomrule',a)
    ordered=[all_rows[i]for i in [0,2,4,1,3,5]]
    sup=sup[:a]+'\n'.join(ordered)+sup[b:]
    (PAPER/'supplement.tex').write_text(sup,encoding='utf-8')
    # Keep exact full numbers for checks, manuscript tables use stated rounding.
    write_data=dict(physical_parameters_SI=asdict(p),comparison_metrics_percent=metrics,
        continuous_pv_nm=continuous,normal_incidence_scaling=ratio,
        normal_incidence_derivation='At theta0 the refracted angle and pupil are unchanged; positive gamma scales both optical solutions and their reference PV equally. Relative optical discrepancies are invariant.',
        temperature_refinement_source=STAGES['transient_fem_convergence']['run_directory'],
        normal_incidence_cross_method_source=STAGES['transient_fem_theta0']['run_directory'])
    (run/'paper_values.json').write_text(json.dumps(write_data,indent=2),encoding='utf-8')
    manifest=dict(run_directory=run.relative_to(ROOT).as_posix(),physical_parameters_SI=asdict(p),
        optical_model=dict(authors=['D. B. Leviton','B. J. Frey'],year=2006,doi='10.1117/12.672853',table='3',
            interpretation='Absolute n and its analytic temperature derivative at the same reference wavelength and temperature; constant linearization thereafter.'),
        scope=dict(historical_thermal_refinement='unchanged thermal problem; temperature metrics only',
            current_optical_checks='finite-reference SA checks plus independently assembled FEM method-to-method comparison',
            new_fem_optical_refinement_certificate='not claimed',zero_angle='exact positive scaling of historical optical fields'),
        historical_input_sha256=INPUTS,new_input_sha256={str(f.relative_to(ROOT)).replace('\\','/'):sha(f)
            for f in run.rglob('*.json') if 'paper_values'not in f.name})
    (PAPER/'optical_revision_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(metrics,indent=2));print('PAPER_VALUES_INTEGRATED')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-directory',required=True,type=Path)
    integrate(p.parse_args().run_directory.resolve())


if __name__=='__main__':main()

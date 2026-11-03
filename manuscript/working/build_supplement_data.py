"""Export already computed supplementary data; no solver or figure changes.

Run when rebuilding supplement.tex or supplement_zh.tex. Reuses the checked
historical reader in integrate_optical_results.py and the optical manifest.
Writes shared numeric TeX rows, CSV tables and provenance under manuscript/data/
supplement_details. A hash mismatch, missing/duplicate time or invalid metric
raises an exception rather than silently exporting an incomplete table.
Controls cover percent conversion, the maximum over thermal planes and rejection
of negative/nonfinite/misaligned input. Independently check exported CSV values
against the named JSON metric arrays before trusting the rendered tables.
Added 2026-11-03 for the data-focused supplement and its Chinese version.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
from integrate_optical_results import historical, INPUTS

ROOT=Path(__file__).resolve().parents[2]
PAPER=ROOT/'manuscript'
OUT=PAPER/'data/supplement_details'
MANIFEST=json.loads((PAPER/'optical_revision_manifest.json').read_text(encoding='utf-8'))
RUN=ROOT/MANIFEST['run_directory']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def metric_row(time,error):
    if not math.isfinite(time) or time<=0:
        raise ValueError('Expected a positive observation time')
    keys=('opd_piston_removed_rms_over_reference_pv',
          'opd_piston_removed_linf_over_reference_pv','opd_relative_pv')
    vals=[max(error['temperature_relative_linf'])]+[error[k]for k in keys]
    if not all(math.isfinite(v) and v>=0 for v in vals):
        raise ValueError('Invalid relative discrepancy')
    return [time]+[100*v for v in vals]

def time_rows(times,errors):
    if len(times)!=len(errors) or len(set(times))!=len(times):
        raise ValueError('Mismatched or duplicate observation times')
    return [metric_row(t,e)for t,e in zip(times,errors)if t>0]

def controls():
    e={'temperature_relative_linf':[.001,.003,.002],
       'opd_piston_removed_rms_over_reference_pv':.004,
       'opd_piston_removed_linf_over_reference_pv':.005,'opd_relative_pv':.006}
    assert metric_row(2,e)==[2,.3,.4,.5,.6]
    for times,errors in [([1,1],[e,e]),([1,2],[e])]:
        try:time_rows(times,errors)
        except ValueError:pass
        else:raise AssertionError('Corrupted time record was accepted')
    for bad in (-1,float('nan')):
        corrupted=dict(e,opd_relative_pv=bad)
        try:metric_row(1,corrupted)
        except ValueError:pass
        else:raise AssertionError('Corrupted metric was accepted')

def checked(path):
    name=path.relative_to(ROOT).as_posix()
    group=MANIFEST['new_input_sha256']
    assert name in group and sha(path)==group[name],name
    INPUTS[name]=sha(path)
    return json.loads(path.read_text(encoding='utf-8'))

def output(stem,columns,rows,tex_rows):
    with (OUT/(stem+'.csv')).open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f);writer.writerow(columns);writer.writerows(rows)
    if stem in ('sa_checks','pupil_pv'):
        # Keep booktabs' noalign immediately after the last row inside the input.
        tex_rows[-1]+=r'\bottomrule'
    (OUT/(stem+'.tex')).write_text('\n'.join(tex_rows)+'\n',encoding='utf-8')

def main():
    controls()
    OUT.mkdir(exist_ok=True,parents=True)
    zero=historical('transient_fem_theta0')
    runs={0:dict(times_s=zero['times_s'],errors_vs_sa=zero['fem_runs'][0]['errors_vs_sa'])}
    for a in (30,60):runs[a]=checked(RUN/f'new_opd_fem_theta{a}.json')
    steady_source=historical('steady_fem')
    steady_cases={int(c['angle_deg']):c for c in steady_source['cases']}
    rows=[];tex=[]
    for a,run in runs.items():
        if a==0:
            error=steady_cases[a]['meshes'][-1]['errors_vs_sa']
        else:
            d=checked(RUN/f'steady_fem_theta{a}.json')
            error=dict(temperature_relative_linf=d['temperature_relative_linf'],
                opd_piston_removed_rms_over_reference_pv=d['optical']['rms_over_reference_pv'][0],
                opd_piston_removed_linf_over_reference_pv=d['optical']['linf_over_reference_pv'][0],
                opd_relative_pv=d['optical']['pv_relative'][0])
        sv=metric_row(1,error)[1:]
        rows.append([a,'steady']+sv)
        tex.append(f'${a}^\\circ$ & '+r'\SteadyState & '+' & '.join(f'{v:.6f}'for v in sv)+r'\\')
        values=time_rows(run['times_s'],run['errors_vs_sa'])
        assert [v[0]for v in values]==[.001,.0025,.01,.05,.1,.625,2.5,20.,50.,75.]
        for value in values:
            rows.append([a]+value)
            tex.append(f'${a}^\\circ$ & {value[0]:g} & '+' & '.join(f'{v:.6f}'for v in value[1:])+r'\\')
        if a!=60:tex.append(r'\midrule')
    output('time_comparison',['angle_deg','time_s','temperature_linf_percent',
        'opd_rms_over_sa_pv_percent','opd_linf_over_sa_pv_percent','opd_pv_percent'],rows,tex)

    rows=[];tex=[]
    checks=[('steady_reference',r'\SteadyCheck'),('modal_reference',r'\TransientCheck'),('ray_quadrature',r'\RayCheck')]
    sa={a:checked(RUN/f'sa_check_theta{a}.json')for a in (0,30,60)}
    for key,label in checks:
        for a in (0,30,60):
            d=sa[a]['checks'][key]
            values=[100*max(d[k])for k in ('pv_relative','rms_over_reference_pv','linf_over_reference_pv')]
            rows.append([key,a]+values)
            tex.append(label+f' & ${a}^\\circ$ & '+' & '.join(f'{v:.6f}'for v in values)+r'\\')
        if key!='ray_quadrature':tex.append(r'\midrule')
    output('sa_checks',['comparison','angle_deg','pv_percent','rms_over_ref_pv_percent','linf_over_ref_pv_percent'],rows,tex)

    rows=[];tex=[]
    summaries={a:checked(RUN/f'sa_theta{a}/summary.json')for a in (0,30,60)}
    for t,label in [('steady',r'\SteadyState'),('0.625','0.625 s'),('75','75 s')]:
        for a,d in summaries.items():
            sampled=d['steady_sampled_pv_nm'] if t=='steady' else next(v['sampled_pv_nm']for v in d['history']if v['time_s']==float(t))
            cont=d['continuous_extrema'][t]['pv']*1e9
            delta=cont-sampled
            assert delta>=-1e-9
            rows.append([t,a,sampled,cont,delta])
            tex.append(label+f' & ${a}^\\circ$ & {sampled:.4f} & {cont:.4f} & {delta:.4f}'+r'\\')
        if t!='75':tex.append(r'\midrule')
    output('pupil_pv',['state_or_time_s','angle_deg','sampled_pv_nm','continuous_pv_nm','difference_nm'],rows,tex)

    # Save all four mesh configurations directly from the historical records.
    steady=steady_source['cases'][0]['meshes']
    conv=historical('transient_fem_convergence')
    verification=historical('transient_fem_convergence','verification.json')
    (OUT/'fem_configuration_source.json').write_text(json.dumps({
        'steady':[{'xy_dofs':d['xy_dofs'],'axial_elements':d['axial_layers'],
                   'total_dofs':d['degrees_of_freedom']}for d in steady],
        'transient_config':conv['config']},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    assert verification['status']=='passed'
    (OUT/'provenance.json').write_text(json.dumps({
        'source_files_sha256':INPUTS,
        'scope':'Retrieval and percent conversion only; no new numerical solution',
        'normal_incidence':'Relative discrepancies invariant under common positive thermo-optic scaling; straight paths and pupil unchanged',
        'table_time_comparison':'Temperature maximum over three recorded planes; full cross-method norms as defined in saved results',
        'sa_checks':'Steady one case; transient maxima over the ten positive observation times',
        'controls':'Known percent conversion passed; duplicate/truncated times, negative and NaN metrics rejected'
    },ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Exported 33 method-comparison rows, 9 SA-check rows and 9 pupil-PV rows; input hashes checked.')

if __name__=='__main__':main()

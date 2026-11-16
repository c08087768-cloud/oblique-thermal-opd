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
    if stem in ('sa_checks','pupil_pv','fem_optical_steady','fem_optical_transient'):
        # Keep booktabs' noalign immediately after the last row inside the input.
        tex_rows[-1]+=r'\bottomrule'
    (OUT/(stem+'.tex')).write_text('\n'.join(tex_rows)+'\n',encoding='utf-8')

def main():
    global historical, INPUTS
    from integrate_optical_results import historical, INPUTS
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

# 2026-11-17: separate export of completed FEM optical refinements.  The
# --optical-refinement branch never rewrites the historical tables above.
# This is a retrieval/aggregation tool, not a solver or scientific certificate.
OPTICAL_KEYS=('opd_relative_pv','opd_piston_removed_rms_over_reference_pv',
              'opd_piston_removed_linf_over_reference_pv')
OPTICAL_TIMES=(.001,.0025,.01,.05,.1,.625,2.5,20.,50.,75.)
REPLAY=ROOT/'outputs/optical_server_replay/optics_20261109T220402Z_eb35d82'
COMPLETION=ROOT/'outputs/optical_completion/20261116T205530Z_ab7d2d6'

def optical_percent(metrics):
    values=[metrics[k] for k in OPTICAL_KEYS]
    if not all(isinstance(v,(int,float)) and math.isfinite(v) and v>=0 for v in values):
        raise ValueError('Invalid optical refinement metric')
    return [100*v for v in values]

def transient_optical_rows(records,direction):
    rows=[r for r in records if r['direction']==direction]
    if tuple(r['time_s']for r in rows)!=OPTICAL_TIMES:
        raise ValueError('Missing, duplicate or reordered optical time')
    threshold=.001 if direction=='time' else .005
    for r in rows:
        gate=r['time_s']!=.001
        values=optical_percent(r['metrics'])
        if r['gate']!=gate or r['threshold']!=threshold:
            raise ValueError('Incorrect optical refinement scope or threshold')
        expected='PASS' if gate else 'NOT APPLICABLE'
        if r['status']!=expected or (gate and max(values)>100*threshold):
            raise ValueError('Optical refinement target not met')
    maxima=[max(optical_percent(r['metrics'])[i]for r in rows[1:])for i in range(3)]
    return rows,maxima,optical_percent(rows[0]['metrics'])

def optical_controls():
    # Written before extracting the real records: .001 fraction = .1 percent;
    # a deliberately larger 1 ms diagnostic must never enter the gated maximum.
    from copy import deepcopy
    plant=dict(zip(OPTICAL_KEYS,(.001,.002,.003)))
    assert optical_percent(plant)==[.1,.2,.3]
    assert optical_percent(dict.fromkeys(OPTICAL_KEYS,0))==[0,0,0]
    records=[dict(direction='xy',time_s=t,gate=t!=.001,threshold=.005,
                  metrics=dict(plant),status='PASS' if t!=.001 else 'NOT APPLICABLE')
             for t in OPTICAL_TIMES]
    records[0]['metrics']=dict.fromkeys(OPTICAL_KEYS,.2)
    assert transient_optical_rows(records,'xy')[1:]==([.1,.2,.3],[20.,20.,20.])
    corruptions=[]
    for value in (-1,float('nan'),float('inf')):
        bad=deepcopy(records);bad[1]['metrics'][OPTICAL_KEYS[0]]=value;corruptions.append(bad)
    bad=deepcopy(records);bad[1]['time_s']=.001;corruptions.append(bad)
    corruptions.append(deepcopy(records[:-1]))
    bad=deepcopy(records);bad[0]['gate']=True;corruptions.append(bad)
    bad=deepcopy(records);bad[1]['threshold']=.01;corruptions.append(bad)
    bad=deepcopy(records);bad[1]['metrics'][OPTICAL_KEYS[0]]=.006;corruptions.append(bad)
    for bad in corruptions:
        try:transient_optical_rows(bad,'xy')
        except ValueError:pass
        else:raise AssertionError('Corrupted optical record was accepted')
    return dict(percent_and_null='passed',diagnostic_exclusion='passed',
                corruptions_rejected=len(corruptions))

def optical_display(value):
    # A display floor avoids attaching physical significance to roundoff-sized
    # changes. Full saved precision is retained in CSV and the evidence JSON.
    return r'$<10^{-5}$' if value<1e-5 else f'{value:.5f}'

def export_optical_refinement():
    controls_result=optical_controls()
    sources={}
    def read(path):
        name=path.relative_to(ROOT).as_posix()
        sources[name]=sha(path)
        return json.loads(path.read_text(encoding='utf-8'))
    OUT.mkdir(exist_ok=True,parents=True)
    evidence={'steady':[],'transient':[],'steady_load_quadrature':[]}
    steady=[];steady_tex=[]
    for group in ('steady_0_30','steady_60'):
        summary=read(REPLAY/group/'summary.json');config=read(REPLAY/group/'config.json')
        if summary['status']!='PASS' or summary['thermal_solves']!=0 or config['mesh_tolerance']!=.001:
            raise ValueError('Steady replay scope mismatch')
        for case in summary['cases']:
            angle=case['angle_deg']
            if case['old_fields_replay']!='PASS' or case['saved_new_metric_recheck']!='PASS':
                raise ValueError('Steady optical replay/readback did not pass')
            for direction in ('xy','z','joint'):
                check=case['new_directional_checks'][direction]
                values=optical_percent(check['errors'])
                if not (check['status']=='PASS' and check['configuration_valid'] and check['passed']
                        and check['candidate_index']==0 and check['reference_index']=={'xy':1,'z':2,'joint':3}[direction]
                        and max(values)<=.1):
                    raise ValueError('Steady directional target/scope mismatch')
                steady.append([angle,direction]+values+[.1])
                evidence['steady'].append(dict(angle_deg=angle,direction=direction,metrics=check['errors'],
                                               source=f'{group}/summary.json#cases/new_directional_checks/{direction}'))
                steady_tex.append(f'${angle}^\\circ$ & '+{'xy':'XY','z':'Z','joint':r'\JointRefinement'}[direction]+
                                  ' & '+' & '.join(map(optical_display,values))+r'\\')
            if angle!=60:steady_tex.append(r'\midrule')
    if [r[0]for r in steady]!=[0]*3+[30]*3+[60]*3:
        raise ValueError('Steady angle coverage mismatch')
    output('fem_optical_steady',['angle_deg','direction','pv_percent','rms_over_refined_fem_pv_percent',
                                'maximum_over_refined_fem_pv_percent','target_percent'],steady,steady_tex)
    completion=read(COMPLETION/'run_summary.json')
    readback=completion['existing_transient_optical_refinement_readback']
    if readback['status']!='PASS' or readback['thermal_solves']!=0:
        raise ValueError('Transient archived-array readback did not pass')
    transient=[];transient_tex=[];timeseries=[]
    for angle in (0,30,60):
        group=f'transient_{angle}'
        summary=read(REPLAY/group/'summary.json');config=read(REPLAY/group/'config.json')
        verification_path=REPLAY/group/'verification.json';verification=read(verification_path)
        rb=next(r for r in readback['cases']if r['angle']==angle)
        if not (summary['status']=='PASS' and summary['readback']=='PASS' and summary['thermal_solves']==0
                and summary['angle_deg']==angle and verification['status']=='passed'
                and verification['eligible_for_standard_certificate_loader'] is False
                and sha(verification_path)==summary['verification_file_sha256']
                and rb['status']=='PASS' and rb['gated_rows']==36
                and rb['array_sha256']==summary['data_file_sha256']):
            raise ValueError('Transient replay provenance mismatch')
        manifest=config['manifest']
        if (manifest['space_tolerance'],manifest['joint_tolerance'],manifest['time_tolerance'],
            manifest['diagnostic_only_times_s'])!=(.005,.005,.001,[.001]):
            raise ValueError('Transient configured targets mismatch')
        if len(summary['optical_metrics'])!=40:
            raise ValueError('Transient directional coverage mismatch')
        for direction in ('xy','z','time','joint'):
            rows,maxima,diagnostic=transient_optical_rows(summary['optical_metrics'],direction)
            check=next(c for c in verification['checks']if c['name']==f'theta{angle}_{direction}_convergence')
            m=check['metrics']
            if (check['status']!='passed' or check['threshold']!=rows[0]['threshold']
                or m['times_s']!=list(OPTICAL_TIMES) or m['used_for_gate_by_time']!=[False]+[True]*9):
                raise ValueError('Transient verification scope mismatch')
            for row,verified in zip(rows,m['errors_by_time']):
                if any(row['metrics'][k]!=verified[k]for k in OPTICAL_KEYS):
                    raise ValueError('Summary/verification optical metric mismatch')
                timeseries.append([angle,direction,row['time_s'],row['gate']]+optical_percent(row['metrics'])+
                                  [100*row['threshold'],row['status']])
            evidence['transient'].append(dict(angle_deg=angle,direction=direction,records=rows,
                archived_optical_array_sha256=summary['data_file_sha256'],
                source_node_sha256=summary['source_array_sha256']))
            transient.append([angle,direction]+maxima+diagnostic+[100*rows[0]['threshold']])
            transient_tex.append(f'${angle}^\\circ$ & '+{'xy':'XY','z':'Z','time':r'\TimeRefinement','joint':r'\JointRefinement'}[direction]+
                                 ' & '+' & '.join(map(optical_display,maxima+diagnostic))+r'\\')
        if angle!=60:transient_tex.append(r'\midrule')
    output('fem_optical_transient',['angle_deg','direction','gated_pv_percent','gated_rms_percent','gated_maximum_percent',
        'diagnostic_1ms_pv_percent','diagnostic_1ms_rms_percent','diagnostic_1ms_maximum_percent','target_percent'],transient,transient_tex)
    with (OUT/'fem_optical_transient_timeseries.csv').open('w',encoding='utf-8',newline='')as f:
        writer=csv.writer(f);writer.writerow(['angle_deg','direction','time_s','used_for_target',
            'pv_percent','rms_over_refined_fem_pv_percent','maximum_over_refined_fem_pv_percent','target_percent','status'])
        writer.writerows(timeseries)
    loads=[]
    load_producer=COMPLETION/'complete_optical_evidence.py'
    sources[load_producer.relative_to(ROOT).as_posix()]=sha(load_producer)
    for angle in (30,60):
        result=read(COMPLETION/f'load_{angle}/summary.json')
        metrics=dict(zip(OPTICAL_KEYS,(result['optical_metrics'][k]for k in ('pv','rms','maximum'))))
        values=optical_percent(metrics)
        if not (result['status']=='PASS' and result['angle_deg']==angle and result['quadrature_orders']==[10,14]
                and result['threshold_fraction']==.0001 and max(values)<=.01
                and result==completion['added_steady_load_quadrature_checks'][str(angle)]):
            raise ValueError('Steady load quadrature check mismatch')
        loads.append([angle,10,14]+values+[.01])
        evidence['steady_load_quadrature'].append(result)
    with (OUT/'fem_optical_load_quadrature.csv').open('w',encoding='utf-8',newline='')as f:
        writer=csv.writer(f);writer.writerow(['angle_deg','base_order','refined_order','pv_percent',
            'rms_over_baseline_q10_fem_pv_percent','maximum_over_baseline_q10_fem_pv_percent','target_percent']);writer.writerows(loads)
    evidence.update(source_files_sha256=sources,controls=controls_result,
        scope='Retrieval and percent conversion of derived saved-node FEM optical re-evaluation; no new transient thermal solves. The Q10-to-Q14 load checks used additional completed steady solves.',
        certification='Derived evidence only; not a standard-loader-compatible full scientific certificate or a new physical validation.',
        denominator='Directional refinement: corresponding refined FEM OPD-PV; piston removed from field differences.',
        steady_load_denominator='Baseline order-10 FEM OPD-PV. Preserved from complete_optical_evidence.py: errors = metrics(optical, original), where original is baseline and optical is order-14. Distinct from directional-refinement denominators.',
        gate_times_s=list(OPTICAL_TIMES[1:]),diagnostic_only_times_s=[.001],
        targets_percent=dict(steady=.1,transient_space_and_joint=.5,transient_time=.1,steady_load_quadrature=.01),
        display='Values below 1e-5 percent shown as <1e-5 percent; CSV/JSON retain saved precision without an accuracy claim.',
        provenance_limit='Large archived node/optical arrays remain remote; their hashes and completed readback are retained. Local export checks summaries/configurations/derived verification only. Steady replay directories have no separate verification.json.',
        independence='Summary-versus-verification and completed readback use shared archived evidence; these are consistency/integrity checks, not independent validation of the underlying numerical method.')
    (OUT/'fem_optical_refinement_provenance.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'Exported {len(steady)} steady rows, {len(transient)} transient maxima/diagnostic rows, {len(timeseries)} time-resolved rows and {len(loads)} load checks; controls passed.')

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--optical-refinement',action='store_true',help='Export only the completed FEM optical-refinement tables; preserve historical tables')
    args=parser.parse_args()
    export_optical_refinement() if args.optical_refinement else main()

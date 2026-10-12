"""Run existing scientific experiments in dependency order, with honest evidence states.

This collects current experiments; COMPLETE means execution and saved-file checks
finished. It does not certify missing Gaussian FEM convergence or failed early-time
checks. Full runs return 2 when scientific evidence remains failed or incomplete.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import perf_counter

from thermal_opd.experiments import EXPERIMENTS, experiment_configuration
from thermal_opd.run_io import create_run, provenance, read_json, resolve_run_path, write_summary

ROOT=Path(__file__).resolve().parent


def build_plan(mode, angles):
    if mode not in ('smoke','full') or not angles or len(set(angles))!=len(angles) or not set(angles)<={0.,30.,60.}:
        raise ValueError('Choose smoke/full and distinct angles from the current selection scope: 0, 30, 60')
    stages=[]
    def add(name,experiment=None,args=(),category=None,needs=(),validate=True,target=None):
        stages.append(dict(name=name,experiment=experiment,args=list(args),category=category,
                           needs=list(needs),validate=validate,target=target))
    add('software_tests',validate=False)
    add('fem_reference','fem-validation',category='fem_verification')
    if mode=='full':
        add('steady_selection','convergence',category='convergence')
        add('transient_selection','transient-convergence',category='convergence')
        for angle in angles:
            add(f'semi_reference_theta{angle:g}','semi-validation',('--angle',str(angle),'--selection','@transient_selection'),
                category='semi_analytical_verification',needs=('transient_selection',))
        for angle in angles:
            add(f'steady_theta{angle:g}','steady',('--angle',str(angle)),category='steady_case',needs=('steady_selection',))
            add(f'transient_fields_theta{angle:g}','transient',('--mode','fields','--angle',str(angle),
                '--selection','@transient_selection'),category='transient_case',needs=('transient_selection',))
        add('steady_fem','fem-comparison',('--mode','full','--selection','@steady_selection'),
            category='fem_comparison',needs=('steady_selection',))
    for angle in (angles if mode=='full' else [max(angles)]):
        args=['--mode',mode,'--angle',str(angle),'--history','step']
        needs=[]
        if mode=='full':
            args+=['--selection','@transient_selection'];needs=['transient_selection']
        add(f'transient_fem_theta{angle:g}','transient',args,category='fem_comparison',needs=needs)
    if mode=='full':
        for name,category in (('azimuthal','azimuthal_inspection'),('source','heat_source_inspection'),
                              ('temperature','temperature_inspection'),('opd','opd_inspection'),
                              ('sweeps','convergence_diagnostic')):
            add(name,name,category=category,validate=False)
    fields=[s['name'] for s in stages if s['name'].startswith(('transient_fields_','transient_fem_'))]
    for name in fields:
        add('boundaries_'+name,'boundaries',('--results','@'+name),needs=(name,),target=name)
    if mode=='full':
        for name in [s['name'] for s in stages if s['name'].startswith(('steady_theta','transient_fields_'))]:
            add('replot_'+name,'replot',('--results','@'+name),category='replot',needs=(name,),validate=False)
    return stages


def _file_hash(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return dict(bytes=path.stat().st_size,sha256=digest.hexdigest())


def _code_hashes():
    files=[ROOT/'run_all.py',ROOT/'run.py',ROOT/'pyproject.toml']
    for folder in ('src','scripts','tests'):
        files.extend((ROOT/folder).rglob('*.py'))
    files.extend((ROOT/'tests/fixtures').glob('*.json'))
    return {p.relative_to(ROOT).as_posix():_file_hash(p) for p in sorted(set(files)) if p.is_file()}


def _artifacts(directory):
    # Boundary stages refresh sidecars/review timestamps on their own new runs.
    # Preserve hashes of all stable numerical data, configuration and figures.
    return {p.relative_to(directory).as_posix():_file_hash(p) for p in sorted(directory.rglob('*'))
            if p.is_file() and p.name not in ('verification.json','check_report.json')
            and not p.name.startswith('boundary_')}


def _run_command(command, log_file, *, append=False):
    with log_file.open('a' if append else 'x',encoding='utf-8') as stream:
        stream.write(json.dumps(command,ensure_ascii=False)+'\n');stream.flush()
        proc=subprocess.Popen(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,
            env={**os.environ,'PYTHONUNBUFFERED':'1','MPLBACKEND':'Agg'},start_new_session=os.name=='posix')
        try:
            return proc.wait()
        except KeyboardInterrupt:
            if os.name=='posix':os.killpg(proc.pid,signal.SIGTERM)
            else:proc.terminate()
            try:proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                if os.name=='posix':os.killpg(proc.pid,signal.SIGKILL)
                else:proc.kill()
                proc.wait()
            raise


def _check_steady_selection(directory):
    selected=read_json(directory/'summary.json')
    defaults=experiment_configuration('steady')['defaults']
    settings=dict(M=defaults['max_azimuthal_order'],N=defaults['num_radial_modes'],
        projection_r=defaults['num_projection_r'],projection_phi=defaults['num_projection_phi'],
        ray_points=defaults['num_opd_ray_points'])
    if selected['selected_settings']!=settings or selected['selected_figure_xy']>defaults['num_opd_xy']:
        raise ValueError('New steady selection differs from steady case settings; update the case settings before running it')


def _save(directory, report):
    write_summary(directory,report)
    with (directory/'data/stages.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=['name','status','scientific_status','run_directory','log_file','error'])
        writer.writeheader()
        for row in report['stages']:
            writer.writerow({key:row.get(key,'') for key in writer.fieldnames})


def run_suite(mode, angles, *, output_dir=None, resume=None, runner=None):
    runner=_run_command if runner is None else runner
    plan=build_plan(mode,angles)
    fingerprint=_code_hashes()
    if resume is not None:
        if output_dir is not None:raise ValueError('Use --resume or --output-dir, not both')
        directory=Path(resume).resolve()
        if directory.parent!=(ROOT/'outputs/experiment_suite').resolve():
            raise ValueError('Resume must name an outputs/experiment_suite/<ID> directory')
        report=read_json(directory/'summary.json')
        if report.get('mode')!=mode or report.get('angles_deg')!=list(angles):
            raise ValueError('Resume mode/angles differ from the recorded batch')
        if report.get('code_sha256')!=fingerprint:raise ValueError('Resume code changed; use a new batch for a new code version')
        if [r['name'] for r in report['stages']]!=[r['name'] for r in plan]:raise ValueError('Resume plan changed')
        for row in report['stages']:
            if row['status']=='COMPLETE' and row.get('run_directory'):
                path=resolve_run_path(directory,row['run_directory'])
                if not path.is_dir() or _artifacts(path)!=row['artifacts']:
                    raise ValueError(f"Completed data changed or is corrupt: {row['name']}")
    else:
        directory=create_run(ROOT,'experiment_suite',output_dir)
        (directory/'data/logs').mkdir()
        report=dict(kind='experiment_suite',status='RUNNING',mode=mode,angles_deg=list(angles),
            started_at_utc=datetime.now(timezone.utc).isoformat(),code_sha256=fingerprint,
            provenance=provenance(ROOT,ROOT/'run_all.py'),paper_ready=False,scientific_status='not_run',
            stages=[dict(name=s['name'],status='PENDING',attempts=[]) for s in plan],
            limits=['Existing experimental settings/thresholds are unchanged.',
                    'Current transient full mode does not establish final FEM space/time convergence.',
                    'Independent semi-analytical checks use this batch selection; unresolved checks retain failed/not_run states.',
                    'Intermediate transient history times are not individually precision-certified.',
                    'Development baselines are not needed or available on a fresh scientific checkout.'])
    print(f'Batch: {directory}',flush=True)
    report['status']='RUNNING';_save(directory,report)
    by_name={row['name']:row for row in report['stages']}
    for index,stage in enumerate(plan,1):
        row=by_name[stage['name']]
        if row['status']=='COMPLETE':
            print(f'[{index}/{len(plan)}] already complete: {stage["name"]}',flush=True);continue
        attempt=len(row['attempts'])+1
        log=directory/f'data/logs/{stage["name"]}_a{attempt}.log'
        result=None
        if stage['category']:
            result=ROOT/'outputs'/stage['category']/f'{directory.name}_{stage["name"]}_a{attempt}'
        elif stage['target']:
            result=resolve_run_path(directory,by_name[stage['target']]['run_directory'])
        row.update(status='RUNNING',scientific_status='not_run',log_file=str(log.relative_to(directory)),error=None)
        if result is not None:row['run_directory']=os.path.relpath(result,directory)
        row['attempts'].append(dict(number=attempt,started_at_utc=datetime.now(timezone.utc).isoformat(),log_file=row['log_file']))
        _save(directory,report)
        started=perf_counter()
        print(f'[{index}/{len(plan)}] running: {stage["name"]}; log={log}',flush=True)
        try:
            for name in stage['needs']:
                if by_name[name]['status']!='COMPLETE':raise ValueError(f'Prerequisite not complete: {name}')
            if stage['experiment']=='steady':
                _check_steady_selection(resolve_run_path(directory,by_name['steady_selection']['run_directory']))
            if stage['experiment'] is None:
                command=[sys.executable,'-m','pytest','-q','-p','no:cacheprovider',
                         '--basetemp',str(directory/f'data/pytest_a{attempt}')]
            else:
                args=[str(resolve_run_path(directory,by_name[a[1:]]['run_directory'])) if a.startswith('@') else a
                      for a in stage['args']]
                command=[sys.executable,'-u',str(ROOT/'scripts'/EXPERIMENTS[stage['experiment']][1]),*args]
                if stage['category']:command+=['--output-dir',str(result)]
            row['command']=command
            code=runner(command,log)
            row['returncode']=code
            if code:raise RuntimeError(f'Command exited {code}; inspect {log}')
            if result is not None:
                summary=read_json(result/'summary.json')
                if summary.get('status')!='COMPLETE':raise ValueError('Experiment did not record COMPLETE')
                if stage['validate']:
                    command=[sys.executable,'-u',str(ROOT/'scripts'/EXPERIMENTS['validate'][1]),'--results',str(result)]
                    code=runner(command,log,append=True)
                    if code:raise RuntimeError(f'Saved-result validation exited {code}; inspect {log}')
                    check=read_json(result/'check_report.json')
                    if check.get('status')!='PASS':raise ValueError('Saved-result check did not pass')
                    row['scientific_status']=check['verification']['status']
                    row['integrity_status']='PASS'
                else:row['integrity_status']='metadata_only'
                row['artifacts']=_artifacts(result)
                row['mesh_convergence_established']=summary.get('mesh_convergence_established')
                row['fem_convergence_established']=summary.get('fem_convergence_established')
            row['status']='COMPLETE'
        except KeyboardInterrupt:
            row['status']='ABORTED';report['status']='ABORTED';_save(directory,report)
            return 130,directory
        except Exception as error:
            row.update(status='FAILED',error=str(error));report['status']='FAILED'
            _save(directory,report)
            print(f'FAILED: {stage["name"]}: {error}',flush=True)
            return 1,directory
        finally:
            row['attempts'][-1].update(elapsed_seconds=perf_counter()-started,status=row['status'])
            _save(directory,report)
    science=[r['scientific_status'] for s,r in zip(plan,report['stages']) if s['validate'] and s['experiment']]
    convergence=(mode=='full' and by_name['steady_fem'].get('mesh_convergence_established') is True
                 and all(by_name[f'transient_fem_theta{a:g}'].get('fem_convergence_established') is True for a in angles))
    report['scientific_status']='failed' if 'failed' in science else ('passed' if science and all(v=='passed' for v in science) and convergence else 'not_run')
    report.update(status='COMPLETE',finished_at_utc=datetime.now(timezone.utc).isoformat(),
                  paper_ready=mode=='full' and report['scientific_status']=='passed')
    _save(directory,report)
    print(f'COMPLETE execution: {directory}; scientific={report["scientific_status"]}; paper_ready={report["paper_ready"]}',flush=True)
    return (2 if mode=='full' and not report['paper_ready'] else 0),directory


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('smoke','full'),default='smoke')
    parser.add_argument('--angles',type=float,nargs='+',default=[0.,30.,60.])
    parser.add_argument('--plan',action='store_true',help='Print the command/dependency plan without computing or creating outputs')
    parser.add_argument('--output-dir',type=Path,help='New outputs/experiment_suite/<ID>; child runs stay in their experiment categories')
    parser.add_argument('--resume',type=Path,help='Resume this exact code and mode; completed data are hash-checked, failed attempts get new output IDs')
    args=parser.parse_args(argv)
    try:
        plan=build_plan(args.mode,args.angles)
        if args.plan:
            print(json.dumps(dict(mode=args.mode,angles_deg=args.angles,stages=plan),ensure_ascii=False,indent=2));return 0
        code,_=run_suite(args.mode,args.angles,output_dir=args.output_dir,resume=args.resume)
        return code
    except (ValueError,OSError,KeyError) as error:
        print(f'ERROR: {error}',file=sys.stderr);return 1


if __name__=='__main__':
    def interrupted(signum,frame):raise KeyboardInterrupt
    signal.signal(signal.SIGTERM,interrupted)
    raise SystemExit(main())

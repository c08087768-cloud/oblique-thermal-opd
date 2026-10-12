"""Batch orchestration must preserve failures and use this batch's selections."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture
def suite():
    spec=importlib.util.spec_from_file_location('_suite_test',ROOT/'run_all.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_full_plan_has_valid_dependency_order_and_all_scientific_entries(suite):
    stages=suite.build_plan('full',[0.,30.,60.])
    completed=set()
    for stage in stages:
        assert set(stage['needs'])<=completed
        completed.add(stage['name'])
    assert len(completed)==len(stages)
    assert {s['experiment'] for s in stages if s['experiment']}=={
        'convergence','transient-convergence','steady','transient','fem-comparison',
        'semi-validation','fem-validation','boundaries','source','azimuthal',
        'temperature','opd','sweeps','replot'}
    assert sum(s['name'].startswith('transient_fields_') for s in stages)==3
    assert sum(s['name'].startswith('transient_fem_') for s in stages)==3
    assert not any('baseline' in str(s) for s in stages)


def test_smoke_plan_does_not_launch_selection_full_fem_or_dense_fields(suite):
    stages=suite.build_plan('smoke',[60.])
    assert {s['experiment'] for s in stages if s['experiment']}=={'fem-validation','transient','boundaries'}
    assert not any('full' in s['args'] or 'fields' in s['args'] for s in stages)


def test_plan_only_does_not_create_outputs_or_launch_commands(suite,tmp_path,monkeypatch,capsys):
    monkeypatch.setattr(suite,'ROOT',tmp_path)
    monkeypatch.setattr(suite,'_run_command',lambda *a,**kw:pytest.fail('plan started a command'))
    assert suite.main(['--mode','full','--plan'])==0
    assert not (tmp_path/'outputs').exists()
    assert json.loads(capsys.readouterr().out)['mode']=='full'


def fake_project(suite,tmp_path,monkeypatch):
    monkeypatch.setattr(suite,'ROOT',tmp_path)
    (tmp_path/'run_all.py').write_text('batch source snapshot')
    monkeypatch.setattr(suite,'provenance',lambda *a:dict(git_revision='test',working_tree_dirty=False))
    monkeypatch.setattr(suite,'_check_steady_selection',lambda *a:None)


def fake_executor(suite,calls,*,science='not_run',fail_name=None):
    def run(command,log_file,*,append=False):
        calls.append(command)
        if '-m' in command:return 0
        script=Path(command[2]).name
        if fail_name and script==fail_name:return 1
        if script=='validate_run.py':
            path=Path(command[command.index('--results')+1])
            (path/'check_report.json').write_text(json.dumps(dict(status='PASS',verification={'status':science})))
            return 0
        if '--output-dir' in command:
            path=Path(command[command.index('--output-dir')+1]);path.mkdir(parents=True)
            (path/'data').mkdir();(path/'figures').mkdir()
            summary=dict(status='COMPLETE',kind='mock',selected_settings=dict(M=6,N=8,projection_r=400,projection_phi=32,ray_points=16))
            (path/'summary.json').write_text(json.dumps(summary))
            (path/'verification.json').write_text(json.dumps(dict(status=science)))
            (path/'data/fields.npz').write_bytes(b'stored fields')
        return 0
    return run


def test_smoke_reports_execution_success_without_scientific_certification(suite,tmp_path,monkeypatch):
    fake_project(suite,tmp_path,monkeypatch)
    calls=[]
    code,path=suite.run_suite('smoke',[60.],runner=fake_executor(suite,calls))
    report=json.loads((path/'summary.json').read_text())
    assert code==0 and report['status']=='COMPLETE'
    assert report['scientific_status']=='not_run' and report['paper_ready'] is False
    assert all('legacy' not in str(command) for command in calls)


def test_failed_scientific_evidence_does_not_stop_full_data_collection_or_become_passed(suite,tmp_path,monkeypatch):
    fake_project(suite,tmp_path,monkeypatch)
    calls=[]
    code,path=suite.run_suite('full',[60.],runner=fake_executor(suite,calls,science='failed'))
    report=json.loads((path/'summary.json').read_text())
    assert code==2 and report['status']=='COMPLETE'
    assert report['scientific_status']=='failed' and report['paper_ready'] is False
    assert all(s['status']=='COMPLETE' for s in report['stages'])
    selections={s['name']:path/s['run_directory'] for s in report['stages'] if s['name'].endswith('_selection')}
    for command in calls:
        if '--selection' in command:
            selected=Path(command[command.index('--selection')+1])
            assert selected.resolve() in {p.resolve() for p in selections.values()}


def test_resume_skips_completed_stages_and_keeps_artifacts(suite,tmp_path,monkeypatch):
    fake_project(suite,tmp_path,monkeypatch)
    calls=[]
    _,path=suite.run_suite('smoke',[60.],runner=fake_executor(suite,calls))
    calls.clear()
    code,same=suite.run_suite('smoke',[60.],resume=path,runner=fake_executor(suite,calls))
    assert code==0 and same==path and not calls


def test_resume_rejects_changed_code_or_corrupt_completed_data(suite,tmp_path,monkeypatch):
    fake_project(suite,tmp_path,monkeypatch)
    _,path=suite.run_suite('smoke',[60.],runner=fake_executor(suite,[]))
    report=json.loads((path/'summary.json').read_text())
    artifact=next(s for s in report['stages'] if s.get('run_directory'))
    (path/artifact['run_directory']/'data/fields.npz').write_bytes(b'corruption')
    with pytest.raises(ValueError,match='changed|corrupt'):
        suite.run_suite('smoke',[60.],resume=path,runner=fake_executor(suite,[]))
    (tmp_path/'run_all.py').write_text('modified code')
    with pytest.raises(ValueError,match='code'):
        suite.run_suite('smoke',[60.],resume=path,runner=fake_executor(suite,[]))


def test_failed_process_stops_dependents_and_resume_uses_new_output(suite,tmp_path,monkeypatch):
    fake_project(suite,tmp_path,monkeypatch)
    code,path=suite.run_suite('smoke',[60.],runner=fake_executor(suite,[],fail_name='run_transient_case.py'))
    before=json.loads((path/'summary.json').read_text())
    assert code==1 and before['status']=='FAILED'
    failed=next(s for s in before['stages'] if s['status']=='FAILED')
    old=failed['run_directory']
    calls=[]
    code,_=suite.run_suite('smoke',[60.],resume=path,runner=fake_executor(suite,calls))
    after=json.loads((path/'summary.json').read_text())
    retried=next(s for s in after['stages'] if s['name']==failed['name'])
    assert code==0 and retried['run_directory']!=old
    assert len(retried['attempts'])==2
    assert all(Path(c[2]).name!='verify_fem.py' for c in calls)


def test_steady_case_must_match_new_selected_parameters(suite,tmp_path,monkeypatch):
    from thermal_opd.experiments import experiment_configuration
    config=experiment_configuration('steady')['defaults']
    settings=dict(M=config['max_azimuthal_order'],N=config['num_radial_modes'],projection_r=config['num_projection_r'],
                  projection_phi=config['num_projection_phi'],ray_points=config['num_opd_ray_points'])
    selected=tmp_path/'summary.json'
    selected.write_text(json.dumps(dict(selected_settings=settings,selected_figure_xy=121)))
    suite._check_steady_selection(tmp_path)
    settings['N']+=1
    selected.write_text(json.dumps(dict(selected_settings=settings,selected_figure_xy=121)))
    with pytest.raises(ValueError,match='steady|settings'):
        suite._check_steady_selection(tmp_path)

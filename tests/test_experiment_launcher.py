import pytest
from thermal_opd.experiments import main,experiment_configuration


def test_interactive_cancel_never_launches_an_experiment(monkeypatch):
    import thermal_opd.experiments as experiments
    def forbidden(*args,**kwargs):
        pytest.fail('Cancelled launcher started computation')
    monkeypatch.setattr(experiments.subprocess,'run',forbidden)
    replies=iter(['1','n'])
    output=[]
    assert main([],input_fn=lambda prompt:next(replies),print_fn=output.append)==0
    assert any('取消' in str(line) for line in output)


def test_direct_launcher_uses_same_script_and_keeps_arguments(monkeypatch):
    import thermal_opd.experiments as experiments
    calls=[]
    class Result:
        returncode=0
    monkeypatch.setattr(experiments.subprocess,'run',lambda command,**kw:(calls.append(command) or Result()))
    assert main(['--experiment','transient','--yes','--mode','fields','--angle','30','--selection','somewhere'],print_fn=lambda *a:None)==0
    assert calls[0][1].endswith('run_transient_case.py')
    assert calls[0][2:]==['--mode','fields','--angle','30','--selection','somewhere']


def test_configuration_preview_does_not_create_output():
    config=experiment_configuration('steady')
    assert config['physical_parameters_SI']['radius']==.0127
    assert config['defaults']['max_azimuthal_order']==6

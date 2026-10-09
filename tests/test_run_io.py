import json
import numpy as np
import pytest
from thermal_opd.run_io import create_run,write_summary,data_directory,figure_directory


def test_default_runs_are_classified_unique_and_never_overwrite(tmp_path):
    a=create_run(tmp_path,'steady_case')
    b=create_run(tmp_path,'steady_case')
    assert a!=b
    assert a.parent==b.parent==tmp_path/'outputs/steady_case'
    assert data_directory(a)==a/'data'
    assert figure_directory(a)==a/'figures'
    assert data_directory(a).is_dir() and figure_directory(a).is_dir()
    angled=create_run(tmp_path,'steady_refracted_theta30.5')
    assert angled.parent==tmp_path/'outputs/steady_case'
    with pytest.raises(FileExistsError):
        create_run(tmp_path,'steady_case',a)
    with pytest.raises(ValueError):
        create_run(tmp_path,'../escape')


def test_metadata_keeps_configuration_metrics_and_unexecuted_status_separate(tmp_path):
    path=create_run(tmp_path,'transient_case')
    summary=dict(kind='example',status='RUNNING',physical_parameters_SI={'radius':.01},numerical_settings={'M':6},provenance={'git_revision':'abc'})
    write_summary(path,summary)
    config=json.loads((path/'config.json').read_text())
    assert config['physical_parameters_SI']=={'radius':.01}
    assert config['numerical_settings']=={'M':6}
    summary.update(status='COMPLETE',opd_pv_m=1e-7)
    write_summary(path,summary)
    assert json.loads((path/'metrics.json').read_text())['opd_pv_m']==1e-7
    assert json.loads((path/'verification.json').read_text())['status']=='not_run'
    assert json.loads((path/'config.json').read_text())['provenance']['git_revision']=='abc'


def test_legacy_flat_artifact_layout_remains_readable(tmp_path):
    assert data_directory(tmp_path)==tmp_path
    assert figure_directory(tmp_path)==tmp_path

from dataclasses import replace
import numpy as np
import pytest
from thermal_opd.parameters import PhysicalParameters
from thermal_opd.verification import verification_check,verification_report,finite_difference_pde_residual


def test_reports_do_not_promote_missing_scientific_evidence_to_pass():
    passed=verification_check('boundary','passed',metrics={'max':1e-12},threshold=1e-8)
    missing=verification_check('mesh','not_run',scope='Requires independent mesh refinement')
    assert verification_report([passed,missing],config={},run_directory='run')['status']=='not_run'
    assert verification_report([passed],config={},run_directory='run')['status']=='passed'
    failed=verification_check('pde','failed',metrics={'max':1.})
    assert verification_report([passed,missing,failed],config={},run_directory='run')['status']=='failed'
    with pytest.raises(ValueError):
        verification_check('mesh','PASS')


def test_independent_pde_check_detects_wrong_time_response():
    p=replace(PhysicalParameters(),radius=1.,thickness=1.,conductivity=1.,density=2.,specific_heat=3.)
    r,phi,z=np.array([.2,.5]),np.array([.1,1.2]),np.array([.3,.6])
    correct=lambda r,phi,z,t:r*r+z*z+t
    wrong=lambda r,phi,z,t:r*r+z*z+2*t
    args=(p,r,phi,z,1.)
    residual=finite_difference_pde_residual(correct,*args,space_step=1e-3,time_step=1e-3)
    np.testing.assert_allclose(residual,0.,atol=2e-8)
    np.testing.assert_allclose(finite_difference_pde_residual(wrong,*args,space_step=1e-3,time_step=1e-3),6.,atol=2e-8)


def test_independent_one_dimensional_reference_has_correct_steady_and_initial_limits():
    from thermal_opd.verification import uniform_flux_1d_reference
    z,T=uniform_flux_1d_reference(1.,2.,3.,4.,5.,6.,[0.,1e5],num_intervals=40)
    np.testing.assert_array_equal(T[0],np.zeros_like(z))
    denominator=2.*(5.+6.)+5.*6.*1.
    expected=4.*(2.+6.*(1.-z))/denominator
    np.testing.assert_allclose(T[-1],expected,rtol=2e-12)


def test_scientific_artifact_validator_recomputes_evidence_and_preserves_failed_status(tmp_path):
    import json
    from thermal_opd.run_validation import validate_saved_run
    field=np.array([.01,.02])
    report=verification_report([verification_check('pde','failed',metrics={'max_error':.02},threshold=.001,evidence='pde')],
                               config={},run_directory=tmp_path)
    summary=dict(kind='scientific_verification',status='COMPLETE',verification_status='failed')
    np.savez(tmp_path/'verification_fields.npz',pde=field)
    (tmp_path/'verification.json').write_text(json.dumps(report))
    (tmp_path/'summary.json').write_text(json.dumps(summary))
    result=validate_saved_run(tmp_path)
    assert result['status']=='PASS'  # File integrity passes, scientific test has failed.
    assert result['verification']['status']=='failed'
    report['checks'][0]['status']='passed'
    (tmp_path/'verification.json').write_text(json.dumps(report))
    with pytest.raises(ValueError,match='Scientific|scientific'):
        validate_saved_run(tmp_path)

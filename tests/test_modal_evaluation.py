import numpy as np


def test_shared_modal_evaluation_broadcasts_and_preserves_derivatives():
    from thermal_opd import modal_evaluation
    c,roots = {0:np.array([1.]),2:np.array([2.])},{0:np.array([80.]),2:np.array([300.])}
    r,phi,z = np.linspace(0,.01,5)[:,None,None],np.arange(8)[None,:,None]*np.pi/4,np.array([0,.002,.005])[None,None,:]
    value = lambda m,n,z:np.exp(-roots[m][n]*z)
    response = lambda m,n,z:(value(m,n,z),-roots[m][n]*value(m,n,z))
    actual = modal_evaluation.modal_temperature(r,phi,z,c,roots,value)
    dr,dz = modal_evaluation.modal_temperature_gradient(r,phi,z,c,roots,response)
    assert actual.shape==dr.shape==dz.shape==(5,8,3)
    eps = 1e-7
    expected = (modal_evaluation.modal_temperature(r,phi,z+eps,c,roots,value)-
                modal_evaluation.modal_temperature(r,phi,z-eps,c,roots,value))/(2*eps)
    np.testing.assert_allclose(dz,expected,rtol=2e-8,atol=1e-8)

import pytest
from lateral_tools.codes.tier1_checks import capacity_shear_stress_rcmf, check_shear_stress_rcmf

def test_capacity_shear_stress_rcmf():
    assert capacity_shear_stress_rcmf('IO', 4000) == pytest.approx(max(2*(4000**0.5),100))
    assert capacity_shear_stress_rcmf('LS', 4000) == 3*(4000**0.5)
    assert capacity_shear_stress_rcmf('CP', 4000) == 4*(4000**0.5)

def test_check_shear_stress_rcmf():
    # Mock inputs
    res = check_shear_stress_rcmf(V_j=100, A_c=200, n_c=10, n_f=2, perf='IO', f_c=4000)
    assert 'demand_psi' in res
    assert 'capacity_psi' in res
    assert res['dc_ratio'] > 0
    assert isinstance(res['pass'], bool)


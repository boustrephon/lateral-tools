"""Tests for IS 1893 (Part 1) — 2016 edition and 2023 draft."""

import pytest

from lateral_tools.codes.IS1893_2016 import (
    T_aa,
    get_A_H_2016,
    get_A_H_2016_with_minima,
    get_A_H_2023_with_minima,
    get_A_NH5,
    get_d_ksi,
    get_rho,
    get_RP_pairs,
    get_S_a_g,
    get_site_class,
    get_T_a,
    get_V_B_2016,
    get_V_BDH,
    get_Z_2016,
    get_Z_2023,
)


# ----------------------------------------------------------------------
# IS 1893-2016
# ----------------------------------------------------------------------
def test_get_Z_2016():
    assert get_Z_2016('II') == 0.10
    assert get_Z_2016('III') == 0.16
    assert get_Z_2016('IV') == 0.24
    assert get_Z_2016('V') == 0.36
    assert get_Z_2016('nope') is None


def test_get_rho():
    assert get_rho('II') == 0.007
    assert get_rho('III') == 0.011
    assert get_rho('IV') == 0.016
    assert get_rho('V') == 0.024


def test_get_S_a_g_ESM():
    # Plateau
    assert get_S_a_g(0.2, 'A', 'ESM') == pytest.approx(2.5)
    # 1/T decay
    assert get_S_a_g(2.0, 'A', 'ESM') == pytest.approx(0.5)
    # Constant tail
    assert get_S_a_g(5.0, 'A', 'ESM') == pytest.approx(0.25)


def test_get_S_a_g_RSM():
    # Rising branch
    assert get_S_a_g(0.05, 'B', 'RSM') == pytest.approx(1 + 15 * 0.05)


def test_get_A_H_2016():
    assert get_A_H_2016(0.16, 2.5, 1.0, 4.0) == pytest.approx(
        0.5 * 0.16 * 2.5 / 4.0
    )


def test_get_A_H_2016_with_minima():
    # Without minima, plain Ah
    assert get_A_H_2016_with_minima(0.2, 'II', 'A', I_e=1.0, R=1.0,
                                    method='ESM', with_minima=False) == pytest.approx(
        0.5 * 0.10 * 2.5
    )
    # With minima, floored at rho
    assert get_A_H_2016_with_minima(6.0, 'II', 'A', I_e=1.0, R=8.0,
                                    method='ESM', with_minima=True) == pytest.approx(0.007)


def test_get_V_B_2016():
    W = 1000.0
    V = get_V_B_2016(W, 0.2, 'III', 'A', I_e=1.0, R=5.0, method='ESM')
    A_H = 0.5 * 0.16 * 2.5 / 5.0
    assert V == pytest.approx(max(A_H, 0.011) * W)


# ----------------------------------------------------------------------
# IS 1893-2023 draft
# ----------------------------------------------------------------------
def test_get_RP_pairs_int():
    name, pair = get_RP_pairs(1)
    assert name == 'Normal Structures'
    assert pair == (475, 73)


def test_get_RP_pairs_str():
    name, pair = get_RP_pairs('Important Structures')
    assert name == 'Important Structures'
    assert pair == (975, 225)


def test_get_Z_2023():
    # ULS (475-yr) zone factor for zone II
    assert get_Z_2023('II', 1, 'ULS') == pytest.approx(0.075)
    # SLS (73-yr)
    assert get_Z_2023('II', 1, 'SLS') == pytest.approx(0.0375)


def test_get_site_class():
    assert get_site_class(2000) == 'A'
    assert get_site_class(1000) == 'B'
    assert get_site_class(500) == 'C'
    assert get_site_class(190) == 'D'


def test_get_A_NH5_LESM():
    assert get_A_NH5(0.2, 'A', 'LESM') == pytest.approx(2.5)


def test_get_A_NH5_LDM():
    assert get_A_NH5(0.005, 'A', 'LDM') == pytest.approx(1.0)


def test_get_A_H_2023_with_minima():
    # Plain design value (no minima)
    ahd = get_A_H_2023_with_minima(0.2, 'II', 'A', I_e=1.0, R=1.0,
                                   category=1, LS='ULS', method='LESM',
                                   with_minima=False)
    Z = 0.075  # zone II, ULS
    assert ahd == pytest.approx(Z * 2.5)
    # Minima floor
    ahd = get_A_H_2023_with_minima(0.2, 'II', 'A', I_e=1.0, R=100.0,
                                   category=1, LS='ULS', method='LESM',
                                   with_minima=True)
    assert ahd == pytest.approx(0.015)


def test_get_V_BDH():
    W = 1000.0
    V = get_V_BDH(W, 0.2, 'II', 'A', I_e=1.0, R=1.0, category=1, LS='ULS',
                  method='LESM')
    Z = 0.075
    assert V == pytest.approx(max(Z * 2.5, 0.625 * Z, 0.015) * W)


def test_get_d_ksi():
    # eta_max for 5% damping
    eta_max = (7 / (2 + 0.05)) ** 0.6
    # Below 0.01 s -> unity
    assert get_d_ksi(0.05, 0.005, 'B') == pytest.approx(1.0)
    # On the plateau (0.1 <= T <= T_AV) -> eta_max
    assert get_d_ksi(0.05, 0.4, 'B') == pytest.approx(eta_max)
    # At 6.0 s -> decays back to unity
    assert get_d_ksi(0.05, 6.0, 'B') == pytest.approx(1.0)
    # Beyond 6.0 s -> unity
    assert get_d_ksi(0.05, 10.0, 'B') == pytest.approx(1.0)


def test_get_d_ksi_errors():
    with pytest.raises(ValueError):
        get_d_ksi(0.05, 1.0, 'Z')
    with pytest.raises(ValueError):
        get_d_ksi(0.40, 1.0, 'B')


# ----------------------------------------------------------------------
# Periods
# ----------------------------------------------------------------------
def test_get_T_a():
    assert get_T_a(10.0, 0) == pytest.approx(0.075 * 10.0 ** 0.75)


def test_T_aa():
    assert T_aa(10.0, 10.0) == pytest.approx(0.09 * 10.0 * 10.0 ** -0.5)

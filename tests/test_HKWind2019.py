"""Tests for HK Wind Code 2019 functions (lateral_tools.codes.HKWind2019)."""

import math

import pytest

from lateral_tools.codes.HKWind2019 import (
    freestanding_wall_pressure,
    get_C_f,
    get_Q_oz,
    get_Q_z,
    get_S_q,
    get_Ss,
    get_table_4_2,
)


def test_get_Q_oz_eq_3_2():
    assert get_Q_oz(3.0) == pytest.approx(3.7 * (3.0 / 500.0) ** 0.16)
    assert get_Q_oz(500.0) == pytest.approx(3.7)
    with pytest.raises(ValueError):
        get_Q_oz(501.0)


def test_get_Q_z_eq_3_1():
    q_oz = get_Q_oz(3.0)
    assert get_Q_z(q_oz, 1.0, 0.85) == pytest.approx(q_oz * 1.0 * 0.85)


def test_get_C_f_eq_4_1():
    he_d = 3.0 / 6.058
    b_d = 2.438 / 6.058
    expected = 1.1 + 0.055 * he_d / math.exp(
        abs(math.log(0.6 * b_d * (1 - 0.011 * he_d))) ** (1.7 - 0.0013 * he_d ** 2)
    )
    assert get_C_f(he_d, b_d) == pytest.approx(expected)
    with pytest.raises(ValueError):
        get_C_f(13.0, 1.0)


def test_get_Ss_overall_eq_c1_1a():
    l05p = 3.0 / 1.5 + 2 * 2.438
    assert get_Ss(l05p) == pytest.approx(math.exp(0.17 - 0.07 * l05p ** 0.32))


def test_get_Ss_edge_corner():
    l05p = 3.0 + 2.438  # opening-frame / element L0.5p = H + B
    s_overall = math.exp(0.17 - 0.07 * l05p ** 0.32)
    assert get_Ss(l05p, "Edge") == pytest.approx(max(s_overall, 1.3 - math.log(l05p) / 9.0))
    assert get_Ss(l05p, "Corner") == pytest.approx(max(s_overall, 1.5 - math.log(l05p) / 5.4))


def test_get_S_q_eq_5_3():
    h, b = 3.0, 2.438
    l05p = h / 1.5 + 2 * b
    assert get_S_q(h, b) == pytest.approx(1.1 * get_Ss(l05p))
    with pytest.raises(ValueError):
        get_S_q(51.0, b)


def test_get_table_4_2_table_values():
    table = {0.01: 2.0, 0.1: 1.8, 0.2: 1.7, 0.3: 1.6,
             0.5: 1.5, 0.8: 1.5, 0.9: 1.6, 1.0: 2.0}
    for phi, cp in table.items():
        assert get_table_4_2(phi) == pytest.approx(cp)


def test_get_table_4_2_interpolation():
    # 0.4 lies between 0.3 (1.6) and 0.5 (1.5)
    assert get_table_4_2(0.4) == pytest.approx(1.55)
    # 0.65 lies between 0.5 (1.5) and 0.8 (1.5)
    assert get_table_4_2(0.65) == pytest.approx(1.5)
    # bookends outside the table
    assert get_table_4_2(0.005) == pytest.approx(2.0)
    assert get_table_4_2(1.5) == pytest.approx(2.0)


def test_freestanding_wall_pressure_no_return():
    # phi=1, l/h=7, no return corners -> interpolate l/h between 5 and 10
    cp = freestanding_wall_pressure(1.0, 7.0, 0.0)
    assert cp['A'] == pytest.approx(3.1)
    assert cp['B'] == pytest.approx(1.92)
    assert cp['C'] == pytest.approx(1.52)
    assert cp['D'] == pytest.approx(1.2)


def test_freestanding_wall_pressure_full_return():
    cp = freestanding_wall_pressure(1.0, 7.0, 1.0)
    assert cp == {'A': 2.1, 'B': 1.8, 'C': 1.4, 'D': 1.2}


def test_freestanding_wall_pressure_phi_08():
    cp = freestanding_wall_pressure(0.8, 7.0, 0.0)
    assert cp == {'A': 1.2, 'B': 1.2, 'C': 1.2, 'D': 1.2}


def test_freestanding_wall_pressure_raises():
    with pytest.raises(ValueError):
        freestanding_wall_pressure(0.5, 7.0, 0.0)

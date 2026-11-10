"""Tests for unit conversion utilities.

Covers `scale_to_SI_units` and the `units_SI_scaling_dict` lookup table
in `lateral_tools.math_utils`.
"""

import pytest

from lateral_tools.math_utils import scale_to_SI_units, units_SI_scaling_dict


# -------------------------------------------------------------------
# Length (to metres)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("mi", 1609.344),
        ("yd", 0.9144),
        ("ft", 0.3048),
        ("in", 0.0254),
        ("km", 1000.0),
        ("m", 1.0),
        ("cm", 0.01),
        ("mm", 0.001),
    ],
)
def test_length_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


# -------------------------------------------------------------------
# Force (to newtons)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("GN", 1.0e9),
        ("MN", 1.0e6),
        ("kN", 1000.0),
        ("N", 1.0),
    ],
)
def test_force_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


# -------------------------------------------------------------------
# Moment (to newton-metres)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("GNm", 1.0e9),
        ("MNm", 1.0e6),
        ("kNm", 1000.0),
        ("Nm", 1.0),
        ("Nmm", 0.001),
    ],
)
def test_moment_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


# -------------------------------------------------------------------
# Distributed load (to newtons per metre)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("kN/m", 1000.0),
        ("N/m", 1.0),
        ("N/mm", 1000.0),
    ],
)
def test_distributed_load_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


# -------------------------------------------------------------------
# Pressure / stress (to pascals)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("GPa", 1.0e9),
        ("MPa", 1.0e6),
        ("kPa", 1000.0),
        ("Pa", 1.0),
        ("kN/m2", 1000.0),
        ("N/m2", 1.0),
        ("N/mm2", 1.0e6),
        ("kN/m^2", 1000.0),
        ("N/m^2", 1.0),
        ("N/mm^2", 1.0e6),
    ],
)
def test_pressure_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


# -------------------------------------------------------------------
# US customary pressure units
# -------------------------------------------------------------------
def test_ksf_psf_values():
    """1 ksf (kip/sq ft) and 1 psf (lb/sq ft) scaled to pascals."""
    assert scale_to_SI_units("ksf") == pytest.approx(47880.26, rel=1e-6)
    assert scale_to_SI_units("psf") == pytest.approx(47.88026, rel=1e-6)


def test_ksf_psf_keys_present():
    """Regression test: 'psf' previously had no value, silently merging
    with 'GPa' to form a bogus 'psfGPa' key."""
    assert "psf" in units_SI_scaling_dict
    assert "ksf" in units_SI_scaling_dict
    assert "psfGPa" not in units_SI_scaling_dict


# -------------------------------------------------------------------
# Velocity (to metres per second)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("km/s", 1000.0),
        ("m/s", 1.0),
        ("cm/s", 0.01),
        ("mm/s", 0.001),
        ("ft/s", 0.3048),
        ("in/s", 0.0254),
    ],
)
def test_velocity_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


# -------------------------------------------------------------------
# Acceleration (to metres per second squared)
# -------------------------------------------------------------------
@pytest.mark.parametrize(
    "unit, expected",
    [
        ("km/s2", 1000.0),
        ("m/s2", 1.0),
        ("cm/s2", 0.01),
        ("mm/s2", 0.001),
        ("ft/s2", 0.3048),
        ("in/s2", 0.0254),
        ("km/s^2", 1000.0),
        ("m/s^2", 1.0),
        ("cm/s^2", 0.01),
        ("mm/s^2", 0.001),
        ("ft/s^2", 0.3048),
        ("in/s^2", 0.0254),
    ],
)
def test_acceleration_units(unit, expected):
    assert scale_to_SI_units(unit) == pytest.approx(expected)


def test_gravity_units():
    assert scale_to_SI_units("g") == pytest.approx(9.80665)


# -------------------------------------------------------------------
# Lookup behaviour
# -------------------------------------------------------------------
def test_case_insensitive_lookup():
    assert scale_to_SI_units("KN") == pytest.approx(1000.0)
    assert scale_to_SI_units("Mpa") == pytest.approx(1.0e6)
    assert scale_to_SI_units("PSF") == pytest.approx(47.88026, rel=1e-6)


def test_unknown_unit_returns_none():
    assert scale_to_SI_units("not_a_unit") is None


def test_unknown_unit_default():
    assert scale_to_SI_units("not_a_unit", default=1.0) == 1.0
    assert scale_to_SI_units("not_a_unit", default=0) == 0

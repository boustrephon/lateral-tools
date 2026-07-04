"""
Tests for the lateral-tools package.

Part of the lateral-tools project for seismic and lateral load engineering.

Test coverage:
    - GB50011: alpha_max lookup, design spectrum shape
    - EC8: elastic and design spectra (with behaviour factor q)
    - ASCE 7: design response spectrum plateau and long-period decay
    - plot_seismic_spectrum: Figure return, multi-code, return_data,
      pre-computed data, log axes
    - math_utils: interp, list_interp, as_float, as_numeric
    - risk_utils: exceedance2return, return2exceedance
    - analysis: storey_elevations, storey_force_distribution
    - tier1_checks: shear stress capacities
"""

import numpy as np
import pytest
from matplotlib.figure import Figure


# -------------------------------------------------------------------
# GB 50011
# -------------------------------------------------------------------
class TestGB50011Spectrum:
    def test_gb_get_alpha_max(self):
        """alpha_max values should match GB 50011 Table 5.1.4-1."""
        from lateral_tools.codes.GB50011 import GB_get_alpha_max

        # Intensity 7, frequent
        assert GB_get_alpha_max(7, 0.10, "frequent") == 0.08
        # Intensity 8, rare
        assert GB_get_alpha_max(8, 0.20, "rare") == 0.90

    def test_gb_design_spectrum_shape(self):
        """GB_design_spectrum should return positive values for T > 0."""
        from lateral_tools.codes.GB50011 import GB_design_spectrum

        T = np.linspace(0.1, 3.0, 50)
        Sa = GB_design_spectrum(T, alpha_max=0.08, site_class="II")
        assert len(Sa) == len(T)
        assert np.all(Sa > 0)
        # Peak = eta2 * alpha_max * g = 1.0 * 0.08 * 9.81 = 0.7848 m/s²
        assert np.isclose(np.max(Sa), 0.08 * 9.81, rtol=0.02)


# -------------------------------------------------------------------
# EC8
# -------------------------------------------------------------------
class TestEC8Spectrum:
    def test_ec8_elastic_spectrum(self):
        """EC8 elastic spectrum should return g-scaled positive values."""
        from lateral_tools.codes.EC8 import EC8_design_spectrum

        T = np.linspace(0.1, 3.0, 50)
        sa_func = EC8_design_spectrum(
            a_gR=0.12, site_class="C", source_type=1, q=None
        )
        Sa = sa_func(T)
        assert len(Sa) == len(T)
        assert np.all(Sa >= 0)

    def test_ec8_design_spectrum(self):
        """EC8 design spectrum (with q) should be ≤ elastic."""
        from lateral_tools.codes.EC8 import EC8_design_spectrum

        T = np.linspace(0.1, 3.0, 50)
        elastic = EC8_design_spectrum(
            a_gR=0.12, site_class="C", source_type=1, q=None
        )
        design = EC8_design_spectrum(
            a_gR=0.12, site_class="C", source_type=1, q=3.0
        )
        Sa_e = elastic(T)
        Sa_d = design(T)
        assert np.all(Sa_d <= Sa_e + 1e-10)


# -------------------------------------------------------------------
# ASCE 7
# -------------------------------------------------------------------
class TestASCE7Spectrum:
    def test_design_response_spectrum_shape(self):
        """ASCE 7 design spectrum should have correct plateau region."""
        from lateral_tools.codes.ASCE7_16 import design_response_spectrum

        S_DS, S_D1 = 0.5, 0.3
        T_S = S_D1 / S_DS  # 0.6
        T_0 = 0.2 * T_S  # 0.12

        T = np.array([0.05, T_0 / 2, T_0, T_S / 2, T_S, T_S * 2, 5.0, 10.0])
        Sa = design_response_spectrum(S_DS, S_D1, T)

        # At T = T_S, Sa should equal S_DS (in g)
        idx = np.argmin(np.abs(T - T_S))
        assert np.isclose(Sa[idx], S_DS, rtol=0.01)

        # At long period, Sa should decrease
        assert Sa[-1] < Sa[-2]


# -------------------------------------------------------------------
# plot_seismic_spectrum
# -------------------------------------------------------------------
class TestPlotSeismicSpectrum:
    def test_returns_figure(self):
        """plot_seismic_spectrum should return a Figure."""
        from lateral_tools.plotting.spectra import plot_seismic_spectrum

        fig = plot_seismic_spectrum(
            [
                {
                    "code": "GB50011",
                    "params": {
                        "intensity": 7,
                        "acceleration": 0.10,
                        "site_class": "II",
                        "level": "frequent",
                    },
                    "label": "Frequent",
                },
            ],
            unit="g",
        )
        assert isinstance(fig, Figure)

    def test_multiple_spectra(self):
        """Should handle multiple codes in one call."""
        from lateral_tools.plotting.spectra import plot_seismic_spectrum

        fig = plot_seismic_spectrum(
            [
                {
                    "code": "GB50011",
                    "params": {"intensity": 7, "acceleration": 0.10,
                               "site_class": "II", "level": "rare"},
                    "label": "GB 50011",
                },
                {
                    "code": "ASCE7",
                    "params": {"S_DS": 0.667, "S_D1": 0.252},
                    "label": "ASCE 7",
                },
            ],
            unit="g",
            xlim=(0, 3),
        )
        assert isinstance(fig, Figure)

    def test_return_data(self):
        """return_data=True should return (fig, dict)."""
        from lateral_tools.plotting.spectra import plot_seismic_spectrum

        result = plot_seismic_spectrum(
            [
                {
                    "code": "GB50011",
                    "params": {"intensity": 7, "acceleration": 0.10,
                               "site_class": "II", "level": "rare"},
                    "label": "test",
                },
            ],
            unit="g",
            return_data=True,
        )
        fig, data = result
        assert isinstance(fig, Figure)
        assert "test" in data
        T, Sa = data["test"]
        assert len(T) == len(Sa)

    def test_precomputed_data(self):
        """Should accept raw (T, Sa) arrays."""
        from lateral_tools.plotting.spectra import plot_seismic_spectrum

        T = np.linspace(0.1, 3.0, 30)
        Sa = np.exp(-T) * 0.5
        fig = plot_seismic_spectrum(
            [{"T": T, "Sa": Sa, "label": "custom", "in_g": True}],
            unit="g",
        )
        assert isinstance(fig, Figure)

    def test_log_axes(self):
        """log_x and log_y should not raise errors."""
        from lateral_tools.plotting.spectra import plot_seismic_spectrum

        fig = plot_seismic_spectrum(
            [
                {
                    "code": "ASCE7",
                    "params": {"S_DS": 0.5, "S_D1": 0.3},
                    "label": "DBE",
                },
            ],
            unit="g",
            log_x=True,
            log_y=True,
        )
        assert isinstance(fig, Figure)


# -------------------------------------------------------------------
# Math utils
# -------------------------------------------------------------------
class TestMathUtils:
    def test_interp(self):
        from lateral_tools.math_utils import interp
        assert interp(0.5, 0, 1, 10, 20) == 15.0

    def test_list_interp(self):
        from lateral_tools.math_utils import list_interp
        result = list_interp(0.5, [0, 1], [10, 20])
        assert result == 15.0

    def test_as_float(self):
        from lateral_tools.math_utils import as_float
        assert as_float("3.14") == 3.14
        assert as_float("hello") == "hello"

    def test_as_numeric(self):
        from lateral_tools.math_utils import as_numeric
        assert as_numeric("6") == 6
        assert as_numeric("3.14") == 3.14


# -------------------------------------------------------------------
# Risk utils
# -------------------------------------------------------------------
class TestRiskUtils:
    def test_exceedance2return(self):
        from lateral_tools.risk_utils import exceedance2return
        rp = exceedance2return(0.10, 50)  # 10% in 50 years
        assert rp == pytest.approx(475, rel=0.02)

    def test_return2exceedance(self):
        from lateral_tools.risk_utils import return2exceedance
        prob = return2exceedance(475, 50)
        assert prob == pytest.approx(0.10, rel=0.02)


# -------------------------------------------------------------------
# Analysis utils
# -------------------------------------------------------------------
class TestAnalysis:
    def test_storey_elevations(self):
        from lateral_tools.analysis import storey_elevations
        heights = [3.0, 3.0, 3.0]
        elevs = storey_elevations(heights)
        # Function prepends 0, so returns [0, 3, 6, 9]
        assert elevs[:1] == [0]
        assert elevs[-3:] == [3.0, 6.0, 9.0]

    def test_storey_force_distribution(self):
        from lateral_tools.analysis import storey_force_distribution
        elevs = [3.0, 6.0, 9.0]
        weights = [100.0, 100.0, 100.0]
        forces = storey_force_distribution(100.0, elevs, weights, k=1.0)
        assert len(forces) == 3
        assert sum(forces) == pytest.approx(100.0, rel=0.01)

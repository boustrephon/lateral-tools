"""A library for carrying out calculations to IS 1893 (Part 1).

Indian Standard — Criteria for Earthquake Resistant Design of Structures,
Part 1: General Provisions and Buildings.

Covers the issued **2016** edition (Sixth Revision) plus the **2023 draft**
(the draft functions are identified in their docstrings and by their names).

- 2016 functions use the Response Spectrum Method (RSM) and Equivalent
  Static Method (ESM).
- 2023 draft functions use the Linear Dynamic Method (LDM) and Linear
  Equivalent Static Method (LESM).
"""
from __future__ import annotations

# ----------------------------------------------------------------------
# IS 1893-2016 — zone factors and design spectra
# ----------------------------------------------------------------------

Z_dict_2016 = {'II': 0.10, 'III': 0.16, 'IV': 0.24, 'V': 0.36}


def get_Z_2016(zone: str) -> float | None:
    """Zone factor Z for a seismic zone (IS 1893-2016, Table 2).

    Args:
        zone (str): 'II', 'III', 'IV' or 'V'.

    Returns:
        float | None: zone factor, or None if the zone is not recognised.

    Examples:
        >>> get_Z_2016('II')
        0.1
        >>> get_Z_2016('V')
        0.36
    """
    return Z_dict_2016.get(zone)


# Elastic response spectra (Sa/g) for IS 1893-2016.
# Each site-class dict maps a period breakpoint to either a constant or a
# callable giving the ordinate for periods up to that breakpoint.
ESM_spectrum = {
    'A': {0.4: 2.5, 4.0: lambda T: 1 / T, 6.0: 0.25},
    'B': {0.55: 2.5, 4.0: lambda T: 1.36 / T, 6.0: 0.34},
    'C': {0.67: 2.5, 4.0: lambda T: 1.67 / T, 6.0: 0.42},
}

RSM_spectrum = {
    'A': {0.1: lambda T: 1 + 15 * T, 0.4: 2.5, 4.0: lambda T: 1 / T, 6.0: 0.25},
    'B': {0.1: lambda T: 1 + 15 * T, 0.55: 2.5, 4.0: lambda T: 1.36 / T, 6.0: 0.34},
    'C': {0.1: lambda T: 1 + 15 * T, 0.67: 2.5, 4.0: lambda T: 1.67 / T, 6.0: 0.42},
}


def _spectrum_ordinate(T: float, s_dict: dict) -> float | None:
    """Ordinate at period ``T`` from a spectrum given as an ordered
    ``{breakpoint: value-or-callable}`` dict.

    Returns None when ``T`` is beyond the last breakpoint.
    """
    for k, v in s_dict.items():
        if T <= k:
            return v if isinstance(v, (int, float)) else v(T)
    return None


def get_S_a_g(T: float, site_class: str, method: str = 'ESM') -> float | None:
    """Design acceleration coefficient (Sa/g) — IS 1893-2016.

    Args:
        T (float): building period (s).
        site_class (str): soil class ('A', 'B' or 'C').
        method (str, optional): 'RSM' or 'ESM'. Defaults to 'ESM'.

    Returns:
        float | None: design acceleration coefficient (Sa/g), or None if the
        period is beyond the defined range.

    Examples:
        >>> get_S_a_g(0.2, 'A', 'ESM')
        2.5
        >>> round(get_S_a_g(2.0, 'A', 'ESM'), 6)
        0.5
    """
    if method == 'RSM':
        s_dict = RSM_spectrum.get(site_class, {})
    else:
        s_dict = ESM_spectrum.get(site_class, {})
    return _spectrum_ordinate(T, s_dict)


def get_A_H_2016(Z: float, S_a_g: float, I_e: float, R: float) -> float:
    """Design horizontal seismic coefficient Ah — IS 1893-2016.

    Ah = (Z / 2) * (Sa/g) * I_e / R

    Args:
        Z (float): zone factor.
        S_a_g (float): design acceleration coefficient (Sa/g).
        I_e (float): importance factor.
        R (float): response reduction factor.

    Returns:
        float: design horizontal seismic coefficient.
    """
    return 0.5 * Z * S_a_g * I_e / R


rho_dict = {'II': 0.007, 'III': 0.011, 'IV': 0.016, 'V': 0.024}


def get_rho(zone: str) -> float | None:
    """Minimum lateral force coefficient (IS 1893-2016)."""
    return rho_dict.get(zone)


def get_A_H_2016_with_minima(
    T: float,
    zone: str,
    site_class: str,
    I_e: float = 1.0,
    R: float = 1.0,
    method: str = 'ESM',
    with_minima: bool = True,
) -> float:
    """Elastic maximum horizontal PSA for IS 1893-2016, with optional minima.

    Args:
        T (float): building period (s).
        zone (str): 'II', 'III', 'IV' or 'V'.
        site_class (str): 'A', 'B' or 'C'.
        I_e (float, optional): importance factor. Defaults to 1.0.
        R (float, optional): response reduction factor. Defaults to 1.0.
        method (str, optional): 'RSM' or 'ESM'. Defaults to 'ESM'.
        with_minima (bool, optional): apply the minimum force requirement.
            Defaults to True.

    Returns:
        float: elastic maximum horizontal PSA (g).
    """
    Z = get_Z_2016(zone)
    S_a_g = get_S_a_g(T, site_class, method=method)
    A_H = 0.5 * Z * S_a_g * I_e / R
    if with_minima:
        return max(A_H, get_rho(zone))
    return A_H


def get_V_B_2016(
    W: float, T: float, zone: str, site_class: str, I_e: float, R: float,
    method: str = 'ESM',
) -> float:
    """Base shear Vb — IS 1893-2016.

    Vb = max(Ah, rho) * W

    Args:
        W (float): seismic weight.
        T (float): building period (s).
        zone (str): 'II', 'III', 'IV' or 'V'.
        site_class (str): 'A', 'B' or 'C'.
        I_e (float): importance factor.
        R (float): response reduction factor.
        method (str, optional): 'RSM' or 'ESM'. Defaults to 'ESM'.

    Returns:
        float: base shear.
    """
    A_H = get_A_H_2016(get_Z_2016(zone), get_S_a_g(T, site_class, method), I_e, R)
    rho = get_rho(zone)
    return max(A_H, rho) * W


# ----------------------------------------------------------------------
# IS 1893-2023 draft
# ----------------------------------------------------------------------

Table_1_2023 = {
    'Normal Structures': (475, 73),
    'Important Structures': (975, 225),
    'Critical & Lifeline Structures': (2475, 475),
    'Special Structures': (4995, 975),
    'Nuclear Power Plant Structures': (None, None),
}


def get_RP_pairs(category: int | str) -> tuple[str, tuple[int, int] | None]:
    """Return-period pairs for ULS and SLS by structure category (Table 1).

    Args:
        category (int | str): integer 1-4, or one of the category names
            ('Normal Structures', 'Important Structures',
            'Critical & Lifeline Structures', 'Special Structures').

    Returns:
        tuple[str, tuple[int, int] | None]: category name and (ULS, SLS)
        return periods.
    """
    if isinstance(category, int):
        return [(k, v) for i, (k, v) in enumerate(Table_1_2023.items())
                if i + 1 == category][0]
    return category, Table_1_2023.get(category)


Table_2_2023 = {
    'II': {73: 0.0375, 225: 0.05, 475: 0.075, 975: 0.1,
           2475: 0.15, 4975: 0.2, 9975: 0.25},
    'III': {73: 0.075, 225: 0.1, 475: 0.15, 975: 0.2,
            2475: 0.3, 4975: 0.4, 9975: 0.5},
    'IV': {73: 0.18, 225: 0.225, 475: 0.3, 975: 0.36,
           2475: 0.45, 4975: 0.54, 9975: 0.675},
    'V': {73: 0.24, 225: 0.3, 475: 0.4, 975: 0.48,
          2475: 0.6, 4975: 0.75, 9975: 0.9},
    'VI': {73: 0.3, 225: 0.375, 475: 0.5, 975: 0.6,
           2475: 0.75, 4975: 0.9375, 9975: 1.125},
}


def get_Z_2023_dict(zone: str) -> dict:
    """Zone factor Z by return period for a zone (IS 1893-2023 draft)."""
    return Table_2_2023.get(zone, {})


def get_Z_2023(zone: str, category: int | str = 1, LS: str = 'ULS') -> float | None:
    """Zone factor Z for a zone, category and limit state (2023 draft)."""
    name, (ULS, SLS) = get_RP_pairs(category)
    T_RP = ULS if LS == 'ULS' else SLS
    return get_Z_2023_dict(zone).get(T_RP)


Table_3_2023 = {'A': 1500, 'B': 760, 'C': 360, 'D': 180, 'E': 0}


def get_site_class(V_s: float) -> str | None:
    """Site class from the average shear-wave velocity Vs (IS 1893-2023 draft).

    The class is determined from the first (highest) threshold that ``V_s``
    meets, i.e. Table 3 is traversed in descending order of Vs.
    """
    for k, v in Table_3_2023.items():
        if V_s >= v:
            return k
    return None


# Normalised horizontal response spectra for IS 1893-2023 draft.
LDM_spectrum = {
    'A': {0.01: 1.0, 0.1: lambda T: 1 + (50 / 3) * (T - 0.01), 0.4: 2.5,
          6.0: lambda T: 1 / T, 10.0: lambda T: 6 / T ** 2},
    'B': {0.01: 1.0, 0.1: lambda T: 1 + (50 / 3) * (T - 0.01), 0.4: 2.5,
          6.0: lambda T: 1 / T, 10.0: lambda T: 6 / T ** 2},
    'C': {0.01: 1.0, 0.1: lambda T: 1 + (50 / 3) * (T - 0.01), 0.6: 2.5,
          6.0: lambda T: 1.5 / T, 10.0: lambda T: 9 / T ** 2},
    'D': {0.01: 1.0, 0.1: lambda T: 1 + (50 / 3) * (T - 0.01), 0.8: 2.5,
          6.0: lambda T: 2 / T, 10.0: lambda T: 12 / T ** 2},
}

LESM_spectrum = {
    'A': {0.4: 2.5, 6.0: lambda T: 1 / T, 10.0: lambda T: 6 / T ** 2},
    'B': {0.4: 2.5, 6.0: lambda T: 1 / T, 10.0: lambda T: 6 / T ** 2},
    'C': {0.6: 2.5, 6.0: lambda T: 1.5 / T, 10.0: lambda T: 9 / T ** 2},
    'D': {0.8: 2.5, 6.0: lambda T: 2 / T, 10.0: lambda T: 12 / T ** 2},
}


def get_A_NH5(T: float, site_class: str, method: str = 'LESM') -> float | None:
    """Normalised horizontal PSA for 5% damping — IS 1893-2023 draft.

    Args:
        T (float): building period (s).
        site_class (str): soil class ('A', 'B', 'C' or 'D').
        method (str, optional): 'LDM' or 'LESM'. Defaults to 'LESM'.

    Returns:
        float | None: normalised horizontal PSA (A_NH,5%).
    """
    if method == 'LDM':
        s_dict = LDM_spectrum.get(site_class, {})
    else:
        s_dict = LESM_spectrum.get(site_class, {})
    return _spectrum_ordinate(T, s_dict)


def get_A_H(T: float, site_class: str, zone: str, I_e: float = 1.0,
            category: int | str = 1, LS: str = 'ULS',
            method: str = 'LDM') -> float | None:
    """Elastic horizontal PSA Ah (before division by R) — 2023 draft."""
    name, (ULS, SLS) = get_RP_pairs(category)
    T_RP = ULS if LS == 'ULS' else SLS
    Z = get_Z_2023_dict(zone).get(T_RP)
    return Z * I_e * get_A_NH5(T, site_class, method)


def get_A_HD_2023(Z: float, A_NH5: float, I_e: float = 1.0, R: float = 1.0,
                  d_ksi: float = 1.0) -> float:
    """Design horizontal PSA Ahd — IS 1893-2023 draft, Section 6.2.4.

    Ahd = Z * I_e * A_NH,5% * d_ksi / R
    """
    A_H = Z * I_e * A_NH5 * d_ksi
    return A_H / R


def get_A_H_2023_with_minima(
    T: float, zone: str, site_class: str,
    I_e: float = 1.0, R: float = 1.0, d_ksi: float = 1.0,
    category: int | str = 1, LS: str = 'ULS', method: str = 'LESM',
    with_minima: bool = True,
) -> float:
    """Elastic maximum horizontal PSA (g) — IS 1893-2023 draft.

    Args:
        T (float): building period (s).
        zone (str): 'II', 'III', 'IV', 'V' or 'VI'.
        site_class (str): 'A', 'B', 'C' or 'D'.
        I_e (float, optional): importance factor. Defaults to 1.0.
        R (float, optional): response reduction factor. Defaults to 1.0.
        d_ksi (float, optional): damping modification factor. Defaults to 1.0.
        category (int | str, optional): structure category. Defaults to 1.
        LS (str, optional): limit state ('ULS' or 'SLS'). Defaults to 'ULS'.
        method (str, optional): 'LDM' or 'LESM'. Defaults to 'LESM'.
        with_minima (bool, optional): apply minimum force requirements.
            Defaults to True.

    Returns:
        float: elastic maximum horizontal PSA (g).
    """
    name, (ULS, SLS) = get_RP_pairs(category)
    T_RP = ULS if LS == 'ULS' else SLS
    Z = get_Z_2023_dict(zone).get(T_RP)
    A_NH = get_A_NH5(T, site_class, method) * d_ksi
    A_H = Z * I_e * A_NH
    A_HD = A_H / R
    A_H_min = 0.625 * Z * I_e / R
    if with_minima:
        return max(A_HD, A_H_min, 0.015)
    return A_HD


def get_V_BDH(W: float, T: float, zone: str, site_class: str,
              I_e: float = 1.0, R: float = 1.0, category: int | str = 1,
              LS: str = 'ULS', method: str = 'LESM') -> float:
    """Design base shear Vb — IS 1893-2023 draft (with minima)."""
    name, (ULS, SLS) = get_RP_pairs(category)
    T_RP = ULS if LS == 'ULS' else SLS
    Z = get_Z_2023_dict(zone).get(T_RP)
    A_H = Z * I_e * get_A_NH5(T, site_class, method)
    A_HD = A_H / R
    A_H_min = 0.625 * Z * I_e / R
    return max(A_HD, A_H_min, 0.015) * W


def V_BDH_min(Z: float, I_e: float, R: float, W: float) -> tuple[float, float]:
    """Minimum design base shear terms — IS 1893-2023 draft."""
    return 0.625 * Z * I_e / R * W, 0.015 * W


def get_d_ksi(ksi: float, T: float, site_class: str) -> float:
    """Damping modification factor ``d_ksi`` — IS 1893-2023 draft.

    Modifies the 5%-damped spectrum for other damping ratios: it is 1.0 up to
    T = 0.01 s, ramps up to ``eta_max`` at T = 0.1 s, stays at ``eta_max``
    over the acceleration plateau (up to ``T_AV``), then decays back to 1.0
    at T = 6.0 s.

    Args:
        ksi (float): damping ratio as a fraction (e.g. 0.05 for 5%).
        T (float): building period (s).
        site_class (str): 'A', 'B', 'C' or 'D'.

    Returns:
        float: damping modification factor.

    Raises:
        ValueError: if ``site_class`` is not 'A'-'D' or ``ksi`` >= 30%.
    """
    T_AV = {'A': 0.4, 'B': 0.4, 'C': 0.6, 'D': 0.8}.get(site_class, 0)
    if T_AV == 0:
        raise ValueError('site_class should be one of A, B, C or D')

    if ksi < 0.0:
        eta_max = 3.2
    elif ksi < 0.005:
        eta_max = 3.2 + ksi * (((7 / 2.005) ** 0.6 - 3.2) / 0.005)
    elif ksi <= 0.05:
        eta_max = (7 / (2 + ksi)) ** 0.6
    elif ksi <= 0.30:
        eta_max = (10 / (5 + ksi)) ** 0.5
    else:
        raise ValueError('ksi should be less than 30%')

    if T <= 0.01:
        return 1.0
    elif T <= 0.1:
        return 1.0 + (T - 0.01) / 0.09 * (eta_max - 1.0)
    elif T <= T_AV:
        return eta_max
    elif T <= 6.0:
        return eta_max - (T - T_AV) / (6.0 - T_AV) * (eta_max - 1.0)
    else:
        return 1.0


# ----------------------------------------------------------------------
# Fundamental natural periods (IS 1893-2016)
# ----------------------------------------------------------------------

T_a_Cd = (0.075, 0.080, 0.085)  # (steel, RC, steel-braced) frame types


def get_T_a(h: float, n: int) -> float:
    """Approximate fundamental natural period Ta for frame buildings.

    Ta = Cd * h^0.75

    Args:
        h (float): height of the building (m).
        n (int): frame type index — 0 steel MRF, 1 RC MRF, 2 steel braced.

    Returns:
        float: approximate fundamental natural period (s).
    """
    return T_a_Cd[n] * h ** 0.75


def get_T_aw(h: float, d: float, Aw_list, Lw_list) -> float:
    """Approximate fundamental natural period for RC shear-wall buildings.

    Args:
        h (float): height of the building (m).
        d (float): base dimension parallel to the applied force (m).
        Aw_list: list of effective wall areas (m^2).
        Lw_list: list of wall lengths (m).

    Returns:
        float: approximate fundamental natural period (s).
    """
    A_w = sum(A * (0.2 + min(0.9, L / h) ** 2) for A, L in zip(Aw_list, Lw_list))
    return min(0.075 * h ** 0.75 * A_w ** -0.5, 0.09 * d ** -0.5)


def T_aa(h: float, d: float) -> float:
    """Alternative approximate fundamental natural period (IS 1893-2016).

    Taa = 0.09 * h / sqrt(d)
    """
    return 0.09 * h * d ** -0.5

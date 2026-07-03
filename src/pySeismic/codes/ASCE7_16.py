"""
A library for carrying out calculations to ASCE7-16

"This Standard provides requirements for 
general structural design and includes 
means for determining dead, live, soil, 
flood, wind, snow, rain, atmospheric 
ice, and earthquake loads, as well as 
their combinations, which are suitable 
for inclusion in building codes and 
other documents."

Previous editions are:
* ASCE 7-10
* ASCE 7-05
* ASCE 7-02
* ASCE 7-98
* ASCE 7-95
* ASCE 7-93
* ASCE 7-88
"""

from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10

import numpy as np
# from itertools import accumulate
import doctest
from sys import path as sys_path
from os.path import abspath, join as path_join, dirname

# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

from pySeismic.utilities import list_interp

#### Seismic Hazard ####

# Functions
# T_code_get

# F_a and F_v taken from ASCE7-16, tables 11.4-1 & 11.4-2
F_a_headers = (0.25,0.5,0.75,1.0,1.25,1.5)
F_a_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
           'B': (0.9,0.9,0.9,0.9,0.9,0.9),
           'C': (1.3,1.3,1.2,1.2,1.2,1.2),
           'D': (1.6,1.4,1.2,1.1,1.0,1.0),
           'E': (2.4,1.7,1.3,None,None,None),
           'F': (None,None,None,None,None,None)}

F_v_headers = (0.1,0.2,0.3,0.4,0.5,0.6)
F_v_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
           'B': (0.8,0.8,0.8,0.8,0.8,0.8),
           'C': (1.5,1.5,1.5,1.5,1.5,1.4),
           'D': (2.4,2.2,2.0,1.9,1.8,1.7),
           'E': (4.2,None,None,None,None,None),
           'F': (None,None,None,None,None,None)}


def site_coeffs(site_class:str, S_s:float, S_1:float, 
                debug:bool=False)->tuple[float|None,float|None]:
    """Returns F_a, F_v when provided with 
    site class (A-F) and S_s and S_1 values (Tables 11.4-1 & 2, ASCE 7-16)
    """
    F_a_list = F_a_dict[site_class]
    F_v_list = F_v_dict[site_class]
    if debug:
        print(f'site_class = {site_class}; S_s = {S_s}; S_1 = {S_1}; F_a = {F_a_list}; F_v = {F_v_list}')
    return list_interp(S_s, F_a_headers, F_a_list), list_interp(S_1, F_v_headers, F_v_list)


def calculate_k_factor(T: float) -> float:
    """
    Calculate the period-dependent exponent k for vertical force distribution.
    Per ASCE 7-16 Section 12.8.3:
    - k = 1 for T <= 0.5 s
    - k = 2 for T >= 2.5 s
    - Linear interpolation for intermediate periods
    """
    if T <= 0.5:
        return 1.0
    elif T >= 2.5:
        return 2.0
    else:
        # Linear interpolation
        return 1.0 + (T - 0.5) * (1.0 / 2.0)

def S_MS_get(F_a, S_s):
    """Eq 11.4-1 ASCE7-16"""
    S_MS = F_a * S_s
    return S_MS


def S_M1_get(F_v, S_1):
    """Eq 11.4-2 ASCE7-16"""
    S_M1 = F_v * S_1
    return S_M1


def S_DS_get(S_MS):
    """Eq 11.4-3 ASCE7-16"""
    S_DS = 2 / 3 * S_MS
    return S_DS


def S_D1_get(S_M1):
    """Eq 11.4-4 ASCE7-16"""
    S_D1 = 2 / 3 * S_M1
    return S_D1

def site_accels_D(site_class:str, S_s:float, S_1:float, debug:bool=False)->tuple[float,float]:
    """Returns site spectral accelerations for D class sites
    Eq 11.4-3 & 4 ASCE7-16

    Args:
        site_class (str): Site class
        S_s (float): Short period spectral acceleration
        S_1 (float): 1-second spectral acceleration
        debug (bool, optional): Debug flag. Defaults to False.

    Raises:
        ValueError: _description_

    Returns:
        tuple[float,float]: _description_
    """
    F_a, F_v = site_coeffs(site_class, S_s, S_1)
    if F_a is None or F_v is None:
        raise ValueError(f'Invalid site class {site_class}')
    S_DS = 2 / 3 * F_a * S_s
    S_D1 = 2 / 3 * F_v * S_1
    if debug:
        print(f'site_class = {site_class}')
        print(f'S_s = {S_s}; S_1 = {S_1}')
        print(f'S_DS = {S_DS}; S_D1 = {S_D1}')
    return S_DS, S_D1
    

def seismic_design_category_get(S_s, S_1, site_class, risk_category, debug=False):
    """Returns the ASCE7-10 seismic design category based on the 
    site spectral accelerations, site class and risk category 
    Table 11.6-1 - Short Period Response
    Table 11.6-2 - 1-sec Period Response
    >>> seismic_design_category_get(0.75, 0.2, 'D', 'II')
    'D'
    >>> seismic_design_category_get(0.3, 0.15, 'C', 'II')
    'C'
    >>> seismic_design_category_get(0.18, 0.09, 'B', 'II')
    'A'
    """
    S_DS, S_D1 = site_accels_D(site_class, S_s, S_1, debug=debug)

    accel_limits_S = [0.0, 0.167, 0.33, 0.50]
    accel_limits_1 = [0.0, 0.067, 0.133, 0.20]
    
    if risk_category == 'IV':
        design_cat_list = ['A', 'C', 'D', 'D']
    else:
        design_cat_list = ['A', 'B', 'C', 'D']
    
    seismic_design_category_S = 'A'
    for S_lim, S_val in zip(accel_limits_S, design_cat_list):
        if S_DS >= S_lim:
            seismic_design_category_S = S_val
    
    seismic_design_category_1 = 'A'
    for S_lim, S_val in zip(accel_limits_1, design_cat_list):
        if S_D1 >= S_lim:
            seismic_design_category_1 = S_val
    
    seismic_design_category = max(seismic_design_category_S, seismic_design_category_1)

    if debug:
        print(f'SDC = max({seismic_design_category_S}, {seismic_design_category_1}) = {seismic_design_category}')

    if risk_category in ['I', 'II', 'III'] and S_1 >= 0.75:
        seismic_design_category = 'E'
    
    if risk_category == 'IV' and S_1 >= 0.75:
        seismic_design_category = 'F'
    
    return seismic_design_category

def design_response_spectrum(
    S_DS: float,
    S_D1: float,
    T: np.ndarray,
    T_L: float = 12.0,
) -> np.ndarray:
    """ASCE 7-16 design response spectrum (Section 11.4.5, Figure 11.4-1).

    Parameters
    ----------
    S_DS : float
        Design spectral acceleration at short periods (g),
        from Eq. 11.4-3.
    S_D1 : float
        Design spectral acceleration at 1-second period (g),
        from Eq. 11.4-4.
    T : np.ndarray
        Array of periods (s) at which to evaluate the spectrum.
    T_L : float, optional
        Long-period transition period (s) from Figures 22-12 to 22-16.
        Default 12.0.

    Returns
    -------
    np.ndarray
        Spectral acceleration Sa (g) at each period in *T*.

    Notes
    -----
    The spectral shape follows:

    * T < T_0:   Sa = S_DS * (0.4 + 0.6 * T / T_0)
    * T_0 ≤ T ≤ T_S:  Sa = S_DS
    * T_S < T ≤ T_L:  Sa = S_D1 / T
    * T > T_L:   Sa = S_D1 * T_L / T²

    where *T_0* = 0.2 * S_D1 / S_DS and *T_S* = S_D1 / S_DS.

    For the MCE_R spectrum, pass 1.5 × the design S_DS and S_D1 values.
    The returned values are in g — multiply by 9.81 for m/s² if needed.
    """
    T_S = S_D1 / S_DS
    T_0 = 0.2 * T_S

    Sa = np.where(
        T < T_0,
        S_DS * (0.4 + 0.6 * T / T_0),
        np.where(
            T <= T_S,
            S_DS,
            np.where(
                T <= T_L,
                S_D1 / T,
                S_D1 * T_L / T ** 2,
            ),
        ),
    )
    return Sa


def main():
    S_s = 2.5 * 0.12
    S_1 = 1.0 * 0.12
    site_class = 'B'
    risk_category = 'IV'
    # S_DS, S_D1 = site_accels_D(site_class, S_s, S_1, debug=True)
    sdc = seismic_design_category_get(S_s, S_1, site_class, risk_category, debug=True)
    print(f'Seismic Design Category: {sdc}')


if __name__ == "__main__":
    doctest.testmod()
    print('== Doctest Complete ==')
    main()

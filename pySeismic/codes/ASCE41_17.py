"""A library for carrying out calculations to ASCE41-17

ASCE 41-17 updates ASCE 41-13 which merges ASCE 31-03 (evaluation) with 
ASCE 41-06 (rehabilitation and retrofit)

* ASCE/SEI 41-06, Seismic Rehabilitation of Existing Buildings
* ASCE/SEI 31-03, Seismic Evaluation of Existing Buildings.

The predecessor to ASCE 31 was FEMA 310, and the 
predecessor to ASCE 41 was FEMA 356 (in turn based
on FEMA 273).

Note that all references to ASCE 7 are related generically 
to ASCE 7-16

The full list of references may be found in section C1.3.

* ASCE 41-13
* ASCE 41-06

TODO:
* Fix the M_s calculations

"""

from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
from itertools import accumulate

from pySeismic.utilities import printt, interp, list_interp, \
    storey_force_coeffs, storey_force_distribution, storey_forces, \
    storey_shear_moment, shear_coords
from pySeismic.importing import k_gen # Eq 4-2a
from pySeismic.codes.ASCE7_16 import site_coeffs as ASCE7_16_site_coeffs


#=== The following functions may be found in the utilities module ===#
# def printt(table, h1 = 20, h2 = 12, h3 = 4):
#     """
#     Printing tables where first column is text and the remainder are floats
#     h1 is the width of the text column, 
#     h2 is the width of all other columns
#     h3 is the number of decimal places
#     """
#     for line in table:
#         print(f'{line[0]:{h1}s}', end = ' ')
#         [print(f'{item:{h2-h3}.{h3}f}', end = ' ') for item in line[1:]]
#     return None


# def interp(x, x1, x2, y1, y2):
#     """Returns an interpolated y-value for a line between (x1, y1) and (x2, y2)
#     for a given x-value"""
#     if y1 is None or y2 is None:
#         return None
#     elif x1 == x2:
#         return 0.5 * (y1 + y2)
#     else:
#         return y1 + (x - x1) / (x2 - x1) * (y2 - y1)


# def list_interp(val, list_1, list_2):
#     """
#     Returns an interpolated value based on a list of reference values and a list of target values
#     val is the lookup value
#     list_1 is the lookup list
#     list_2 is the result list with values corresponding to those on list_1
#     """
#     i_list = [i for i, n in enumerate(list_1) if n == val]
#     if i_list:
#         return list_2[i_list[0]]
#     elif val <= list_1[0]:
#         return list_2[0]
#     elif val >= list_1[-1]:
#         return list_2[-1]
#     else:
#         index_list = [i for i, (x1, x2) in enumerate(zip(list_1[:-1], list_1[1:])) if val >= x1 and val <= x2]
#         if len(index_list) == 1:
#             j = index_list[0]
#             return interp(val, list_1[j], list_1[j+1], list_2[j], list_2[j+1])
#         else:
#             return None
#
#
# def storey_force_coeffs(elev_list, weight_list, k = 1.0):
#     """
#         Inputs:
#             elev_list - list of storey elevations (m)
#           weight_list - list of storey weights (kN)
#     """
#     divisor = sum(w * h**k for w, h in zip(weight_list, elev_list))
#     return [w * h**k / divisor for w, h in zip(weight_list, elev_list)]


# def storey_force_distribution(V, elev_list, weight_list, k = 1.0):
#     """
#         Inputs:
#             V  - base shear in kN
#             elev_list - list of storey elevations (m)
#           weight_list - list of storey weights (kN)
#     """
#     divisor = sum(w * h**k for w, h in zip(weight_list, elev_list))
#     return [w * h**k * V / divisor for w, h in zip(weight_list, elev_list)]


# def storey_forces(force_list, rev = True):
#     """Returns a list of accumulated storey forces - i.e. storey shears
#     NB lists have to be reversed before applying the accumulate function"""
#     shear_list_rev = accumulate(force_list[::-1] if rev else force_list)
#     return list(shear_list_rev)[::-1] if rev else list(shear_list_rev)


# def storey_shear_moment(force_list, height_list, rev = True):
#     """Returns a tuple of lists containing storey shear and storey moment
#     (from ground up, increasing z if 'rev' is True)"""
#     shear_list_rev = list(accumulate(force_list[::-1] if rev else force_list))
#     v_h_rev = [v * h for v, h in zip(shear_list_rev, height_list[::-1] if rev else height_list)]
#     mom_list_rev = list(accumulate([0]+v_h_rev[:-1]))
#     return shear_list_rev[::-1] if rev else shear_list_rev, mom_list_rev[::-1] if rev else mom_list_rev


# def shear_coords(elevs, shears):
#     """Returns a list of (z,V) tuples for plotting storey shear diagrams
#     Inputs are lists of elevations and shear values
#     NB. Interleaving is required because of the steps in the shear diagram"""
#     l_1 = list(zip(elevs[1:],shears[1:]))
#     l_2 = list(zip(elevs,shears[1:] + [0.0]))
#     out = []
#     while l_2:
#         out.append(l_2.pop()) if l_2 else None
#         out.append(l_1.pop()) if l_1 else None
#     return out[::-1]


# -----------------------------------------------------------------
# ----- Seismic Hazard --------------------------------------------
# -----------------------------------------------------------------

# T_code_get, B_1_func, site_coeffs, code_period, S_a_func, 

# Eq 4-5 & 7-18 ASCE41-17
# `S2` & `S2a refers to eccentric braced frames
#   'S2b': 0.030, 
C_t_dict = {'S1':0.035, 'S1a':0.035, 'S2':0.030, 'S2a': 0.030, 
            'C1':0.018, 'other':0.020}

def C_t_get(x):
    """Empirical building period parameter (C_t) for imperial units of ft (default of 0.020)
    Eq 4-5 & 7-18 ASCE41"""
    return C_t_dict.get(x, 0.020)  

# Eq 4-5 & 7-18 ASCE41
beta_dict = {'S1':0.80, 'S1a':0.80, 'C1':0.90, 'other':0.75}

def beta_get(x):
    """Empirical building period parameter, x (default of 0.75) Eq 4-5 & 7-18 ASCE41"""
    return beta_dict.get(x, 0.75)  

unit_l_dict = {'m': 1.0, 'cm': 0.01, 'mm':0.001, 'in':0.0254, 'ft':0.3048, 'yd':0.9144}

def unit_l_get(x):
    """a units scaling lookup function that returns default of None"""
    return unit_l_dict.get(x, None)  


def T_code_get(s_type: str, height: float, units: str = 'm') -> float:
    """Returns the code period in seconds for a given structure type and height 
    (Eq 4-4) & (Eq 7-18) ASCE41
    
    Args:
        s_type: The building type (e.g. 'S1' or 'C2a')
            :latex: \\type
        height: The building height
            :units: m
            :latex: \\S_S
        units: The height units (e.g. 'm')
            :latex: \\units
    
    Returns:
        T_code: The building period
            :units: s
            :latex: \\T_code
    
    .. codeauthor:: Andrew Mole <andrew.mole@arup.com>
    
    .. levelofreview:: Work in progress

    """
    C_t = C_t_get(s_type)
    beta = beta_get(s_type)
    u_fac = unit_l_get(units)
    ft_fac = unit_l_get('ft')
    return C_t * height ** beta if units == 'ft' else C_t * (u_fac / ft_fac * height) ** beta


# def code_period(h_n, s_type):
#     """Returns period in seconds when given height (h_n) in metres and ASCE41-17 structural type (e.g. 'C1')"""
#     return C_t_get(s_type) * (h_n / 0.3048) ** beta_get(s_type)


# Table 2-4 Level of Seismicity Definitions
def seismicity_level_get(S_DS, S_D1):
    """Returns level of seismicity, Table 2-4, ASCE41-17"""
    S_DS_line = (0.0, 0.167, 0.33, 0.5)
    S_D1_line = (0.0, 0.067, 0.133, 0.2)
    seismicity_levels = ('Very Low', 'Low', 'Moderate', 'High')
    i_1 = [i for i, s_ds in enumerate(S_DS_line) if s_ds <= S_DS][-1]
    i_2 = [i for i, s_d1 in enumerate(S_D1_line) if s_d1 <= S_D1][-1]
    return seismicity_levels[max(i_1, i_2)]


# F_a and F_v taken from ASCE7-16, tables 11.4-1 & 11.4-2
# F_a_headers = (0.25,0.5,0.75,1.0,1.25,1.5)
# F_a_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
#            'B': (0.9,0.9,0.9,0.9,0.9,0.9),
#            'C': (1.3,1.3,1.2,1.2,1.2,1.2),
#            'D': (1.6,1.4,1.2,1.1,1.0,1.0),
#            'E': (2.4,1.7,1.3,None,None,None),
#            'F': (None,None,None,None,None,None)}

# F_v_headers = (0.1,0.2,0.3,0.4,0.5,0.6)
# F_v_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
#            'B': (0.8,0.8,0.8,0.8,0.8,0.8),
#            'C': (1.5,1.5,1.5,1.5,1.5,1.4),
#            'D': (2.4,2.2,2.0,1.9,1.8,1.7),
#            'E': (4.2,None,None,None,None,None),
#            'F': (None,None,None,None,None,None)} 


# def site_coeffs(site_class, S_s, S_1):
#     """Returns F_a, F_v when provided with site class (A-F) and S_s and S_1 values"""
#     F_a_list = F_a_dict[site_class]
#     F_v_list = F_v_dict[site_class]
#     return list_interp(S_s, F_a_headers, F_a_list), list_interp(S_1, F_v_headers, F_v_list)


def B_1_func(damp_ratio):
    """Spectral response damping coefficient from 
    damping ratio (not percentage) Eq 2-11 ASCE41"""
    return 4.0 / (5.6 - log(100 * damp_ratio)) 


def S_a_func(T, S_s, S_1, site_class, T_L = 12.0, damp = 0.05):
    """Returns the spectral acceleration when provided 
    with period (T), S_s, S_1 and site class, Figure 2-1, ASCE41-17
    """
    F_a, F_v = ASCE7_16_site_coeffs(site_class, S_s, S_1)
    S_xs = F_a * S_s  # Eq 2-1
    S_x1 = F_v * S_1  # Eq 2-2
    T_s = S_x1 / S_xs
    # T_0 = 0.2 * T_s
    B_1 = B_1_func(damp)
    return min(S_xs / B_1, 
               S_xs * ((5.0 / B_1 - 2.0) * T / T_s + 0.4), 
               S_x1 / B_1 / T) if T < T_L else T_L * S_x1 / (B_1 * T**2.0)


# -----------------------------------------------------------------
# ----- Tier One Data and Functions -------------------------------
# -----------------------------------------------------------------

# M_s_func, v_j_avg_mf, D_r, v_j_avg_sw


# Simplified Rehabilitation
# Table 3-4  ASCE41-17 Limitations on Use of the Tier 1 & Tier 2 Procedures (number of storeys)
# No Limit (NL) = -1, Not Permitted (NP) = 0
struct_types_full_dict = {
'W1': {'Description':'Light Wood Frames', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):4, 
       ('Moderate', 'S-5'):4, ('Moderate', 'S-1'):4, ('High', 'S-5'):4, ('High', 'S-1'):4},
'W1a': {'Description':'Multistory Multi-Unit Residential Wood Frames', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):6, ('High', 'S-5'):6, ('High', 'S-1'):4},
'W2': {'Description':'Commercial and Industrial Wood Frames', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):6, ('High', 'S-5'):6, ('High', 'S-1'):4},

'S1': {'Description':'Steel MRF Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'S1a': {'Description':'Steel MRF Flexible diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'S2': {'Description':'Steel Braced Frames Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-5'):8, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'S2a': {'Description':'Steel Braced Frames Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-5'):8, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'S3': {'Description':'Metal Building Frames', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):1, ('Low', 'S-5'):1, ('Low', 'S-1'):1, 
       ('Moderate', 'S-5'):1, ('Moderate', 'S-1'):1, ('High', 'S-5'):1, ('High', 'S-1'):1},
'S4': {'Description':'Dual Systems with Backup Steel MRF', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'S5': {'Description':'Steel Frame Infill URM Rigid Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):4},
'S5a': {'Description':'Steel Frame Infill URM Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):4},
'S6': {'Description':'Steel Plate Shear Wall', 
       ('Very Low', 'S-5'):0, ('Very Low', 'S-1'):0, ('Low', 'S-5'):0, ('Low', 'S-1'):0, 
       ('Moderate', 'S-5'):0, ('Moderate', 'S-1'):0, ('High', 'S-5'):0, ('High', 'S-1'):0},

'CFS1': {'Description':'Cold-Formed Steel Light-Frame Construction, Shear Wall', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):6, ('High', 'S-5'):6, ('High', 'S-1'):4},
'CFS2': {'Description':'Cold-Formed Steel Light-Frame Construction, Strap-braced wall', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):6, ('High', 'S-5'):6, ('High', 'S-1'):4},
'C1': {'Description':'RC MRF', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'C2': {'Description':'RC Shear Wall Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'C2a': {'Description':'RC Shear Wall Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'C3': {'Description':'RC Frame URM Infill Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):4},
'C3a': {'Description':'RC Frame URM Infill Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-5'):12, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):4},
'PC1': {'Description':'Precast/Tilt-Up RC SW Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):3, ('Low', 'S-1'):2, 
       ('Moderate', 'S-5'):2, ('Moderate', 'S-1'):2, ('High', 'S-5'):2, ('High', 'S-1'):2},
'PC1a': {'Description':'Precast/Tilt-Up RC SW Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):3, ('Low', 'S-1'):2, 
       ('Moderate', 'S-5'):2, ('Moderate', 'S-1'):2, ('High', 'S-5'):2, ('High', 'S-1'):2},
'PC2': {'Description':'Precast RC Frame with Shear Walls', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):0, ('High', 'S-5'):4, ('High', 'S-1'):0},
'PC2a': {'Description':'Precast RC Frame without Shear Walls', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):0, ('High', 'S-5'):4, ('High', 'S-1'):0},

'RM1': {'Description':'Reinforced Masonry Bearing Walls Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-5'):8, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},
'RM2': {'Description':'Reinforced Masonry Bearing Walls Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-5'):8, ('Moderate', 'S-1'):8, ('High', 'S-5'):8, ('High', 'S-1'):6},

'URM': {'Description':'Unreinforced Masonry Bearing Walls Rigid diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):6, ('Low', 'S-1'):4, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):0, ('High', 'S-5'):4, ('High', 'S-1'):0},
'URMa': {'Description':'Unreinforced Masonry Bearing Walls Flexible Diaphragm', 
       ('Very Low', 'S-5'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-5'):6, ('Low', 'S-1'):4, 
       ('Moderate', 'S-5'):6, ('Moderate', 'S-1'):0, ('High', 'S-5'):4, ('High', 'S-1'):0},

'SI_PD': {'Description':'Seismic Isolation or Passive Dissipation', 
       ('Very Low', 'S-5'):0, ('Very Low', 'S-1'):0, ('Low', 'S-5'):0, ('Low', 'S-1'):0, 
       ('Moderate', 'S-5'):0, ('Moderate', 'S-1'):0, ('High', 'S-5'):0, ('High', 'S-1'):0},
       }


def S_a_Tier_1_func(S_s: float, S_1: float, T: float, site_class: str) -> float:
    """Returns spectral acceleration for Tier 1 Assessment

    4.4.2.3 Spectral Acceleration ASCE41-17

    Args:
        S_s: The short period site spectral acceleration
            :units: g
            :latex: \\S_S
        S_1: The 1-second site spectral acceleration
            :units: g
            :latex: \\S_1
        T: Building period
            :units: s
        site_class: The site classification (A-E)
            :latex: \\class
    
    Returns:
        S_a: The building spectral acceleration
            :units: g
            :latex: \\S_a
    
    .. codeauthor:: Andrew Mole <andrew.mole@arup.com>
    
    .. levelofreview:: Work in progress
    """
    F_a, F_v = ASCE7_16_site_coeffs(site_class, S_s, S_1)
    S_XS = F_a * S_s  # Eq 2-1
    S_X1 = F_v * S_1  # Eq 2-2
    return min(S_XS, S_X1 / T)  # Eq 4-3


# ASCE41-17 Table 4-2 Default compressive strengths (f'c) in ksi
table_4_2_dict_17 = {
    (1900, 1919): {'Beams': 2.0, 'Slabs and Columns': 1.5, 'Walls': 1.0},
    (1920, 1949): {'Beams': 2.0, 'Slabs and Columns': 2.0, 'Walls': 2.0},
    (1950, 1969): {'Beams': 3.0, 'Slabs and Columns': 3.0, 'Walls': 2.5},
    (1970, 2100): {'Beams': 3.0, 'Slabs and Columns': 3.0, 'Walls': 3.0},
}


# Table 4-3. Default Yield Strengths (fy) of Reinforcing Steel (ksi)
# The dictionary keys are the start & end dates as a tuple,
# and the sub-keys are the steel grade in ksi
# The terms `structural` (33), `intermediate` (40), and `hard (50)` 
#   became obsolete in 1968.
table_4_3_dict_17 = {
             (1911, 1959): {33:True, 40:True, 50:True, 60:False, 65:True, 70:False, 75:False},
             (1959, 1966): {33:True, 40:True, 50:True, 60:True, 65:True, 70:True, 75:True},
             (1966, 1987): {33:False, 40:True, 50:True, 60:True, 65:True, 70:True, 75:False},
             (1987, 3000): {33:False, 40:True, 50:True, 60:True, 65:True, 70:True, 75:True},
             }


# Table 4-4. Default Yield Strengths (Fy) of Archaic Materials (ksi)
table_4_4_dict_17 = {('Pre-1900', 'Steel'): 24}


# ASCE41-17 Table 4-7 Modification Factor, C
# For number of storeys
table_4_7_dict_17 = {
    ('W1', 'W1a', 'W2', 'CFS1'): {1: 1.3, 2: 1.1, 3: 1.0, 4: 1.0},
    ('S1', 'S3', 'C1', 'PC2a'): {1: 1.3, 2: 1.1, 3: 1.0, 4: 1.0},
    ('S4', 'S5', 'C2', 'C3', 'PC1a', 'PC2', 'RM2', 'URM'): {1: 1.4, 2: 1.2, 3: 1.1, 4: 1.0},
    ('S2', ): {1: 1.4, 2: 1.2, 3: 1.1, 4: 1.0},
    ('CFS2', ): {1: 1.4, 2: 1.2, 3: 1.1, 4: 1.0},
    ('URM', ) : {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0},
    ('S1a', 'S2a', 'S5a', 'C2a', 'C3a', 'PC1', 'RM1') : {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0},
}


def table_4_7_dict_17_get(building_type, n_storeys):
    """Return Modification Factor, C (Table 4-7) ASCE41-17 
    For building type and number of storeys"""
    for k, v in table_4_7_dict_17.items():
        if building_type in k:
            return v.get(min(n_storeys, 4), None)
        else:
            return None


# ASCE41-17 Table 4-8 Ms Factors for Shear Walls
M_s_4_8_dict = {'C1': {'CP':4.5, 'LS': 3.0, 'IO': 1.5},
            'URM': {'CP':1.75, 'LS': 1.25, 'IO': 1.0}}


M_s_4_9_dict = {
    'Tube': {'CP': (7.0,3.5), 'LS': (4.5,2.5), 'IO': (2.0,1.25), 'func': lambda x: x**(-0.5), 'lims': (90,190)},
    'Pipe': {'CP': (7.0,3.5), 'LS': (4.5,2.5), 'IO': (2.0,1.25), 'func': lambda x: x**(-1), 'lims': (1500,6000)},
    'Tension-only': {'CP': (3.5,3.5), 'LS': (2.5,2.5), 'IO': (1.25,1.25)},
    'Cold-formed steel strap-braced wall': {'CP': (3.5,3.5), 'LS': (2.5,2.5), 'IO': (1.25,1.25)},
    'Others': {'CP': (7.0,7.0), 'LS': (4.5,4.5), 'IO': (2.0,2.0)},
}


# M_s_table_dict = {
#     '4_8'
# }


def table_4_9_func(brace_type, d_t, F_y, perf):
    M_s_dict = M_s_4_9_dict.get(brace_type, M_s_4_9_dict.get('Others'))
    func = M_s_dict.get('func', None)
    lims = M_s_dict.get('lims', None)
    if lims is None:
        return M_s_dict.get(perf)[0]
    else:
        F_ye = 1.25 * F_y
        FF = func(F_ye)
        limits = (lims[0] * FF, lims[1] * FF)
        if d_t <= limits[0]:
            return M_s_dict.get(perf)[0]
        elif d_t >= limits[1]:
            return M_s_dict.get(perf)[1]
        else:
            M_s = M_s_dict.get(perf)
            return interp(d_t, limits[0], limits[1], M_s[0], M_s[1])


# def M_s_func(building_system, table, perf):
#     """Returns data from a dictionary of M_s (Tier 1 system modification factor,
#     ch 4, ASCE41-17)
    
#     inputs:
#         name - the name of the system (currently a choice of 'C1', 'C2', 'URM')
#         perf - performance level:
#             'IO' - Immediate Occupancy
#             'LS' - Life Safety
#             'CP' - Collapse Prevention
#     """
#     return M_s_dict.get(name).get(perf)



def v_j_avg_mf(n_c, n_f, V_j, A_c, M_s):
    """Returns 
    Inputs:
        n_c - number of columns
        n_f - number of frames in direction of loading
        V_j - storey shear (sect 4.4.2.2)
        A_c - sum of cross-sectional area of """
    return n_c/(n_c - n_f) * V_j / A_c / M_s


def D_r(k_b, k_c, h, E, V_c, base_is_pinned=False):
    """Story drift calculation for moment frames(Eq 4-6)
    k_b = I_b / L_b
    k_c = I_c / h_c
    h story height (in)
    I moment of inertia (in^4)
    E modulus of elasticity (ksi)
    V_c column shear (kip)
    base_is_pinned
    NB This should be doubled if base of columns is pinned
    """
    base_factor = 2 if base_is_pinned else 1
    return V_c * (k_b + k_c) / k_b / k_c * h / 12 / E * base_factor


def v_j_avg_rcmf(V_j, A_c, M_s, n_c, n_f, perf=None):
    """Shear stress in concrete frame columns (Eq 4-7) ASCE41-17
    n_c = Total number of columns;
    n_f = Total number of frames in the direction of loading;
    A_c = Summation of the cross-sectional area of all columns in the
    story under consideration;
    V_j = Story shear computed in accordance with Section 4.4.2.2; and
    M_s = System modification factor from dictionary below"""
    if perf is not None:
        M_s_dict = {'CP':2.0, 'LS':1.5, 'IO':1.0}
        M_s = M_s_dict.get(perf,1.0)
    return (n_c / (n_c -n_f)) * (V_j / A_c) / M_s


def v_j_avg_sw(V_j, A_w, M_s):
    """Shear stress in shear walls (Eq 4-8) ASCE41-17
    V_j = Story shear at level j computed in accordance with Section 4.4.2.2;
    A_w = Summation of the horizontal cross-sectional area of all
    shear walls in the direction of loading. Openings shall be
    taken into consideration where computing A_w. For masonry walls, the net area shall be used. For wood-framed
    walls, the length shall be used rather than the area; and
    M_s = System modification factor - shall be taken from Table 4-8."""
    return V_j / A_w / M_s


def v_j_avg_br(V_j, s, N_br, L_br, A_br, M_s):
    """Axial stress in diagonal bracing (Eq 4-9) ASCE41-17
    L_br = Average length of the braces (ft);
    N_br = Number of braces in tension and compression if the braces
    are designed for compression, number of diagonal braces
    in tension if the braces are designed for tension only;
    s = Average span length of braced spans (ft);
    A_br = Average area of a diagonal brace (in.2);
    V_j = Maximum story shear at each level (kip); and
    M_s = System modification factor; Ms shall be taken from Table 4-9."""
    return V_j / s / N_br * (L_br * A_br) / () / M_s


def M_gj(V_j, n_c, n_f, h, M_s, perf=None):
    """Precast moment frame connection strength (Eq 4-10) ASCE41-17
    n_c = Total number of columns;
    n_f = Total number of frames in the direction of loading;
    V_j = Story shear at the level directly below the connection
    under consideration;
    h = Typical column story height; and
    M_s = System modification factor from dictionary below"""
    if perf is not None:
        M_s_dict = {'CP':2.5, 'LS':1.5, 'IO':1.0}
        M_s = M_s_dict.get(perf,1.0)
    return  V_j / (n_c - n_f) * h / 2 / M_s


def p_ot(V, h_n, L, n_f, A_col, M_s, perf=None):
    """Column axial stress caused by overturning (Eq 4-11) ASCE41-17
    n_f = Total number of frames in the direction of loading;
    V = Pseudo seismic force;
    h_n = Height (ft) above the base to the roof level;
    L = Total length of the frame (ft);
    M_s = System modification factor from dictionary below
    A_col = Area of the end column of the frame."""
    if perf is not None:
        M_s_dict = {'CP':2.5, 'LS':1.5, 'IO':1.0}
        M_s = M_s_dict.get(perf,1.0)
    return 2 / 3 * (V * h_n) / (L * n_f) / A_col / M_s


def T_c(ksi, S_XS, w_p, A_p, perf=None):
    """Flexible Diaphragm Connection Forces (Eq 4-12) ASCE41-17
    w_p = Unit weight of the wall;
    A_p = Area of wall tributary to the connection;
    ksi from dictionary below
    S_XS = short period spectral acceleration
    """
    if perf is not None:
        ksi_dict = {'CP':1.0, 'LS':1.3, 'IO':1.8}
        ksi = ksi_dict.get(perf,1.8)
    return ksi * S_XS * w_p * A_p


def f_p(f_pe, n_p, A_p):
    """Prestressed element prestress (Eq 4-13) ASCE41-17
    f_pe = Effective force of a prestressed strand;
    n_p = Number of prestressed strands; and
    A_p = Gross area of prestressed concrete elements"""
    return f_pe * n_p / A_p


def f_j_avg_smf(V_j, M_s, n_c, n_f, h, Z, perf=None):
    """Shear stress in coluns and beams of steel moment frames (Eq 4-14) ASCE41-17
    n_c = Total number of frame columns at the level, j, under consideration.
    n_f = Total number of frames in the direction of loading at the level, j, under consideration.
    V_j = Story shear computed in accordance with Section 4.4.2.2.
    h = Story height (in.).
    Z = For columns, the sum of the plastic section moduli of all
    the frame columns at the level under consideration. For
    beams, it is the sum of the plastic section moduli of all the
    frame beams with moment-resisting connections. If a
    beam has moment-resisting connections at both ends, then
    the contribution of that beam to the sum is twice the plastic
    section modulus of that beam (in3).
    M_s = System modification factor; Ms shall be taken from dictionary below
        provided columns and beams satisfy the checklist items for compactness and column axial stress. 
        Otherwise “Noncompliant"""
    if perf is not None:
        M_s_dict = {'CP':9.0, 'LS':6.0, 'IO':2.5}
        M_s = M_s_dict.get(perf,1.0)
    return V_j * (n_c / (n_c -n_f)) * (h / 2 / Z) / M_s


####### Plotting #########################################################
# ASCE41_17_plots, ASCE41_17_tripartite_plot, tripartite_plot_ss

import matplotlib.pyplot as plt
import matplotlib.lines as lines


def ASCE41_17_plots(S_s, S_1, T_code_s = 0.0, soil_classes = ['A', 'B', 'C', 'D', 'E'], S_a_func = S_a_func, site_class='B', **kwargs):
    """
    """
    T_max = max(3.0, 1.2 * T_code_s)
    lim = round(T_max * 100) + 1
    T_range = [0.01 * float(T) for T in range(1,lim)]

    #soil_classes = ['A', 'B', 'C', 'D', 'E']
    S_a_ranges = [[S_a_func(T, S_s, S_1, s_c) 
                   for T in T_range] for s_c in soil_classes]
    S_v_ranges = [[9.81 * 0.5 / pi * T * S_a 
                   for T, S_a in zip(T_range, S_a_range)] for S_a_range in S_a_ranges]
    S_d_ranges = [[0.5 / pi * T * S_v 
                   for T, S_v in zip(T_range, S_v_range)] for S_v_range in S_v_ranges]

    figs, (ax1, ax2, ax3, ax4) = plt.subplots(4,figsize = (10,15))

    [ax1.plot(T_range, S_a_range, label=s_c) 
     for s_c, S_a_range in zip(soil_classes, S_a_ranges)]
    ax1.grid(True)
    ax1.legend(soil_classes)

    [ax2.plot(T_range, S_v_range, label=s_c) 
     for s_c, S_v_range in zip(soil_classes, S_v_ranges)]
    ax2.grid(True)
    ax2.legend(soil_classes)

    [ax3.plot(T_range, S_d_range, label=s_c) 
     for s_c, S_d_range in zip(soil_classes, S_d_ranges)]
    
    ax3.grid(True)
    ax3.legend(soil_classes)

    [ax4.plot(S_d_range, S_a_range, label=s_c) 
     for s_c, S_d_range, S_a_range in zip(soil_classes, S_d_ranges, S_a_ranges)]
    
    ax4.grid(True)
    ax4.legend(soil_classes)
    
    if T_code_s > 0 and site_class in soil_classes:
        S_a_g = S_a_func(T_code_s, S_s, S_1, site_class)
        ax1.plot(T_code_s, S_a_g, marker = 'o')
        S_v_g = 0.5 / pi * T_code_s *  S_a_g * 9.81  # in m/s
        ax2.plot(T_code_s, S_v_g, marker = 'o')
        S_d_g = 0.5 / pi * T_code_s *  S_v_g  # in m
        ax3.plot(T_code_s, S_d_g, marker = 'o')
        ax4.plot(S_d_g, S_a_g, marker = 'o')
    
    return figs


def ASCE41_17_tripartite_plot(S_s, S_1, T_max = 3.0, soil_classes = ['A', 'B', 'C', 'D', 'E'], S_a_func = S_a_func, T_code_s = 0.0, site_class = 'B', **kwargs):
    """Returns a tripartite plot of ASCE41_17 spectra or some other specified function"""
    lim = round(T_max * 100) + 1
    T_range = [0.01 * float(T) for T in range(1,lim) if 0.01 * float(T) >= 0.1]

    S_a_ranges = [[S_a_func(T, S_s, S_1, s_c) 
                   for T in T_range] 
                   for s_c in soil_classes]
    
    plot_pts = []
    if T_code_s > 0 and site_class in soil_classes:
        S_a_g = S_a_func(T_code_s, S_s, S_1, site_class)
        S_v_g = 0.5 / pi * T_code_s *  S_a_g * 9.81  # in m/s
        plot_pts.append((T_code_s, S_v_g))
    
    return tripartite_plot_ss(
        T_range, 
        S_a_ranges, 
        soil_classes, 
        plot_pts=plot_pts
        )


def tripartite_plot_ss(T_range, S_a_ranges, titles, plot_pts=[], **kwargs):
    """Plots tripartite data
        majorlines is a list of the minor log lines to plot (default [1,2,3,4,5,6,7,8,9])
    
    Args:
        T_range (list[float]: time values - it is assumed that the range is the same for all plots
        S_a_ranges: spectral acceleration (assumed to be in units of g)
        titles (list[str]): one title for each S_a_range
        plot_pts (list[tuple,list]): (T, S_v) points for plotting 
        kwargs: parameters for passing through to? Nowhere at present 
            (can specify figsize and majorlines)
    """
    # Check that data dimensions match
    err_msgs = []
    if not (len(S_a_ranges) == len(titles)):
        err_msgs.append(f'Number of curves does not match: len(S_a_ranges) = {len(S_a_ranges)}, len(titles) = {len(titles)}')
    S_a_ranges_len = [len(S_a_range) for S_a_range in S_a_ranges]
    if not (all(S_a_range_len == len(T_range) 
                for S_a_range_len in S_a_ranges_len)):
        err_msgs.append(f'Number of curves does not match: len(T_range) = {len(T_range)}, len(S_a_ranges) = {S_a_ranges_len}')
    if len(err_msgs) > 0:
        err_msg = '\n'.join(err_msgs)
        raise ValueError(err_msg)
    
    # calculate the spectral pseudo velocities
    S_v_ranges = [[9.81 * 0.5 / pi * T * S_a 
                   for T, S_a in zip(T_range, S_a_range)] 
                   for S_a_range in S_a_ranges]
    
    figsize = kwargs.get('figsize', (10,8))
    fig0, ax0 = plt.subplots(figsize = figsize)
    [ax0.loglog(T_range, S_v_range, label=s_c) 
     for s_c, S_v_range in zip(titles, S_v_ranges)]
    # ax0.loglog(T_code_s, S_v_g, marker = 'o') # marking the results for the building period

    ax0.grid(True)
    #fig = plt.gcf(); fig.set_size_inches(16,8)
    ax0.legend(titles)

    ##### Tripartite Lines ###########
    x0 = fig0.axes[0].get_xlim() #; print(x0)
    y0 = fig0.axes[0].get_ylim() #; print(y0)
    
    T_min, T_max = x0  
    log_T_min, log_T_max = log10(T_min), log10(T_max)
    #print('T_min, T_max', T_min, T_max)
    
    v_min, v_max = y0  
    log_v_min, log_v_max = log10(v_min), log10(v_max)
    #print('v_min, v_max', v_min, v_max)
    
    a_min, a_max = [2.0 * pi / T * v 
                    for v, T in ((v_min, T_max),(v_max, T_min))]
    log_a_min, log_a_max = log10(a_min), log10(a_max)
    #print('a_min, a_max', a_min, a_max)
    
    d_min, d_max = [0.5 / pi * T * v 
                    for v, T in ((v_min, T_min),(v_max, T_max))]
    log_d_min, log_d_max = log10(d_min), log10(d_max)
    #print('d_min, d_max', d_min, d_max)
    
    # Calculate a-values and d-values on the log-log line between a_min & a_max and d_min & d_max respectively
    major_lines = kwargs.get('major_lines', [1, 2, 3, 4, 5, 6, 7, 8, 9])
    
    a_fac, a_unit = (1.0, 'm/s^2')
    #a_fac, a_unit = (9.81, 'g') # does not work!
    
    d_fac, d_unit = (1.0, 'm')
    
    #===========
    # generate the list of acceleration lines provided they are within the limits
    a_list = [(m * 10**float(n)) 
              for n in range(int(log10(a_min))-2, int(log10(a_max))) 
              for m in major_lines if (m * 10**float(n) * a_fac > a_min) and (m * 10**float(n) * a_fac < a_max)]
    a_prec = [max(0,int(1+log10(0.99/a))) for a in a_list]
    a_str = ['{:.{prec}f} {}'.format(a, a_unit, prec=p) 
             for a, p in zip(a_list, a_prec)]
    a_values = [a * a_fac for a in a_list]
    #[print(a,b,s) for a,b,s in zip(a_list, a_values, a_str)]
    
    # generate the list of displacement lines provided they are within the limits
    d_list = [(m * 10**float(n)) 
              for n in range(int(log10(d_min)), int(log10(d_max))+1) for m in major_lines if (m * 10**float(n) * d_fac > d_min) and (m * 10**float(n) / d_fac < d_max)]
    d_prec = [max(0,int(1+log10(0.99/d))) 
              for d in d_list]
    d_str = ['{:.{prec}f} {}'.format(d, d_unit, prec=p) 
             for d, p in zip(d_list, d_prec)]
    d_values = [d * d_fac for d in d_list]
    #print('d_values (min & max)', d_values[0], d_values[-1])
    # 
    log_a_values = [log10(a) for a in a_values]
    log_d_values = [log10(d) for d in d_values]

    #change the following lines
    #a_ls = ['-' if abs(a + log10(a_fac) - int(a + log10(a_fac))) < 0.001 else ':' for a in log_a_values]
    a_ls = ['-' if abs(log10(a) - int(log10(a))) < 0.001 else ':' for a in a_list]
    #d_ls = ['-' if abs(d - int(d)) < 0.001 else ':' for d in log_d_values]
    d_ls = ['-' if abs(log10(d) - int(log10(d))) < 0.001 else ':' for d in d_list]

    ta_values = [(a - log_a_min)/(log_a_max - log_a_min) for a in log_a_values]
    #print('ta_min, ta_max', ta_values[0], ta_values[-1])
    td_values = [(d - log_d_min)/(log_d_max - log_d_min) for d in log_d_values]
    #print('td_min, td_max', td_values[0], td_values[-1])

    va_values = [log_v_min + t * (log_v_max - log_v_min) for t in ta_values]
    #print('va_min, va_max', 10**va_values[0], 10**va_values[-1])
    vd_values = [log_v_min + t * (log_v_max - log_v_min) for t in td_values]
    #print('vd_min, vd_max', 10**vd_values[0], 10**vd_values[-1])

    Ta_values = [log_T_max + t * (log_T_min - log_T_max) for t in ta_values]
    #print('Ta_min, Ta_max', 10**Ta_values[0], 10**Ta_values[-1])
    Td_values = [log_T_min + t * (log_T_max - log_T_min) for t in td_values]
    #print('Td_min, Td_max', 10**Td_values[0], 10**Td_values[-1])

    #s1
    s_a2 = [min(log_v_max - v, log_T_max - T) for v, T in zip(va_values, Ta_values)]
    s_a1 = [max(log_v_min - v, log_T_min - T) for v, T in zip(va_values, Ta_values)]
    #print('s_a1[0], s_a2[0]', s_a1[0], s_a2[0])

    # identifying the index of the longest acceleration line
    (m_a, i_a) = max((v2-v1,i) for i, (v1, v2, als) in enumerate(zip(s_a1, s_a2, a_ls)) if als == '-')
    #print('m_a, i_a', m_a, i_a)

    #s1
    s_d2 = [max(v - log_v_max, log_T_min - T) for v, T in zip(vd_values, Td_values)]
    s_d1 = [min(v - log_v_min, log_T_max - T) for v, T in zip(vd_values, Td_values)]
    #print('s_d1[0], s_d2[0]', s_d1[0], s_d2[0])

    # identifying the index of the longest displacement line
    (m_d, i_d) = max((v1-v2,i) for i, (v1, v2, dls) in enumerate(zip(s_d1, s_d2, d_ls)) if dls == '-')
    #print('m_d, i_d', m_d, i_d)

    # calculate the coordinates of the ends of the acceleration lines
    a_coords = [((10**(T+s1), 10**(T+s2)), (10**(v+s1), 10**(v+s2))) for v, T, s1, s2 in zip(va_values, Ta_values, s_a1, s_a2)]
    a_lines = [lines.Line2D(P1, P2, transform=ax0.transData, axes=ax0, color='grey', ls = ls, lw = 0.5, zorder=0) for (P1, P2), ls in zip(a_coords, a_ls)]

    # calculate the coordinates of the ends of the displacement lines
    d_coords = [((10**(T+s1), 10**(T+s2)), (10**(v-s1), 10**(v-s2))) for v, T, s1, s2 in zip(vd_values, Td_values, s_d1, s_d2)]
    d_lines = [lines.Line2D(P1, P2, transform=ax0.transData, axes=ax0, color='grey', ls = ls, lw = 0.5, zorder=0) for (P1, P2), ls in zip(d_coords, d_ls)]

    #print('a_coords', a_coords[i_a])
    #print('d_coords', d_coords[i_d])

    ax_lines = [lines.Line2D(P1, P2, transform=ax0.transData, axes=ax0, color='magenta', ls = '-', lw = 1.0, zorder=0) for (P1, P2) in [a_coords[i_a], d_coords[i_d]]]

    ax0.lines.extend(a_lines + d_lines + ax_lines) #, l3, l4])

    #ax0.text(0.5*(T_min + T_max), 0.5*(v_min + v_max), "text on plot $m/s^2$")
    [ax0.text(10**T, 10**v, '$'+txt+'$', rotation=45, zorder=0.1) for T, v, txt in zip(Ta_values, va_values, a_str)]
    [ax0.text(10**T, 10**v, '$'+txt+'$', rotation=-45, zorder=0.1) for T, v, txt in zip(Td_values, vd_values, d_str)]
    
    # plot any additionally provided points
    if len(plot_pts) > 0:
        [ax0.plot(T, v, marker = 'o') for T,v in plot_pts]
    
    return fig0



"""A library for carrying out calculations to ASCE41-13

ASCE 41-13 which merges ASCE 31-03 (evaluation) with 
ASCE 41-06 (rehabilitation and retrofit)

* ASCE/SEI 41-06, Seismic Rehabilitation of Existing Buildings
* ASCE/SEI 31-03, Seismic Evaluation of Existing Buildings.

The predecessor to ASCE 31 was FEMA 310, and the 
predecessor to ASCE 41 was FEMA 356 (in turn based
on FEMA 273).

Note that all references to ASCE 7 are related generically 
to ASCE 7-10

The full list of references may be found in section C1.3.

* ASCE 41-13
* ASCE 41-06

TODO:
* Everything...

"""


from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
from itertools import accumulate

from pySeismic.utilities import printt, interp, list_interp, \
    storey_force_coeffs, storey_force_distribution, storey_forces, \
    storey_shear_moment, shear_coords
# from pySeismic.importing import k_gen # Eq 4-2a
from pySeismic.codes.ASCE7_10 import site_coeffs as ASCE7_10_site_coeffs


# -----------------------------------------------------------------
# ----- Tier One Data and Functions -------------------------------
# -----------------------------------------------------------------


# Table 2-5 Level of Seismicity Definitions
def seismicity_level_get(S_s, S_1, site_class):
    """Returns level of seismicity, Table 2-5, ASCE41-13
    
    S_s & S_1 are the short period and 1-sec period spectral
    acceleration associated with BSE-2N (section 2.4.1.1)

    Args:
        S_s: The short period site spectral acceleration
            :units: g
            :latex: \\S_S
        S_1: The 1-second site spectral acceleration
            :units: g
            :latex: \\S_1
        site_class: The site classification (A-E)
            :latex: \\class
    
    Returns:
        Seismicity: The level of seismicity per section 2.5
            :latex: \\Seismicity
    
    .. codeauthor:: Andrew Mole <andrew.mole@arup.com>
    
    .. levelofreview:: Work in progress
    """
    S_DS_line = (0.0, 0.167, 0.33, 0.5)
    S_D1_line = (0.0, 0.067, 0.133, 0.2)
    F_a, F_v = ASCE7_10_site_coeffs(site_class, S_s, S_1)
    S_DS = 2/3 * F_a * S_s  # Eq 2-12
    S_D1 = 2/3 * F_v * S_1  # Eq 2-13
    seismicity_levels = ('Very Low', 'Low', 'Moderate', 'High')
    i_1 = [i for i, s_ds in enumerate(S_DS_line) if s_ds <= S_DS][-1]
    i_2 = [i for i, s_d1 in enumerate(S_D1_line) if s_d1 <= S_D1][-1]
    return seismicity_levels[max(i_1, i_2)]



# -----------------------------------------------------------------
# ----- Ch 3 Evaluation and Retrofit Requirements -----------------
# -----------------------------------------------------------------


# Table 3-4  ASCE41-17 Limitations on Use of the Tier 1 & Tier 2 Procedures (number of storeys)
# No Limit (NL) = -1, Not Permitted (NP) = 0
struct_types_full_dict = {
'W1': {'Description':'Light Wood Frames', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):4, 
       ('Moderate', 'S-3'):4, ('Moderate', 'S-1'):4, ('High', 'S-3'):4, ('High', 'S-1'):4},
'W1a': {'Description':'Multistory Multi-Unit Residential Wood Frames', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-3'):6, ('Moderate', 'S-1'):6, ('High', 'S-3'):6, ('High', 'S-1'):4},
'W2': {'Description':'Commercial and Industrial Wood Frames', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-3'):6, ('Moderate', 'S-1'):6, ('High', 'S-3'):6, ('High', 'S-1'):4},

'S1': {'Description':'Steel MRF Stiff diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'S1a': {'Description':'Steel MRF Flexible diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'S2': {'Description':'Steel Braced Frames Stiff diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-3'):8, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'S2a': {'Description':'Steel Braced Frames Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-3'):8, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'S3': {'Description':'Steel Light Frames', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):1, ('Low', 'S-3'):1, ('Low', 'S-1'):1, 
       ('Moderate', 'S-3'):1, ('Moderate', 'S-1'):1, ('High', 'S-3'):1, ('High', 'S-1'):1},
'S4': {'Description':'Dual Systems with Backup Steel MRF', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'S5': {'Description':'Steel Frame Infill URM Stiff Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):4},
'S5a': {'Description':'Steel Frame Infill URM Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):4},
'S6': {'Description':'Steel Plate Shear Wall', 
       ('Very Low', 'S-3'):0, ('Very Low', 'S-1'):0, ('Low', 'S-3'):0, ('Low', 'S-1'):0, 
       ('Moderate', 'S-3'):0, ('Moderate', 'S-1'):0, ('High', 'S-3'):0, ('High', 'S-1'):0},

'C1': {'Description':'RC MRF', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'C2': {'Description':'RC Shear Wall Stiff Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'C2a': {'Description':'RC Shear Wall Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'C3': {'Description':'RC Frame URM Infill Rigid diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):4},
'C3a': {'Description':'RC Frame URM Infill Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):12, 
       ('Moderate', 'S-3'):12, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):4},
'PC1': {'Description':'Precast/Tilt-Up RC SW Stiff diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):3, ('Low', 'S-1'):2, 
       ('Moderate', 'S-3'):2, ('Moderate', 'S-1'):2, ('High', 'S-3'):2, ('High', 'S-1'):2},
'PC1a': {'Description':'Precast/Tilt-Up RC SW Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):3, ('Low', 'S-1'):2, 
       ('Moderate', 'S-3'):2, ('Moderate', 'S-1'):2, ('High', 'S-3'):2, ('High', 'S-1'):2},
'PC2': {'Description':'Precast RC Frame with Shear Walls', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-3'):6, ('Moderate', 'S-1'):0, ('High', 'S-3'):4, ('High', 'S-1'):0},
'PC2a': {'Description':'Precast RC Frame without Shear Walls', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):6, 
       ('Moderate', 'S-3'):6, ('Moderate', 'S-1'):0, ('High', 'S-3'):4, ('High', 'S-1'):0},

'RM1': {'Description':'Reinforced Masonry Bearing Walls Rigid diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-3'):8, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},
'RM2': {'Description':'Reinforced Masonry Bearing Walls Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):-1, ('Low', 'S-1'):8, 
       ('Moderate', 'S-3'):8, ('Moderate', 'S-1'):8, ('High', 'S-3'):8, ('High', 'S-1'):6},

'URM': {'Description':'Unreinforced Masonry Bearing Walls Stiff diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):6, ('Low', 'S-1'):4, 
       ('Moderate', 'S-3'):6, ('Moderate', 'S-1'):0, ('High', 'S-3'):4, ('High', 'S-1'):0},
'URMa': {'Description':'Unreinforced Masonry Bearing Walls Flexible Diaphragm', 
       ('Very Low', 'S-3'):-1, ('Very Low', 'S-1'):-1, ('Low', 'S-3'):6, ('Low', 'S-1'):4, 
       ('Moderate', 'S-3'):6, ('Moderate', 'S-1'):0, ('High', 'S-3'):4, ('High', 'S-1'):0},

'SI_PD': {'Description':'Seismic Isolation or Passive Dissipation', 
       ('Very Low', 'S-3'):0, ('Very Low', 'S-1'):0, ('Low', 'S-3'):0, ('Low', 'S-1'):0, 
       ('Moderate', 'S-3'):0, ('Moderate', 'S-1'):0, ('High', 'S-3'):0, ('High', 'S-1'):0},
       }


# -----------------------------------------------------------------
# ----- Tier One Data and Functions -------------------------------
# -----------------------------------------------------------------
# Simplified Rehabilitation

# M_s_func, v_j_avg_mf, D_r, v_j_avg_sw



def S_a_Tier_1_func(S_s: float, S_1: float, T: float, site_class: str) -> float:
    """Returns spectral acceleration for Tier 1 Assessment

    4.4.2.3 Spectral Acceleration ASCE41-13
    For BSE-2N: use S_s & S_1 values from MCE_R contour maps in Ch 22 ASCE7-10
    For BSE-1N: use two-thirds of the values for BSE-2N
    For BSE-2E: use S_s & S_1 values from 5%/50-year contour maps
    For BSE-1E: use S_s & S_1 values from 20%/50-year contour maps

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
    F_a, F_v = ASCE7_10_site_coeffs(site_class, S_s, S_1)
    S_XS = F_a * S_s  # Eq 2-1
    S_X1 = F_v * S_1  # Eq 2-2
    return min(S_XS, S_X1 / T)  # Eq 4-3


# ASCE41-13 Table 4-2 Default compressive strengths (f'c) in ksi
table_4_2_dict_13 = {
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
table_4_3_dict_13 = {
             (1911, 1959): {33:True, 40:True, 50:True, 60:False, 65:True, 70:False, 75:False},
             (1959, 1966): {33:True, 40:True, 50:True, 60:True, 65:True, 70:True, 75:True},
             (1966, 1987): {33:False, 40:True, 50:True, 60:True, 65:True, 70:True, 75:False},
             (1987, 3000): {33:False, 40:True, 50:True, 60:True, 65:True, 70:True, 75:True},
             }


# Table 4-4. Default Yield Strengths (Fy) of Archaic Materials (ksi)
table_4_4_dict_13 = {
    ('Pre-1900', 'Cast iron'): 18,
    ('Pre-1900', 'Steel'): 24,
    }



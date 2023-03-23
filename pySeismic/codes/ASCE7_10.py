"""A library for carrying out calculations to ASCE7-10
"""

from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
# from itertools import accumulate

from pySeismic.utilities import list_interp

#### Seismic Hazard ####

# Functions
# T_code_get, site_coeffs, site_accels_M, seismicity_level_get, \
# code_period, S_a_get, seismic_design_category_get, k_gen

# not used
# def C_s(S_D1, T, R, I_e, T_L):
#     res = (S_D1*T_L/T**2 if T > T_L else (S_D1*/T if T > (S_D1/S_DS) else S_DS)) * I_e / R
#     pass


def C_u_interp_func(S_D1): 
    """Coefficient for Upper Limit on Calculated Period (Table 12.8-1) ASCE7-10"""
    return list_interp(
                        S_D1, 
                        [0.1, 0.15, 0.2, 0.3, 0.4], 
                        [1.7, 1.6, 1.5, 1.4, 1.4]
                       )


# Empirical period factors (returns a tuple of imperial and metric)
C_t_dict = {'SMRF':(0.028,0.0724), 'CMRF':(0.016,0.0466), 'SEBF': (0.03,0.0731), 
            'SBRBF': (0.03,0.0731), 'other': (0.02,0.0488)}
C_t_imp_dict = {k: v[0] for k, v in C_t_dict.items()}
C_t_metric_dict = {k: v[1] for k, v in C_t_dict.items()}
# C_t_imp_dict = {'SMRF':0.028, 'CMRF':0.016, 'SEBF': 0.03, 'SBRBF': 0.03, 'other': 0.02}
# C_t_metric_dict = {'SMRF':0.0724, 'CMRF':0.0466, 'SEBF': 0.0731, 'SBRBF': 0.0731, 'other': 0.0488}
x_dict = {'SMRF':0.8, 'CMRF':0.9, 'SEBF': 0.75, 'SBRBF': 0.75, 'other': 0.75}


def x_get(lateral_system):
    """Building period parameter, returns default of 0.75
    Table 12.8-2 (ASCE7-10)"""
    return x_dict.get(lateral_system, 0.75)

def C_t_get(lateral_system):
    """Building period parameter, returns default of 0.020
    Table 12.8-2 (ASCE7-10)"""
    return C_t_dict.get(lateral_system, (0.020, 0.0488))

unit_l_dict = {'mm':0.001, 'cm': 0.01, 'm': 1.0, 
               'in':0.0254, 'ft':0.3048, 'yd':0.9144}

m_unit_l_dict = {'mm':0.001, 'cm': 0.01, 'm': 1.0}

ft_unit_l_dict = {'in':1/12, 'ft':1.0, 'yd':3.0}


def unit_dict_get(unit_dict, lateral_system):
    """a units lookup function that returns default of None"""
    return unit_dict.get(lateral_system, None)  


def T_code_get(lateral_system, height, units = 'm'):
    """Returns the code period in seconds for 
    a given structure type and height Eqn 12.8-7, ASCE7-10"""
    usys, unit_dict = (0, ft_unit_l_dict) if units in ('in', 'ft', 'yd') \
        else (1, m_unit_l_dict)
    opt = 0 if units in ('m', 'ft') else 1
    C_t = C_t_get(lateral_system)[usys]
    x = x_get(lateral_system)
    u_fac = unit_dict_get(unit_dict, units)
    return C_t * height ** x if opt == 0 \
        else C_t * (u_fac * height) ** x


# def B_1_func(x):
#     """for x as damping ratio (not percentage), ASCE7-10"""
#     return 4.0 / (5.6 - log(100 * x))


# F_a and F_v taken from ASCE7-10, tables 11.4-1 & 11.4-2
F_a_headers = (0.25,0.5,0.75,1.0,1.25,1.5)
F_v_headers = (0.1,0.2,0.3,0.4,0.5,0.6)
F_a_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
           'B': (1.0,1.0,1.0,1.0,1.0,1.0),
           'C': (1.2,1.2,1.1,1.0,1.0,1.0),
           'D': (1.6,1.4,1.2,1.1,1.0,1.0),
           'E': (2.5,1.7,1.2,0.9,0.9,0.9),
           'F': (None,None,None,None,None,None)}
F_v_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
           'B': (1.0,1.0,1.0,1.0,1.0,1.0),
           'C': (1.7,1.6,1.5,1.4,1.3,1.3),
           'D': (2.4,2.0,1.8,1.6,1.5,1.5),
           'E': (3.5,3.2,2.8,2.4,2.4,2.4),
           'F': (None,None,None,None,None,None)}


# Importance Factors - Table 1.5-2 ASCE7-10
importance_factors_headers = ('I_s', 'I_i', 'I_w', 'I_e')
importance_factors_dict = {'I': (0.80, 0.80, 1.00, 1.00), 
                            'II': (1.00, 1.00, 1.00, 1.00), 
                            'III': (1.00, 1.10, 1.00, 1.25), 
                            'IV': (1.20, 1.25, 1.00, 1.50)}


def site_coeffs(site_class, S_s, S_1):
    """Returns F_a, F_v when provided with site class (A-F) and S_s and S_1 values"""
    F_a_list = F_a_dict[site_class]
    F_v_list = F_v_dict[site_class]
    return list_interp(S_s, F_a_headers, F_a_list), list_interp(S_1, F_v_headers, F_v_list)


def site_accels_M(site_class, S_s, S_1):
    """Returns S_MS, S_M1 when provided with site class (A-F) and S_s and S_1 values
    >>> site_accels_M(site_class, S_s, S_1)
    (0.468, 0.33)
    """
    F_a_list = F_a_dict[site_class]
    F_v_list = F_v_dict[site_class]
    F_a = list_interp(S_s, F_a_headers, F_a_list)
    F_v = list_interp(S_1, F_v_headers, F_v_list)
    return F_a * S_s, F_v * S_1


def seismicity_level_get():
    """Not implemented"""
    pass


def code_period(h_n, s_type):
    """Returns period in seconds when given height (h_n) in metres and ASCE41-17 structural type (e.g. 'C1')"""
    return C_t_get(s_type)[1] * (h_n / 0.3048) ** x_get(s_type)


def S_DS_get(F_a, S_s):
    """Incomplete - not needed"""
    S_DS = 2 / 3 * F_a * S_s
    return S_DS


def S_a_get(T, S_s, S_1, site_class, T_L = 12.0):
    """Returns the spectral acceleration when provided with period (T), S_s, S_1 and site class
    using values for F_a & F_v from tables 11.4-1 & 11.4-2"""
    F_a, F_v = site_coeffs(site_class, S_s, S_1)
    S_DS = 2 / 3 * F_a * S_s
    S_D1 = 2 / 3 * F_v * S_1
    
    T_S = S_D1 / S_DS
    T_0 = 0.2 * T_S
    return min(S_DS, S_DS*(0.4 + 0.6*T/T_0), S_D1 / T) if T < T_L else T_L * S_D1 / (T**2.0)


def seismic_design_category_get(S_s, S_1, site_class, risk_category):
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
    F_a, F_v = site_coeffs(site_class, S_s, S_1)
    S_DS = 2 / 3 * F_a * S_s
    S_D1 = 2 / 3 * F_v * S_1
    print('F_a, F_v:', F_a, F_v)
    print('S_DS, S_D1:', S_DS, S_D1)
    
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
    
    print(seismic_design_category_S, seismic_design_category_1)
    seismic_design_category = max(seismic_design_category_S, seismic_design_category_1)
    
    if risk_category in ['I', 'II', 'III'] and S_1 >= 0.75:
        seismic_design_category = 'E'
    
    if risk_category == 'IV' and S_1 >= 0.75:
        seismic_design_category = 'F'
    
    return seismic_design_category


#### Seismic Hazard ####

# Functions
# k_gen

def k_gen(T):
    """Returns the coefficient for the vertical distribution of lateral force (Eq 4-2a)"""
    return max(1.0, min(2.0, 0.5 * T + 0.75))




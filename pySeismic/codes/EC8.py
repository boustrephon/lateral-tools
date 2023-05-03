"""A module for EC8 calculations based on default values

Eurocode 8: Design of structures for earthquake resistance — 
    Part 1: General rules, seismic actions and rules for buildings

Default data is used, but the options exist to replace these with values
    from the National Application Documents (NAD).

There are two types of spectra: Type 1 and Type 2. If the 
earthquakes that contribute most to the seismic hazard defined 
for the site for the purpose of probabilistic hazard assessment 
have a surface-wave magnitude, Ms, not greater than 5,5, it 
is recommended that the Type 2 spectrum is adopted. 

5.2.2.1(1)P Concrete buildings shall be classified into one of the following structural 
types (see 5.1.2) according to their behaviour under horizontal seismic actions: 
a) frame system; 
b) dual system (frame or wall equivalent); 
c) ductile wall system (coupled or uncoupled); 
d) system of large lightly reinforced walls; 
e) inverted pendulum system; 
f) torsionally flexible system


ductile wall 
    wall fixed at its base so that the relative rotation of this base with respect to the rest of 
    the structural system is prevented, and that is designed and detailed to dissipate energy 
    in a flexural plastic hinge zone free of openings or large perforations, just above its base 

large lightly reinforced wall 
    wall with large cross-sectional dimensions, that is, a horizontal dimension lw at least 
    equal to 4,0 m or two-thirds of the height hw of the wall, whichever is less, which is 
    expected to develop limited cracking and inelastic behaviour under the seismic design 
    situation 

coupled wall 
    structural element composed of two or more single walls, connected in a regular pattern 
    by adequately ductile beams ("coupling beams"), able to reduce by at least 25% the sum 
    of the base bending moments of the individual walls if working separately 

wall system 
    structural system in which both vertical and lateral loads are mainly resisted by vertical 
    structural walls, either coupled or uncoupled, whose shear resistance at the building 
    base exceeds 65% of the total shear resistance of the whole structural system

frame system 
    structural system in which both the vertical and lateral loads are mainly resisted by 
    spatial frames whose shear resistance at the building base exceeds 65% of the total shear 
    resistance of the whole structural system 

dual system 
    structural system in which support for the vertical loads is mainly provided by a spatial 
    frame and resistance to lateral loads is contributed to in part by the frame system and in 
    part by structural walls, coupled or uncoupled

frame-equivalent dual system 
    dual system in which the shear resistance of the frame system at the building base is 
    greater than 50% of the total shear resistance of the whole structural system 

wall-equivalent dual system 
    dual system in which the shear resistance of the walls at the building base is higher than 
    50% of the total seismic resistance of the whole structural system 

torsionally flexible system 
    dual or wall system not having a minimum torsional rigidity (see 5.2.2.1(4)P and (6))

inverted pendulum system 
    system in which 50% or more of the mass is in the upper third of the height of the 
    structure, or in which the dissipation of energy takes place mainly at the base of a single 
    building element

Behaviour factors (q) are defined for different structural types:
- concrete: section 5.2.2 & 5.3.3
- precast concrete: section 5.11.1.4
- steel: section 6.3.2
- composite steel: 7.3.2
- timber: 8.3
- masonry: 9.3
- base isolation: 10.7

"""


concrete_structural_system_dict = {
    'CF1S': {'Description': 'one-storey frame system', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.5, 'alphaM': True, 'alphaH': True, 'alpha_ratio': 1.1, 'k_w': 1.0, }, 
    'CF1B': {'Description': 'one-bay frame system', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.5, 'alphaM': True, 'alphaH': True, 'alpha_ratio': 1.2, 'k_w': 1.0, }, 
    'CFMB': {'Description': 'multi-bay frame system', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.5, 'alphaM': True, 'alphaH': True, 'alpha_ratio': 1.3, 'k_w': 1.0, }, 
    'CFED': {'Description': 'dual system (frame equivalent)', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.5, 'alphaM': True, 'alphaH': True, 'alpha_ratio': 1.3, 'k_w': 1.0, }, 
    'CFTF': {'Description': 'torsionally flexible frame system', 
         'DCL': 1.5, 'DCM': 2.0, 'DCH': 3.0, 'alphaM': False, 'alphaH': False, 'alpha_ratio': 1.0, 'k_w': 1.0, },
    'CWSC': {'Description': 'ductile wall system (coupled)', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.5, 'alphaM': True, 'alphaH': True, 'alpha_ratio': 1.2, 'k_w': None, }, 
    'CW2U': {'Description': 'two-wall ductile wall system (uncoupled)', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.0, 'alphaM': False, 'alphaH': True, 'alpha_ratio': 1.0, 'k_w': None, }, 
    'CWMU': {'Description': 'multi-wall ductile wall system (uncoupled)', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.0, 'alphaM': False, 'alphaH': True, 'alpha_ratio': 1.1, 'k_w': None, }, 
    'CWLR': {'Description': 'system of large lightly reinforced walls', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 3.0, 'alphaM': False, 'alphaH': False, 'alpha_ratio': 1.2, 'k_w': None, }, 
    'CWED': {'Description': 'dual system (wall equivalent)', 
         'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.5, 'alphaM': True, 'alphaH': True, 'alpha_ratio': 1.2, 'k_w': None, }, 
    'CWTF': {'Description': 'torsionally flexible wall system', 
         'DCL': 1.5, 'DCM': 2.0, 'DCH': 3.0, 'alphaM': False, 'alphaH': False, 'alpha_ratio': 1.0, 'k_w': None, },
    'CIPS': {'Description': 'inverted pendulum system', 
         'DCL': 1.5, 'DCM': 2.0, 'DCH': 4.5, 'alphaM': False, 'alphaH': False, 'alpha_ratio': 1.0, 'k_w': None, }, 
    }


table_5_1_q0_dict = {
    'Dual': {'DCM': 3.0, 'DCH': 4.5},
    'Wall': {'DCM': 3.0, 'DCH': 4.0},
    'Flexible': {'DCM': 2.0, 'DCH': 3.0},
    'Pendulum': {'DCM': 1.5, 'DCH': 2.0},
}


# Maximum value is 1.5
alpha_ratio_dict = {
    'Frame Single': 1.1,
    'Frame Multi 1-bay': 1.2,
    'Frame Multi n-bay': 1.3,
    'Frame Dual': 1.3,
    'Two Wall': 1.0,
    'Wall': 1.1,
    'Wall Dual': 1.2,
}


def get_alpha_ratio(
        structural_system: str
    ) -> float:
    """Alpha ratio - alpha_u / alpha_1

    For regular structural systems, an approximate value for 
    alpha_u / alpha_1 can be based on structural systems (section 5.2.2.2(5))
        a) Frames or frame-equivalent dual systems. 
        - 'Frame Single' One-storey buildings: alpha_ratio = 1.1 
        - 'Frame Multi 1-bay' multistorey, one-bay frames: alpha_ratio = 1.2 
        - 'Frame Multi n-bay' multistorey, multi-bay frames: alpha_ratio = 1.3
        - 'Frame Dual' multistorey, multi-bay frame-equivalent dual structures: alpha_ratio = 1.3
        b) Wall- or wall-equivalent dual systems. 
        - 'Two Wall' wall systems with only two uncoupled walls per horizontal direction: alpha_ratio = 1.0
        - 'Wall' other uncoupled wall systems: alpha_ratio = 1.1 
        - 'Wall Dual' wall-equivalent dual, or coupled wall systems: alpha_ratio = 1.2

    Alternatively, alpha_u and alpha_1 may be calculated directly (5.2.2.2(4))
    alpha_1 is the value by which the horizontal seismic design action is multiplied in order 
        to first reach the flexural resistance in any member in the structure, while all 
        other design actions remain constant; 
    alpha_u is the value by which the horizontal seismic design action is multiplied, in order 
        to form plastic hinges in a number of sections sufficient for the development of 
        overall structural instability, while all other design actions remain constant. The 
        factor alpha_u may be obtained from a nonlinear static (pushover) global analysis.
        
    Args:
        structural_system: frame or wall options (defined above)
            frame options - 'Frame Single', 'Frame Multi 1-bay', 'Frame Multi n-bay', 'Frame Dual'
            wall options - 'Two Wall', 'Wall', 'Wall Dual'

    Returns:
        alpha_ratio: ratio of alpha_u over alpha_1
    """
    alpha_ratio = alpha_ratio_dict.get(structural_system, None)
    return alpha_ratio


def get_q_ss(
        structural_system: str, 
        ductility_class: str,
        alpha_ratio: float = 1.0, 
        alpha_0: float = 2.0, 
        regular: bool = False, 
    ) -> float:
    """EC8 Behaviour factor q0

    The basic behaviour factor for regular systems is provided
    in table 5.1.
    Note: q0 is reduced by 80% if the system is not regular
    
    Args:
        structural_system:
            'Dual' - Frame system, dual system, coupled wall system
            'Wall' - Uncoupled wall system
            'Flexible' - Torsionally flexible system
            'Pendulum' - Inverted pendulum system 
        ductility_class: low, medium or high (DCL, DCM or DCH)
        alpha_ratio: this value can be calculated directly 
            (see 5.2.2.2(4)) or taken from 5.2.2.2(5)
        alpha_0: prevailing aspect ratio of the walls of
            the structural system - sum(h_wi) / sum(l_wi)
            h_wi is height of wall i
            l_wi is length of wall i
        regular: whether structural system is regular (default is False)
            (see 4.2.3.1(7) and Table 4.1)
    Returns:
        q: behaviour factor
    """
    ar_limit = min(1.5, alpha_ratio)
    table_5_1_q0_dict = {
        'Dual': {'DCL': 1.5, 'DCM': 3.0 * ar_limit, 'DCH': 4.5 * ar_limit},
        'Wall': {'DCL': 1.5, 'DCM': 3.0, 'DCH': 4.0 * ar_limit},
        'Flexible': {'DCL': 1.5, 'DCM': 2.0, 'DCH': 3.0},
        'Pendulum': {'DCL': 1.5, 'DCM': 1.5, 'DCH': 2.0},
        }
    q0 = table_5_1_q0_dict.get(structural_system).get(ductility_class, 1.5) * 1.0 if regular else 0.8
    k_w = max(0.5, min(1 + alpha_0)/3, 1.0)
    q = max(1.5, q0 * k_w)
    return q


def get_q(
        structural_system: str, 
        ductility_class: str,
        alpha_0: float = 0.5, 
        alpha_ratio: float = None, 
        plan_regular: bool = False, 
        elev_regular: bool = False, 
    ) -> float:
    """EC8 Behaviour factor q0

    The basic behaviour factor for regular systems is provided
    in table 5.1.
    Note: q0 is reduced by 80% if the system is not regular
    
    Args:
        structural_system:
            CF1S - one-storey frame system
            CF1B - one-bay frame system
            CFMB - multi-bay frame system
            CFED - dual system (frame equivalent)
            CFTF - torsionally flexible frame system
            CWSC - ductile wall system (coupled)
            CW2U - two-wall ductile wall system (uncoupled)
            CWMU - multi-wall ductile wall system (uncoupled)
            CWLR - system of large lightly reinforced walls
            CWED - dual system (wall equivalent)
            CWTF - torsionally flexible wall system
            CIPS - inverted pendulum system
        ductility_class: low, medium or high (DCL, DCM or DCH)
        alpha_0: prevailing aspect ratio of the walls of
            the structural system - sum(h_wi) / sum(l_wi)
            h_wi is height of wall i
            l_wi is length of wall i
        alpha_ratio: this value can be calculated directly 
            (see 5.2.2.2(4)) or based on 5.2.2.2(5)
        plan_regular: whether structural system is regular (default is False)
            (see 4.2.3.2 and Table 4.1)
        elev_regular: whether structural system is regular (default is False)
            (see 4.2.3.1(7), 4.2.3.3 and Table 4.1)
    Returns:
        q: behaviour factor
    """
    prop_dict = concrete_structural_system_dict.get(structural_system,{})
    if alpha_ratio is None:
        alpha_ratio = prop_dict.get('alpha_ratio')
        if plan_regular is False:
            alpha_ratio = 0.5 * (1 + alpha_ratio)
    ar_limit = min(1.5, alpha_ratio )
    if ductility_class == 'DCH':
        alphaH = prop_dict.get('alphaH')
        q0 = prop_dict.get('DCH', 1.5) * (ar_limit if alphaH else 1.0)
    elif ductility_class == 'DCM':
        alphaM = prop_dict.get('alphaM')
        q0 = prop_dict.get('DCM', 1.5) * (ar_limit if alphaM else 1.0)
    else:
        q0 = prop_dict.get('DCL', 1.5)
    
    k_w = prop_dict.get('k_w')
    if k_w is None:
        k_w = max(0.5, min((1 + alpha_0)/3, 1.0))
    
    regularity_factor = 1.0 if (elev_regular and plan_regular) else 0.8

    q = max(1.5, q0 * k_w * regularity_factor)
    return q


from math import pi
from typing import Union, Any, Optional, cast


importance_classes = ['I', 'II', 'III', 'IV']

importance_class_dict = {
    'I': 'Buildings of minor importance for public safety, e.g. agricultural buildings, etc.',
    'II':	'Ordinary buildings, not belonging in the other categories.',
    'III':	'Buildings whose seismic resistance is of importance in view of the consequences associated with a collapse, e.g. schools, assembly halls, cultural institutions etc',
    'IV':	'Buildings whose integrity during earthquakes is of vital importance for civil protection, e.g. hospitals, fire stations, power plants, etc.'
}


importance_factor_dict = {
    'I': 0.8,
    'II': 1.0,
    'III': 1.2,
    'IV': 1.4
    }


def get_importance_factor(
        importance_class: str, 
        importance_factor_dict: dict[str, float] = importance_factor_dict
        ) -> float:
    """EC8 importance factor I_E
    
    Returns the EC8 importance factor I_E
    Default values are provided, but can be replaced with NAD values
    by specifying the importance_factor_dict

    Args:
        importance_class: 
        importance_factor_dict:

    Returns:
        importance_factor:  
    """

    importance_factor =  importance_factor_dict.get(importance_class)
    return importance_factor


"""There are two types of spectra: Type 1 and Type 2. If the 
earthquakes that contribute most to the seismic hazard defined 
for the site for the purpose of probabilistic hazard assessment 
have a surface-wave magnitude, Ms, not greater than 5,5, it 
is recommended that the Type 2 spectrum is adopted. 
"""
# Table 3.2 - Type 1
table_3_2_type1_dict = {
    'A': {'S': 1.0, 'T_B': 0.15, 'T_C': 0.4, 'T_D':2.0, }, 
    'B': {'S': 1.2, 'T_B': 0.15, 'T_C': 0.5, 'T_D':2.0, }, 
    'C': {'S': 1.15, 'T_B': 0.2, 'T_C': 0.6, 'T_D':2.0, }, 
    'D': {'S': 1.35, 'T_B': 0.2, 'T_C': 0.8, 'T_D':2.0, }, 
    'E': {'S': 1.4, 'T_B': 0.15, 'T_C': 0.5, 'T_D':2.0, }, 
}

# Table 3.3 - Type 1
table_3_3_type2_dict = {
    'A': {'S': 1.0, 'T_B': 0.05, 'T_C': 0.25, 'T_D':1.2, },
    'B': {'S': 1.35, 'T_B': 0.05, 'T_C': 0.25, 'T_D':1.2, },
    'C': {'S': 1.5, 'T_B': 0.10, 'T_C': 0.25, 'T_D':1.2, },
    'D': {'S': 1.8, 'T_B': 0.10, 'T_C': 0.3, 'T_D':1.2, },
    'E': {'S': 1.6, 'T_B': 0.05, 'T_C': 0.25, 'T_D':1.2, },
}

table_3_4_vert_dict = {
    'type1': {'avg_ag': 0.9, 'T_B': 0.05, 'T_C': 0.15, 'T_D':1.0, },
    'type2': {'avg_ag': 0.45, 'T_B': 0.05, 'T_C': 0.15, 'T_D':1.0, },
}


Table_type_AA = dict[str, dict[str, float]]


def get_eta(zeta=0.05):
    """
    zeta should be a decimal ratio, not as a percentage"""
    return max(0.55, (0.1/(0.05+zeta)))


def get_period_EC8(height_m, building_type='Other'):
    """
    
    """
    C_t_dict = {'Steel_MRF': 0.085, 'RCMRF': 0.075, 'Other':0.05}
    return C_t_dict.get(building_type, 0.05) * height_m ** 0.75



# def get_period_range(T_B:float, T_C:float, T_D:float, T_max:float=10, 
#                     n:int=20, T_add:Any=None) -> tuple[float]:

def get_period_range(T_B:float, 
                     T_C:float, 
                     T_D:float, 
                     T_max:float=10, 
                     n:int=20, 
                     T_add: Union[int, float, 
                           list[Union[int, float], tuple[Union[int, float]]]] = None
                    ) -> tuple[float]:
    """"""
    if any([T_B==0, T_C==0, T_D==0]):
        T_list = []
    else:
        T_list = [0] + \
            [T_B / n * i for i in range(1, n + 1)] + \
            [T_B + (T_C - T_B) / n * i for i in range(1, n + 1)] + \
            [T_C + (T_D - T_C) / n * i for i in range(1, n + 1)] + \
            [T_D + (T_max - T_D) / n * i for i in range(1, n + 1)]
    
    if isinstance(T_add, (int, float)):
        return tuple(sorted(T_list + [T_add]))
    elif isinstance(T_add, (list, tuple)):
        return tuple(sorted(T_list + T_add))
    else:
        return tuple(T_list)
    

def set_period_range(T, T_B=0, T_C=0, T_D=0):
    """Returns a 

    If T is None, the a suitable range of periods will be defined that 
        will plot well (i.e. will hit all the corner points)
    If T is one negative value, then it will be added to the range so 
        that the curve is defined there
    If T is one positive value, then it will be the only value returned
    If T is a collection (tuple or list), then those points will be 
        used to define the range for plotting

    Arg:

    """
    if isinstance(T, (tuple, list)):
        T_tuple = tuple(T)
    elif isinstance(T, (int, float)):
        if T < 0:
            T_tuple = get_period_range(T_B, T_C, T_D, T_add=-T)
        else:
            T_tuple = (T,)
    elif T is None:
        T_tuple = get_period_range(T_B, T_C, T_D)
    else:
        raise ValueError('T should be an int, float or collection of ints and/or floats')
    
    return T_tuple


def S_e(T:Union[float,tuple[float]], 
        a_gR:float, 
        soil_type:str, 
        source_type:int, 
        importance:float=1.0, 
        damp:float=0.05,
        table_3_2_type1_dict: Table_type_AA = table_3_2_type1_dict,
        table_3_3_type2_dict: Table_type_AA = table_3_3_type2_dict
        ) -> Union[float, tuple[tuple[float], tuple[float]]]:
    """EC8 spectral horizontal elastic acceleration 
    
    Returns the EC8 spectral horizontal elastic acceleration (section 3.2.2.2, 
    Table 3.2 & 3.3, eq 3.2 - 3.5).

    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        soil_type: A to E
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        damp: should be as a decimal not as a percentage (default is 0.05 for 5 % damping)
        table_3_2_type1_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 1)
        table_3_3_type2_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 2)
    
    Returns:
        T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 
    """
    if source_type == 1:
        p_dict = table_3_2_type1_dict.get(soil_type, {})
    elif source_type == 2:
        p_dict = table_3_3_type2_dict.get(soil_type, {})
    else:
        return {}
    if p_dict is {}:
        return {}
    
    S, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    eta = 1.0 if damp == 0.05 else get_eta(damp)
    eta_max = 2.5 * eta

    # Create a collection of periods for which S is to be calculated
    T_tuple = set_period_range(T, T_B=T_B, T_C=T_C, T_D=T_D)

    s_list = []
    
    for T in T_tuple:
        if T < T_B:
            s_list.append(a_g * S * (1 + T / T_B * (eta_max - 1.0)))
        elif T <= T_C:
            s_list.append(a_g * S * eta_max)
        elif T <= T_D: 
            s_list.append(a_g * S * eta_max * T_C / T)
        else:
            s_list.append(a_g * S * eta_max * T_C * T_D / (T * T))
    
    T_s_list = (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    return T_s_list
    

def S_De(T, a_gR, soil_type, source_type, importance=1.0, damp=0.05):
    """ elastic displacement response spectrum (EC8, eq 3.7)"""
    return S_e(T, a_gR, soil_type, source_type, 
               importance=importance, damp=damp) * (0.5 * T / pi)**2.0


def S_d(T:Union[float,tuple[float]], 
        a_gR:float, 
        soil_type:str, 
        source_type:int, 
        importance:float=1.0, 
        beta:float=0.2, 
        q:float=1.5, 
        ) -> Union[float, tuple[tuple[float], tuple[float]]]:
    """EC8 spectral horizontal design acceleration

    Returns the EC8 spectral horizontal design acceleration (Table 3.2 & 3.3, eq. 3.13 - 3.16)
    
    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        soil_type: A to E
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        beta: minimum base shear percentage, from NAD (default is 0.2)
        q: 
    
    Returns:
        T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 
    
    source_type is either 1 or 2
    damp should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
    """
    if source_type == 1:
        p_dict = table_3_2_type1_dict.get(soil_type, {})
    elif source_type == 2:
        p_dict = table_3_3_type2_dict.get(soil_type, {})
    else:
        return {}
    if p_dict is {}:
        return {}
    
    S, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    eta = 1.0 / q
    eta_max = 2.5 * eta
    
    # Create a collection of periods for which S is to be calculated
    T_tuple = set_period_range(T, T_B=T_B, T_C=T_C, T_D=T_D)
    
    s_list = []
    
    # Note that `eta_max` is actually 2.5/ `q`
    for T in T_tuple:
        if T < T_B:
            s_list.append(a_g * S * (2/3 + T / T_B * (eta_max - 2/3)))
        elif T <= T_C:
            s_list.append(a_g * S * eta_max)
        elif T <= T_D: 
            s_list.append(a_g * max(beta, S * eta_max * T_C / T))
        else:
            s_list.append(a_g * max(beta, S * eta_max * T_C * T_D / (T * T)))
    
    T_s_list = (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    return T_s_list
    

def S_ve(T, a_gR, source_type, importance=1.0, damp=0.05):
    """EC8 spectral vertical elastic acceleration (S_ve)
    
    EC8 spectral vertical elastic acceleration (S_ve, section 3.2.2.3, 
    EC8 Table 3.4, eq. 3.8 - 3.11)    
    Note that `S` here is actually `a_vg` / `a_g`

    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        damp: inherent damping. Should be as a decimal not as a percentage 
            (default is 0.05 for 5 % damping)
    
    Returns:
        T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 
    """
    if source_type == 1:
        p_dict = table_3_4_vert_dict['type1']
    elif source_type == 2:
        p_dict = table_3_4_vert_dict['type2']
    else:
        p_dict = {}
    
    S, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    eta = 1.0 if damp == 0.05 else get_eta(damp)
    eta_max = 3.0 * eta
    
    # Create a collection of periods for which S is to be calculated
    T_tuple = set_period_range(T, T_B=T_B, T_C=T_C, T_D=T_D)
    
    s_list = []
    
    for T in T_tuple:
        if T < T_B:
            s_list.append(a_g * S * (1 + T / T_B * (eta_max - 1.0)))
        elif T <= T_C:
            s_list.append(a_g * S * eta_max)
        elif T <= T_D: 
            s_list.append(a_g * S * eta_max * T_C / T)
        else:
            s_list.append(a_g * S * eta_max * T_C * T_D / (T * T))
    
    T_s_list = (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    return T_s_list

def S_vd(T, 
         a_gR, 
         source_type, 
         importance=1.0, 
         beta=0.2, 
         q=1.0):
    """EC8 spectral vertical design acceleration
    
    EC8 spectral vertical design acceleration (S_ve, section 3.2.2.5(5-7), 
    EC8 Table 3.4, eq. 3.13 - 3.16)
    
    Note that `S` here is actually `a_vg` / `a_g`
        and `eta` is actually 1 / `q`

    Args:
        T:  is period or collection of periods
        source_type: is either 1 or 2
        q: behaviour factor. Note that this is likely to be different 
            from horizontal q and should not be taken 
            as higher than 1.5 without appropriate analysis (3.2.2.5(6))
        
    """
    if source_type == 1:
        p_dict = table_3_4_vert_dict['type1']
    elif source_type == 2:
        p_dict = table_3_4_vert_dict['type2']
    else:
        return {}
    
    S, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    eta = 1.0 / q
    eta_max = 2.5 * eta
    
    # Create a collection of periods for which S is to be calculated
    T_tuple = set_period_range(T, T_B=T_B, T_C=T_C, T_D=T_D)
    
    s_list = []
    
    for T in T_tuple:
        # Note that `eta_max` is actually 2.5/ `q`
        if T < T_B:
            s_list.append(a_g * S * (2/3 + T / T_B * (eta_max - 2/3)))
        elif T <= T_C:
            s_list.append(a_g * S * eta_max)
        elif T <= T_D: 
            s_list.append(a_g * max(beta, S * eta_max * T_C / T))
        else:
            s_list.append(a_g * max(beta, S * eta_max * T_C * T_D / (T * T)))
    
    return (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    

def get_d_g(a_gR, soil_type, source_type, importance=1.0, damp=0.05):
    """Ground displacement (EC8 eq 3.12)
    source_type is either 1 or 2
    damp should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
    """
    if source_type == 1:
        p_dict = table_3_2_type1_dict.get(soil_type, {})
    elif source_type == 2:
        p_dict = table_3_3_type2_dict.get(soil_type, {})
    else:
        return {}
    if p_dict is {}:
        return {}
    
    S, _, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR

    return 0.025 * a_g * S * T_C * T_D


# Plotting

try:
    import matplotlib.pyplot as plt


    def plot_EC8_spectra(a_gR, site_class, source_type=1, 
                        importance=1.0, q=None, beta=0.2, damp=0.05, T=None, plot_pts=[], 
                        figsize=(12,9), debug=False, **kwargs):
        """Spectral horizontal elastic acceleration (EC8, section 3.2.2.2, 
        Table 3.2 & 3.3, eq 3.2 - 3.5).

        Args:
            T: period or list of periods (self-generated range if T=None)
            a_gR: reference peak ground acceleration on type A ground (in units of g)
            site_class: A to E
            source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
            importance: default is 1.0 (ref:2.1(2-3,P))
            q: force reduction factor (may be an individual value or collection of values)
            damp: should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
            T: a period where a result should be plotted
        
        Returns:
            T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 

        """
        x_limits = kwargs.pop('x_limits', (None, None))
        y_limits = kwargs.pop('y_limits', (None, None))

        Te_range, Se_range = S_e(T, a_gR, site_class, source_type, 
                                 importance=importance, damp=damp)
        legend = ['Elastic Spectrum']
        fig0, ax0 = plt.subplots(figsize=figsize)
        ax0.set_title(f'EC8 Soil Type {site_class}')
        ax0.set_xlabel('Building Period (s)')
        ax0.set_ylabel('Base Shear Coefficient (g)')
        ax0.plot(Te_range, Se_range)
        ax0.grid(True)
        ax0.set_xlim(*x_limits)
        ax0.set_ylim(*y_limits)
        if isinstance(q, (int, float)):
            Td_range, Sd_range = S_d(T, a_gR, site_class, source_type, 
                                     importance=importance, beta=beta, q=q)
            ax0.plot(Td_range, Sd_range) 
            legend.append(f'Design Spectrum (q={q})')
        elif isinstance(q, (tuple, list)):
            for q0 in q:
                Td_range, Sd_range = S_d(T, a_gR, site_class, source_type, 
                                        importance=importance, beta=beta, q=q0)
                ax0.plot(Td_range, Sd_range) 
                legend.append(f'Design Spectrum (q={q0})')
        if len(plot_pts) > 0:
            s_e = S_e(plot_pts, a_gR, site_class, source_type, 
                      importance=importance, damp=damp)
            ax0.plot(plot_pts, s_e)
            if isinstance(q, (int, float)):
                s_d = S_d(plot_pts, a_gR, site_class, source_type, 
                          importance=importance, beta=beta, q=q)
                ax0.plot(plot_pts, s_d)
                if debug: print(f's_d is {s_d} for q={q0}')
            elif isinstance(q, (tuple, list)):
                for q0 in q:
                    s_d = S_d(plot_pts, a_gR, site_class, source_type, 
                            importance=importance, beta=beta, q=q0)
                    ax0.plot(plot_pts, s_d)
                    if debug: print(f's_d is {s_d} for q={q0}')
        
        ax0.legend(legend)

        return fig0

except:
    def plot_EC8_spectra():
        print("MatPlotLib library is not available")
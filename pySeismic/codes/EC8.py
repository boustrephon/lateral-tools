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

from typing import Union, Dict, List, Tuple


def version():
    return 'code: BS EN 1998-1:2004 BS EN 1998-1:2004+A1:2013 \nmodule version: 0.9'

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
        elev_regular: bool = False, 
        plan_regular: bool = False, 
        alpha_0: float = 0.5, 
        alpha_ratio: float = None, 
        debug: bool = False, 
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
        elev_regular: whether structural system is regular (default is False)
            (see 4.2.3.1(7), 4.2.3.3 and Table 4.1)
        plan_regular: whether structural system is regular (default is False)
            (see 4.2.3.2 and Table 4.1)
        alpha_0: prevailing aspect ratio of the walls of
            the structural system - sum(h_wi) / sum(l_wi)
            h_wi is height of wall i
            l_wi is length of wall i
        alpha_ratio: alpha_u / alpha_1 override - this value can be 
            calculated directly (see 5.2.2.2(4)) - default is from 5.2.2.2(5)
    Returns:
        q: behaviour factor
    """
    prop_dict = concrete_structural_system_dict.get(structural_system,{})
    if debug: print(prop_dict.get('Description'))
    if alpha_ratio is None:
        alpha_ratio = prop_dict.get('alpha_ratio')
        if plan_regular is False:
            alpha_ratio = 0.5 * (1 + alpha_ratio)
    ar_limit = min(1.5, alpha_ratio )

    if debug: print(f'- alpha_ratio={alpha_ratio}, ar_limit={ar_limit}, ', end='')
    
    if ductility_class == 'DCH':
        alphaH = prop_dict.get('alphaH')
        q0_basic = prop_dict.get('DCH', 1.5)
        if debug: print(f'alphaH={alphaH}, q0_basic={q0_basic}, ', end='')
        q0 = q0_basic * (ar_limit if alphaH else 1.0)
    elif ductility_class == 'DCM':
        alphaM = prop_dict.get('alphaM')
        q0_basic = prop_dict.get('DCM', 1.5)
        if debug: print(f'alphaM={alphaM}, q0_basic={q0_basic}, ', end='')
        q0 = q0_basic * (ar_limit if alphaM else 1.0)
    else:
        q0_basic = prop_dict.get('DCL', 1.5)
        if debug: print(f'q0_basic={q0_basic}, ', end='')
        q0 = q0_basic
    
    k_w = prop_dict.get('k_w')
    if k_w is None:
        k_w = max(0.5, min((1 + alpha_0)/3, 1.0))
    
    regularity_factor = 1.0 if (elev_regular and plan_regular) else 0.8
    if debug: print(f'regularity_factor={regularity_factor}, ', end='')
    
    q = max(1.5, q0 * k_w * regularity_factor)
    if debug: print(f'q={q:4.3f}')
    
    return q


from math import pi


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
        importance_factor_dict: Dict[str, float] = importance_factor_dict
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


Table_type_AA = Dict[str, Dict[str, float]]


def get_eta(zeta=0.05):
    """
    zeta should be a decimal ratio, not as a percentage"""
    return max(0.55, (0.1/(0.05+zeta)))


def get_period_EC8(height_m, building_type='Other'):
    """
    
    Args:
        height_m: building height in metres
        building_type: (default 'Other')
            'Steel_MRF' - steel moment-resisting frame
            'RCMRF' - reinforced concrete moment-resisting frame
            'Other' - other building types
    
    Returns:
        period: empirical building period (seconds)
    """
    # SMF_group = ('CF1S','CF1B','CFMB')
    CMF_group = ('CF1S','CF1B','CFMB')
    building_type = 'RCMRF' if (building_type in CMF_group) else building_type
    building_type = 'Steel_MRF' if building_type.startswith('S') else building_type
    C_t_dict = {'Steel_MRF': 0.085, 'RCMRF': 0.075, 'Other':0.05}
    return C_t_dict.get(building_type, 0.05) * height_m ** 0.75



# def get_EC8_period_range(T_B:float, T_C:float, T_D:float, T_max:float=10, 
#                     n:int=20, T_add:Any=None) -> tuple[float]:

def get_EC8_period_range(
        T_B:float, 
        T_C:float, 
        T_D:float, 
        T_max:float=10, 
        n:int=20, 
        T_add: Union[
            int, 
            float, 
            List[Union[int, float]], 
            Tuple[Union[int, float]]] = None
    ) -> Tuple[float]:
    """Period values consistent with EC8 response spectra
    
    Generates monotonically increasing values of period that
    are consistent with 

    Periods T_B, T_C & T_D are shown on Figure 3.1 
    "Shape of the elastic response spectrum"

    Args:
        T_B: Corner period at start of constant acceleration zone
        T_C: Corner period between constant acceleration 
            and constant velocity zones
        T_D: Corner period between constant velocity 
            and constant displacement zones
        T_max: Maximum value
        n_int: number of points per zone
        T_add: additional points to add to list
    
    Returns:
        T_tuple: 
    """
    if any([T_B==0, T_C==0, T_D==0]):
        T_list = []
    else:
        T_list = [0]
        
        # Create list of corner points
        corner_list = [T for T in (0, T_B, T_C, T_D) if T < T_max]
        corner_list.append(T_max)
        
        # build up ranges of numbers in each zone
        for T1, T2 in zip(corner_list[:-1], corner_list[1:]):
            T2 = T_max if T2 > T_max else T2
            T_list.extend(
                [T1 + (T2 - T1) / n * i for i in range(1, n + 1)]
            )


    if isinstance(T_add, (int, float)):
        T_tuple = tuple(sorted(T_list + [T_add]))
    elif isinstance(T_add, (list, tuple)):
        T_tuple = tuple(sorted(T_list + T_add))
    else:
        T_tuple = tuple(T_list)
    return T_tuple
    

# def set_period_range(T, T_B=0, T_C=0, T_D=0):
#     """Returns a 

#     If T is None, the a suitable range of periods will be defined that 
#         will plot well (i.e. will hit all the corner points)
#     If T is one negative value, then it will be added to the range so 
#         that the curve is defined there
#     If T is one positive value, then it will be the only value returned
#     If T is a collection (tuple or list), then those points will be 
#         used to define the range for plotting

#     Arg:

#     """
#     if isinstance(T, (tuple, list)):
#         T_tuple = tuple(T)
#     elif isinstance(T, (int, float)):
#         if T < 0:
#             T_tuple = get_EC8_period_range(T_B, T_C, T_D, T_add=-T)
#         else:
#             T_tuple = (T,)
#     elif T is None:
#         T_tuple = get_EC8_period_range(T_B, T_C, T_D)
#     else:
#         raise ValueError('T should be an int, float or collection of ints and/or floats')
    
#     return T_tuple


def S_e(T:Union[float,Tuple[float],List[float]], 
        a_gR:float, 
        site_class:str, 
        source_type:int, 
        importance:float = 1.0, 
        damp:float = 0.05,
        table_3_2_type1_dict: Table_type_AA = table_3_2_type1_dict,
        table_3_3_type2_dict: Table_type_AA = table_3_3_type2_dict,
        ) -> Union[float, Tuple[Tuple[float], Tuple[float]]]:
    """EC8 spectral horizontal elastic acceleration 
    
    Returns the EC8 spectral horizontal elastic acceleration (section 3.2.2.2, 
    Table 3.2 & 3.3, eq 3.2 - 3.5).

    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        site_class: A to E
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
        p_dict = table_3_2_type1_dict.get(site_class, {})
    elif source_type == 2:
        p_dict = table_3_3_type2_dict.get(site_class, {})
    else:
        return {}
    if p_dict is {}:
        return {}
    
    S, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    eta = 1.0 if damp == 0.05 else get_eta(damp)
    eta_max = 2.5 * eta
    
    s_list = []
    
    T = T if isinstance(T, (list, tuple)) else [T]
    for t in T:
        if t < T_B:
            s_list.append(a_g * S * (1 + t / T_B * (eta_max - 1.0)))
        elif t <= T_C:
            s_list.append(a_g * S * eta_max)
        elif t <= T_D: 
            s_list.append(a_g * S * eta_max * T_C / t)
        else:
            s_list.append(a_g * S * eta_max * T_C * T_D / (t * t))
    
    T_s_list = (T, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    
    return T_s_list
    

def S_De(
        T, 
        a_gR, 
        site_class, 
        source_type, 
        importance=1.0, 
        damp=0.05,
        table_3_2_type1_dict: Table_type_AA = table_3_2_type1_dict,
        table_3_3_type2_dict: Table_type_AA = table_3_3_type2_dict,
        ) -> Union[float, Tuple[Tuple[float], Tuple[float]]]:
    """EC8 elastic displacement response spectrum (EC8, eq 3.7)
    
    Returns the EC8 spectral elastic displacement (eq 3.7)

    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        site_class: A to E
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        damp: should be as a decimal not as a percentage (default is 0.05 for 5 % damping)
        table_3_2_type1_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 1)
        table_3_3_type2_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 2)
    
    Returns:
        T_sd_list: tuple of tuple of periods, T & tuple of spectral displacements (m) 
    """
    T_s_list = S_e(T, a_gR, site_class, source_type, 
                importance=importance, damp=damp,
                table_3_2_type1_dict = table_3_2_type1_dict,
                table_3_3_type2_dict = table_3_3_type2_dict,
                )
    if isinstance(T_s_list, (int, float)):
        T_sd_list = T_s_list * (0.5 * T / pi)**2.0
    else:
        sd_list = [s_a * (0.5 * T / pi)**2.0 for T, s_a in zip(*T_s_list)]
        T_sd_list = (T, tuple(sd_list))
    return T_sd_list


def S_d(T:Union[float, List[float], Tuple[float]], 
        a_gR:float, 
        site_class:str, 
        source_type:int, 
        importance:float=1.0, 
        beta:float=0.2, 
        q:float=1.5, 
        table_3_2_type1_dict: Table_type_AA = table_3_2_type1_dict,
        table_3_3_type2_dict: Table_type_AA = table_3_3_type2_dict,
        ) -> Union[float, tuple[tuple[float], tuple[float]]]:
    """EC8 spectral horizontal design acceleration

    Returns the EC8 spectral horizontal design acceleration (Table 3.2 & 3.3, eq. 3.13 - 3.16)
    
    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        site_class: A to E
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        beta: minimum base shear percentage, from NAD (default is 0.2)
        q: force reduction factor which are defined for different structural types:
            - concrete: section 5.2.2 & 5.3.3
            - precast concrete: section 5.11.1.4
            - steel: section 6.3.2
            - composite steel: 7.3.2
            - timber: 8.3
            - masonry: 9.3
            - base isolation: 10.7
        table_3_2_type1_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 1)
        table_3_3_type2_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 2)
    
    Returns:
        T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g)
    """
    if source_type == 1:
        p_dict = table_3_2_type1_dict.get(site_class, {})
    elif source_type == 2:
        p_dict = table_3_3_type2_dict.get(site_class, {})
    else:
        return {}
    if p_dict is {}:
        return {}
    
    S, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    eta = 1.0 / q
    eta_max = 2.5 * eta
    
    T = T if isinstance(T, (list, tuple)) else [T]

    s_list = []
    
    # Note that `eta_max` is actually 2.5/ `q`
    for t in T:
        if t < T_B:
            s_list.append(a_g * S * (2/3 + t / T_B * (eta_max - 2/3)))
        elif t <= T_C:
            s_list.append(a_g * S * eta_max)
        elif t <= T_D: 
            s_list.append(a_g * max(beta, S * eta_max * T_C / t))
        else:
            s_list.append(a_g * max(beta, S * eta_max * T_C * T_D / (t * t)))
    
    T_s_list = (T, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    return T_s_list
    

def S_ve(
        T:Union[float,list[float],tuple[float]], 
        a_gR: float, 
        source_type: int, 
        importance:float = 1.0, 
        damp: float = 0.05,
        table_3_4_vert_dict: Table_type_AA = table_3_4_vert_dict,
        ):
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
        table_3_4_vert_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for source types 1 & 2

    Returns:
        T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 
    """
    if source_type == 1:
        p_dict = table_3_4_vert_dict['type1']
    elif source_type == 2:
        p_dict = table_3_4_vert_dict['type2']
    else:
        p_dict = {}
    
    avg_ag, T_B, T_C, T_D = [p_dict[p] for p in ('avg_ag', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    a_vg = avg_ag * a_g
    eta = 1.0 if damp == 0.05 else get_eta(damp)
    eta_max = 3.0 * eta
    
    T = T if isinstance(T, (list, tuple)) else [T]

    s_list = []
    
    for t in T:
        if t < T_B:
            s_list.append(a_vg * (1 + t / T_B * (eta_max - 1.0)))
        elif t <= T_C:
            s_list.append(a_vg * eta_max)
        elif t <= T_D: 
            s_list.append(a_vg * eta_max * T_C / t)
        else:
            s_list.append(a_vg * eta_max * T_C * T_D / (t * t))
    
    T_s_list = (T, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    return T_s_list


# Create a collection of periods for which S is to be calculated
# T_tuple = get_EC8_period_range(T_B=T_B, T_C=T_C, T_D=T_D, T_add=T)


def S_vd(T: Union[float,list[float],tuple[float]], 
        a_gR: float, 
        source_type: int, 
        importance: float = 1.0, 
        beta:float = 0.2, 
        q:float = 1.0, 
        table_3_4_vert_dict: Table_type_AA = table_3_4_vert_dict,
         ):
    """EC8 spectral vertical design acceleration
    
    EC8 spectral vertical design acceleration (S_ve, section 3.2.2.5(5-7), 
    EC8 Table 3.4, eq. 3.13 - 3.16)
    
    Note that `S` here is actually `a_vg` / `a_g`
        and `eta` is actually 1 / `q`

    Args:
        T: period or collection of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        beta: minimum base shear percentage, from NAD (default is 0.2)
        q: behaviour factor. Note that this is likely to be different 
            from horizontal q and should not be taken 
            as higher than 1.5 without appropriate analysis (3.2.2.5(6))
        table_3_4_vert_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for source types 1 & 2
        
    """
    if source_type == 1:
        p_dict = table_3_4_vert_dict['type1']
    elif source_type == 2:
        p_dict = table_3_4_vert_dict['type2']
    else:
        return None
    
    avg_ag, T_B, T_C, T_D = [p_dict[p] for p in ('avg_ag', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR
    a_vg = avg_ag * a_g
    eta = 1.0 / q
    eta_max = 3.0 * eta
    
    # Create a collection of periods for which S is to be calculated
    # T_tuple = get_EC8_period_range(T_B=T_B, T_C=T_C, T_D=T_D, T_add=T)
    
    T = T if isinstance(T, (list, tuple)) else [T]

    s_list = []
    
    for t in T:
        # Note that `eta_max` is actually 2.5/ `q`
        if t < T_B:
            s_list.append(a_vg * (2/3 + t / T_B * (eta_max - 2/3)))
        elif t <= T_C:
            s_list.append(a_vg * eta_max)
        elif t <= T_D: 
            s_list.append(a_vg * max(beta, eta_max * T_C / t))
        else:
            s_list.append(a_vg * max(beta, eta_max * T_C * T_D / (t * t)))
    
    return (T, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    

def get_d_g(
        a_gR: float, 
        site_class: str, 
        source_type: int, 
        importance: float = 1.0, 
        table_3_2_type1_dict: Table_type_AA = table_3_2_type1_dict,
        table_3_3_type2_dict: Table_type_AA = table_3_3_type2_dict,
            ) -> float:
    """EC8 design horizontal ground displacement (eq 3.12)

    Returns EC8 design horizontal ground displacement in metres (section 3.2.2.4)

    Args:
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        site_class: A to E
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        table_3_2_type1_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 1)
        table_3_3_type2_dict: dictionary of dictionaries 
            containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 2)
    
    Returns:
        d_g: ground displacement (m)
    """
    if source_type == 1:
        p_dict = table_3_2_type1_dict.get(site_class, {})
    elif source_type == 2:
        p_dict = table_3_3_type2_dict.get(site_class, {})
    else:
        return {}
    if p_dict is {}:
        return {}
    
    S, _, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    a_g = importance * a_gR

    d_g = 0.025 * a_g * S * T_C * T_D
    return d_g


# Plotting

try:
    import matplotlib.pyplot as plt


    def plot_EC8_spectra(
            a_gR:float, 
            site_class:str, 
            source_type:int = 1, 
            importance:float = 1.0, 
            q:Union[float, list[float], tuple[float]] = None, 
            beta:float = 0.2, 
            damp:float = 0.05, 
            T:Union[int, float, List[float], Tuple[float]] = None, 
            plot_pts:Union[float, list[float]] = [], 
            figsize:tuple[int, int] = (12,9), 
            table_3_2_type1_dict: Table_type_AA = table_3_2_type1_dict,
            table_3_3_type2_dict: Table_type_AA = table_3_3_type2_dict,
            **kwargs:dict
            ):
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
            plot_pts: adds points to the plot calculating spectral acceleration
                for elastic and every q-value
            figsize: size of plot in points (default is (12,9))
            table_3_2_type1_dict: dictionary of dictionaries 
                containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 1)
            table_3_3_type2_dict: dictionary of dictionaries 
                containing 'S', 'T_B', 'T_C', 'T_D' values for all soil types (source type 2)
    
        
        Returns:
            T_s_list: tuple of tuple of periods, T and tuple of spectral accelerations (g) 

        """
        x_limits = kwargs.pop('x_limits', (0, None))
        y_limits = kwargs.pop('y_limits', (0, None))
        plot_title = kwargs.pop('plot_title', f'EC8 Site Class {site_class}')
        debug = kwargs.pop('debug', False)

        if source_type == 1:
            p_dict = table_3_2_type1_dict.get(site_class, {})
        elif source_type == 2:
            p_dict = table_3_3_type2_dict.get(site_class, {})
        else:
            return {}
        if p_dict is {}:
            return {}
    
        _, T_B, T_C, T_D = [p_dict[p] for p in ('S', 'T_B', 'T_C', 'T_D', )]
    
        T_tuple = get_EC8_period_range(T_B=T_B, T_C=T_C, T_D=T_D, T_add=T)
        Te_range, Se_range = S_e(T_tuple, a_gR, site_class, source_type, 
                                 importance=importance, damp=damp)
        legend = ['Elastic Spectrum']
        fig0, ax0 = plt.subplots(figsize=figsize)
        ax0.set_title(f'Spectral Accelerations')
        ax0.set_xlabel('Building Period (s)')
        ax0.set_ylabel('Base Shear Coefficient (g)')
        ax0.plot(Te_range, Se_range)
        ax0.set_xlim(*x_limits)
        ax0.set_ylim(*y_limits)
        fig0.suptitle(plot_title)
        ax0.minorticks_on()
        ax0.grid(True, which='major', linestyle='-', linewidth=0.5)
        ax0.grid(True, which='minor', linestyle='--', linewidth=0.25)
        # ax0.tick_params(axis='x', rotation=45)

        
        
        if isinstance(q, (int, float)):
            Td_range, Sd_range = S_d(T_tuple, a_gR, site_class, source_type, 
                                     importance=importance, beta=beta, q=q)
            ax0.plot(Td_range, Sd_range) 
            legend.append(f'Design Spectrum (q={q})')
        elif isinstance(q, (tuple, list)):
            for q0 in q:
                Td_range, Sd_range = S_d(T_tuple, a_gR, site_class, source_type, 
                                        importance=importance, beta=beta, q=q0)
                ax0.plot(Td_range, Sd_range) 
                legend.append(f'Design Spectrum (q={q0})')
        
        ax0.legend(legend)

        T = T if isinstance(T, (list, tuple)) else [T]
        
        for t in T:
            s_e = S_e(t, a_gR, site_class, source_type, importance=1.0, damp=damp)
            # print(f's_e is {s_e:.3f}')
            ax0.plot(t, s_e, marker='o', mec = 'C0', mfc = 'C0')
            if isinstance(q, (int, float)):
                s_d = S_d(t, a_gR, site_class, source_type, importance=1.0, beta=beta, q=q)
                ax0.plot(t, s_d, marker='o', mec = 'C1', mfc = 'C1')
            elif isinstance(q, (tuple, list)):
                for i, q0 in enumerate(q):
                    s_d = S_d(t, a_gR, site_class, source_type, importance=1.0, beta=beta, q=q0)
                    ax0.plot(t, s_d, marker='o', mec = f'C{i+1}', mfc = f'C{i+1}')
                # print(f's_d is {s_d:.3f}')
    
        # if len(plot_pts) > 0:
        #     s_e = S_e(plot_pts, a_gR, site_class, source_type, 
        #               importance=importance, damp=damp)
        #     ax0.plot(plot_pts, s_e)
        #     if isinstance(q, (int, float)):
        #         s_d = S_d(plot_pts, a_gR, site_class, source_type, 
        #                   importance=importance, beta=beta, q=q)
        #         ax0.plot(plot_pts, s_d)
        #         if debug: print(f's_d is {s_d} for q={q0}')
        #     elif isinstance(q, (tuple, list)):
        #         for q0 in q:
        #             s_d = S_d(plot_pts, a_gR, site_class, source_type, 
        #                     importance=importance, beta=beta, q=q0)
        #             ax0.plot(plot_pts, s_d)
        #             if debug: print(f's_d is {s_d} for q={q0}')
        #     else:
        #         pass        
        
        return fig0

except:
    def plot_EC8_spectra(**kwargs):
        print("MatPlotLib library is not available")
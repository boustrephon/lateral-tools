"""

"""

from math import pi
from typing import Union, Any, Optional, cast
from typing import List, Set, Dict, Tuple, Optional

"""There are two types of spectra: Type 1 and Type 2. If the 
earthquakes that contribute most to the seismic hazard defined 
for the site for the purpose of probabilistic hazard assessment 
have a surface-wave magnitude, Ms, not greater than 5,5, it 
is recommended that the Type 2 spectrum is adopted. 
"""
# Type 1
type1_dict = {
    'A': {'S': 1.0, 'T_B': 0.15, 'T_C': 0.4, 'T_D':2.0, }, 
    'B': {'S': 1.2, 'T_B': 0.15, 'T_C': 0.5, 'T_D':2.0, }, 
    'C': {'S': 1.15, 'T_B': 0.2, 'T_C': 0.6, 'T_D':2.0, }, 
    'D': {'S': 1.35, 'T_B': 0.2, 'T_C': 0.8, 'T_D':2.0, }, 
    'E': {'S': 1.4, 'T_B': 0.15, 'T_C': 0.5, 'T_D':2.0, }, 
}

# Type 1
type2_dict = {
    'A': {'S': 1.0, 'T_B': 0.05, 'T_C': 0.25, 'T_D':1.2, },
    'B': {'S': 1.35, 'T_B': 0.05, 'T_C': 0.25, 'T_D':1.2, },
    'C': {'S': 1.5, 'T_B': 0.10, 'T_C': 0.25, 'T_D':1.2, },
    'D': {'S': 1.8, 'T_B': 0.10, 'T_C': 0.3, 'T_D':1.2, },
    'E': {'S': 1.6, 'T_B': 0.05, 'T_C': 0.25, 'T_D':1.2, },
}

vert_dict = {
    'type1': {'avg_ag': 0.9, 'T_B': 0.05, 'T_C': 0.15, 'T_D':1.0, },
    'type2': {'avg_ag': 0.45, 'T_B': 0.05, 'T_C': 0.15, 'T_D':1.0, },
}


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

def get_period_range(T_B:float, T_C:float, T_D:float, T_max:float=10, 
                     n:int=20, T_add:Any=None) -> Tuple[float]:
    
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
    if isinstance(T, (tuple,list)):
        T_tuple = tuple(T)
    elif isinstance(T, (int,float)):
        if T < 0:
            T_tuple = get_period_range(T_B, T_C, T_D, T_add=-T)
        else:
            T_tuple = (T,)
    elif T is None:
        T_tuple = get_period_range(T_B, T_C, T_D)
    else:
        raise ValueError('T should be an int, float or collection of ints and floats')
    
    return T_tuple


def S_e(T:Union[float,Tuple[float]], a_gR:float, soil_type:str, 
        source_type:int, importance:float=1.0, damp:float=0.05
        ) -> Union[float,Tuple[Tuple[float],Tuple[float]]]:
    """Spectral horizontal elastic acceleration (EC8, section 3.2.2.2, 
    Table 3.2 & 3.3, eq 3.2 - 3.5).

    Args:
        T: period or list of periods
        a_gR: reference peak ground acceleration on type A ground (in units of g)
        soil_type: A to E
        source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
        importance: default is 1.0 (ref:2.1(2-3,P))
        damp: should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
    
    Returns:
        T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 
    """
    if source_type == 1:
        p_dict = type1_dict.get(soil_type, {})
    elif source_type == 2:
        p_dict = type2_dict.get(soil_type, {})
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
    
    return (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    

def S_De(T, a_gR, soil_type, source_type, importance=1.0, damp=0.05):
    """ elastic displacement response spectrum (EC8, eq 3.7)"""
    return S_e(T, a_gR, soil_type, source_type, 
               importance=importance, damp=damp) * (0.5 * T / pi)**2.0


def S_d(T, a_gR, soil_type, source_type, importance=1.0, beta=0.2, q=1.5, n_pt=20):
    """Spectral horizontal design acceleration (EC8, Table 3.2 & 3.3, eq. 3.13 - 3.16)
    source_type is either 1 or 2
    damp should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
    """
    if source_type == 1:
        p_dict = type1_dict.get(soil_type, {})
    elif source_type == 2:
        p_dict = type2_dict.get(soil_type, {})
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
    
    return (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    

def S_ve(T, a_gR, soil_type, source_type, importance=1.0, damp=0.05):
    """Spectral vertical elastic acceleration (S_ve, section 3.2.2.3, 
    EC8 Table 3.4, eq. 3.8 - 3.11)
    source_type is either 1 or 2
    damp should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
    Note that `S` here is actually `a_vg` / `a_g`
    """
    if source_type == 1:
        p_dict = vert_dict['type1']
    elif source_type == 2:
        p_dict = vert_dict['type2']
    else:
        return {}
    
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
    
    return (T_tuple, tuple(s_list)) if len(s_list) > 1 else s_list[0]
    

def S_vd(T, a_gR, soil_type, source_type, importance=1.0, beta=0.2, q=1.0):
    """Spectral vertical design acceleration (S_ve, section 3.2.2.5(5-7), 
    EC8 Table 3.4, eq. 3.13 - 3.16)
    T is period or collection of periods
    source_type is either 1 or 2
    q is likely to be different from horizontal q and should not be taken 
    as higher than 1.5 withough appropriate analysis
    Note that `S` here is actually `a_vg` / `a_g`
        and `eta` is actually 1 / `q`"""
    if source_type == 1:
        p_dict = vert_dict['type1']
    elif source_type == 2:
        p_dict = vert_dict['type2']
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
        p_dict = type1_dict.get(soil_type, {})
    elif source_type == 2:
        p_dict = type2_dict.get(soil_type, {})
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
            soil_type: A to E
            source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
            importance: default is 1.0 (ref:2.1(2-3,P))
            damp: should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
        
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
        if q is not None:
            Td_range, Sd_range = S_d(T, a_gR, site_class, source_type, 
                                     importance=importance, beta=beta, q=q)
            ax0.plot(Td_range, Sd_range) 
            legend.append(f'Design Spectrum (q={q})')
        if len(plot_pts) > 0:
            s_e = S_e(plot_pts, a_gR, site_class, source_type, 
                      importance=importance, damp=damp)
            ax0.plot(plot_pts, s_e)
            if q is not None:
                s_d = S_d(plot_pts, a_gR, site_class, source_type, 
                          importance=importance, beta=beta, q=q)
                ax0.plot(plot_pts, s_d)
                if debug: print(f's_d is {s_d}')
        
        
        ax0.legend(legend)

        return fig0

except:
    def plot_EC8_spectra():
        print("MatPlotLib library is not available")
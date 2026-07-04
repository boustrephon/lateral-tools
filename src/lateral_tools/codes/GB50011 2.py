"""
These are some old functions that are mainly based on the 2001 version of the code, 
since there is a readily available translation.

GB 18306-2015 "Seismic Ground Motion Parameter Zonation Map of China" 
(which replaced older versions) provides the zonation map using these 
parameters, and GB 50011-2010 "Code for Seismic Design of Buildings" 
uses the intensity levels defined in these standards to determine design 
parameters (e.g., 0.10g, 0.15g, 0.20g, 0.30g, 0.40g).
"""

from typing import Sequence, Any
import numpy as np
from numpy import floating
from numpy.typing import NDArray #, ArrayLike


# GB/T 17742-2008: "Chinese Seismic Intensity Scale"
intensity_dict = {
    'V':{'median': 0.319, 'min': 0.223, 'max': 0.456, 'median_PGV': 0.0259},
    'VI':{'median': 0.653, 'min': 0.457, 'max': 0.936, 'median_PGV': 0.0557},
    'VII':{'median': 1.35, 'min': 0.937, 'max': 1.94, 'median_PGV': 0.120},
    'VIII':{'median': 2.79, 'min': 1.95, 'max': 4.01, 'median_PGV': 0.258},
    'IX':{'median': 5.77, 'min': 4.02, 'max': 8.30, 'median_PGV': 0.555},
    'X':{'median': 11.9, 'min': 8.31, 'max': 17.2, 'median_PGV': 1.19},
}

def get_intensity(value, debug = False):
    """
    Convert a PGA value to the corresponding intensity based on the GB/T 17742-2008 scale
    
    Args:
        value (float): PGA value
    
    Returns:
        str: The intensity corresponding to the given PGA value
    """
    intensity = 'High'
    for key, low, high in [(k, v['min'] / 9.81, v['max'] / 9.81) 
        for k, v in intensity_dict.items()]:
        if debug: 
            print(f'{key}, {value:.4f} | {value < low} | {low <= value < high} | {low:.4f}, {high:.4f}')
        if value < low:
            intensity = 'Low'
        elif low <= value < high:
            intensity = key
            break
        else:
            pass
    return intensity

# China Code
# Table 5.1.4-1
seismic_zone_dict = {
    '6': {'intensity':'6 (0.05g)', 'PGA': 0.05, 'alphamax_L1': 0.04, 'alphamax_L3': None}, 
    '7': {'intensity':'7 (0.10g)', 'PGA': 0.10, 'alphamax_L1': 0.08, 'alphamax_L3': 0.5},
    '7+': {'intensity':'7 (0.15g)', 'PGA': 0.15, 'alphamax_L1': 0.12, 'alphamax_L3': 0.72},
    '8': {'intensity':'8 (0.20g)', 'PGA': 0.20, 'alphamax_L1': 0.16, 'alphamax_L3': 0.9},
    '8+': {'intensity':'8 (0.30g)', 'PGA': 0.30, 'alphamax_L1': 0.24, 'alphamax_L3': 1.2},
    '9': {'intensity':'9 (0.40g)', 'PGA': 0.40, 'alphamax_L1': 0.32, 'alphamax_L3': 1.4},
    }

fortification_dict = {
    '6':0.05, '7':0.10, '7+':0.15, '8':0.20, '8+':0.30, '9':0.40, 
}


def get_fortification(intensity):
    return fortification_dict.get(intensity, None)


def get_fortification_dict(intensity, intensity_dict=seismic_zone_dict):
    return intensity_dict.get(intensity, {})


# Table 5.1.4-2
T_G_dict_GB50011_2001 = {
    1:{'I': 0.25, 'II': 0.35, 'III': 0.45, 'IV': 0.65},
    2:{'I': 0.30, 'II': 0.40, 'III': 0.55, 'IV': 0.75},
    3:{'I': 0.35, 'II': 0.45, 'III': 0.65, 'IV': 0.90}
    }

# Table 5.1.4-2
T_G_dict_GB50011_2010 = {
    1:{'I0': 0.20, 'I1': 0.25, 'II': 0.35, 'III': 0.45, 'IV': 0.65},
    2:{'I0': 0.25, 'I1': 0.30, 'II': 0.40, 'III': 0.55, 'IV': 0.75},
    3:{'I0': 0.30, 'I1': 0.35, 'II': 0.45, 'III': 0.65, 'IV': 0.90}
    }

# Table 3.0.5
T_G_dict_GB50023_2009 = {
    1:{'I': 0.20, 'II': 0.30, 'III': 0.40, 'IV': 0.65},
    2:{'I': 0.20, 'II': 0.30, 'III': 0.40, 'IV': 0.65},
    3:{'I': 0.25, 'II': 0.40, 'III': 0.55, 'IV': 0.85}
}

def T_G_freq_calc(MD_Group:int, soil_class:str, 
                  T_G_dict = T_G_dict_GB50011_2001)->float:
    """GB50011-2001 Table 5.1.4-2

    Args:
        MD_Group (int): _description_
        soil_class (str): _description_
        T_G_dict (_type_, optional): _description_. Defaults to T_G_dict_GB50011_2001.

    Raises:
        Exception: _description_

    Returns:
        float: spectral acceleration

    Example:
        >>> T_G_freq_calc(2, 'IV', T_G_dict_2001)
        0.75
    """
    T_code = T_G_dict.get(MD_Group,{}).get(soil_class, None)
    if T_code is None:
        raise ValueError(f"MD_Group: {MD_Group}, soil_class: {soil_class}")
    return T_code


def T_G_rare(MD_Group:int, soil_class:str, T:float)->float:
    """_summary_

    Args:
        MD_Group (int): _description_
        soil_class (str): _description_
        T (float): _description_

    Returns:
        float: _description_
    """
    return T_G_freq_calc(MD_Group, soil_class) + (0.05 if T > 2 else 0.05)


alpha_max_dict = {
    1: {'6':18, '7':35, '7+':55, '8':70, '8+':110, '9':140, },
    #2: {'6':0.04, '7':None, '7+':None, '8':None, '8+':None, '9':None, },
    3: {'6':None, '7':220, '7+':310, '8':400, '8+':510, '9':620, },
}

alpha_dict = {
    1: {'6':0.04, '7':0.08, '7+':0.12, '8':0.16, '8+':0.24, '9':0.32, },
    #2: {'6':0.04, '7':0.08, '7+':0.12, '8':0.16, '8+':0.24, '9':0.32, },
    3: {'6':None, '7':0.50, '7+':0.72, '8':0.90, '8+':1.2, '9':1.4, },
}

# def alpha_calc(zeta = 0.05):
#     # 5.1.5-1
#     gamma = 0.9 + (0.05 - zeta) / (0.5 + 5 * zeta)
#     # 5.1.5-2
#     eta_1 = max(0, 0.02 + 0.125 * (0.05 - zeta))
#     # 5.1.5-3
#     eta_2 = max(0.55, 1 + (0.05 - zeta) / (0.06 + 1.7 * zeta))

# GB 50011 Table 5.1.4-1: seismic influence coefficient alpha_max
GB_ALPHA_FREQUENT = {6: 0.04, 7: 0.08, 8: 0.16, 9: 0.32}
GB_ALPHA_RARE = {6: 0.28, 7: 0.50, 8: 0.90, 9: 1.40}


def GB_get_alpha_max(
    intensity: int = 7,
    acceleration: float = 0.10,
    level: str = "rare",
) -> float:
    """Return the seismic influence coefficient alpha_max for GB 50011.

    Parameters
    ----------
    intensity : int
        Seismic fortification intensity (6, 7, 8, or 9).
    acceleration : float
        Design basic acceleration (g), e.g. 0.10, 0.15, 0.20.
    level : str
        ``'frequent'`` (多遇), ``'fortification'`` (设防), or ``'rare'`` (罕遇).

    Returns
    -------
    float
        The alpha_max value (dimensionless).

    Notes
    -----
    For the fortification level, alpha_max is taken as the larger of
    *acceleration* × 2.25 and the frequent alpha_max × 2.5.
    """
    if level == "frequent":
        return GB_ALPHA_FREQUENT.get(intensity, 0.08)
    elif level == "fortification":
        return max(acceleration * 2.25, GB_ALPHA_FREQUENT.get(intensity, 0.08) * 2.5)
    else:  # rare
        return GB_ALPHA_RARE.get(intensity, 0.50)


def GB_design_spectrum(
    T: np.ndarray,
    alpha_max: float,
    site_class: str = "II",
    group: int = 1,
    zeta: float = 0.05,
    T_G_dict: dict = T_G_dict_GB50011_2010,
) -> np.ndarray:
    """Convenience wrapper around GB_spectrum returning Sa in m/s\\u00b2.

    Parameters are identical to :func:`GB_spectrum` except that the
    output is multiplied by *g* = 9.81 m/s\\u00b2.
    """
    return GB_spectrum(alpha_max, T, site_class, group, zeta, T_G_dict) * 9.81


def GB_spectrum(alpha_max:float, T:Sequence[float]|NDArray[floating[Any]], 
                site_class:str='I0', 
                group:int=1, zeta:float=0.05,
                T_G_dict:dict=T_G_dict_GB50011_2010)->NDArray[Any]:
    """Returns GB spectral acceleration Sa(g) for given PGA(g) and period array T(s).
    GB soil class I0 - IV, Group 1 spectrum.
    
    Args:
        alpha_max (float): _description_
        T (Sequence[float] | NDArray[floating[Any]]): _description_
        site_class (str, optional): _description_. Defaults to 'I0'.
        group (int, optional): _description_. Defaults to 1.
        zeta (float, optional): _description_. Defaults to 0.05.
        T_G_dict (dict, optional): Soil characteristic period. Defaults to T_G_dict_2010.
            Options are:
                T_G_dict_GB50011_2001
                T_G_dict_GB50011_2010
                T_G_dict_GB50023_2009
    
    Raises:
        ValueError: _description_

    Returns:
        NDArray[Any]: _description_
    """
    gamma = 0.9 + (0.05 - zeta) / (0.3 + 6.0 * zeta)
    eta_1 = max(0.0, 0.02 + (0.05 - zeta) / (4.0 + 32.0 * zeta))
    eta_2 = max(0.55, 1.0 + (0.05 - zeta) / (0.08 + 1.6 * zeta))
    
    T_g = T_G_dict[group][site_class]
    # pga = alpha_max * 0.45
    Sa = np.zeros_like(T)
    for i, t in enumerate(T):
        if t <= 0.1:
            Sa[i] = (0.45 + (eta_2 - 0.45) * 10 * t) * alpha_max
        elif t <= T_g:
            Sa[i] = eta_2 * alpha_max
        elif t <= 5 * T_g:
            Sa[i] = (T_g / t)**gamma * eta_2 * alpha_max
        elif t <= 6.0:
            Sa[i] = (eta_2 * 0.2**gamma - eta_1 * (t - 5 * T_g)) * alpha_max
        else:
            raise ValueError(f'Period ({t:.4f}) is greater than 6.0')
    return Sa

def EC8_spectrum(pga:float, T:Sequence[float],  
             S:float=1.2, TB:float=0.15, 
             TC:float=0.5, TD:float=2.0)->NDArray[Any]:
    """
    Returns spectral acceleration Sa(g) for given PGA(g) and period array T(s).
    Eurocode 8 Type 1 spectrum, ground type B.
    
    Args:
        pga (float): _description_
        T (Sequence[float]): _description_
        S (float, optional): _description_. Defaults to 1.2.
        TB (float, optional): _description_. Defaults to 0.15.
        TC (float, optional): _description_. Defaults to 0.5.
        TD (float, optional): _description_. Defaults to 2.0.

    Returns:
        NDArray[Any]: _description_
    """
    Sa = np.zeros_like(T)
    for i, t in enumerate(T):
        if t <= TB:
            Sa[i] = pga * S * (1 + (t/TB)*(2.5 - 1))
        elif t <= TC:
            Sa[i] = 2.5 * pga * S
        elif t <= TD:
            Sa[i] = 2.5 * pga * S * (TC / t)
        else:
            Sa[i] = 2.5 * pga * S * (TC * TD / t**2)
    return Sa


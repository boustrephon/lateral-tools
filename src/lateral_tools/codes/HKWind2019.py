"""2019 Wind Code

Note that a building is considered to be one with resonant dynamic response if:
  a) the height exceeds five times the least dimension
  b) the height of the building is greater than 100m
  Unless the fundamental natural frequency of the building can be shown to be
  greater than 1.0Hz
"""

from math import log, exp #, log10, pi, sin, cos, tan, acos, asin, atan, atan2
from itertools import accumulate
import matplotlib.pyplot as plt
import doctest
from sys import path as sys_path
from os.path import abspath, join as path_join, dirname

# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

# printt, interp, storey_force_coeffs, storey_force_distribution, storey_forces, storey_shear_moment, shear_coords
from lateral_tools.math_utils import list_interp
# from lateral_tools.math_utils import *
from lateral_tools.analysis import *
from lateral_tools.risk_utils import *
from lateral_tools.print_utils import *
# from lateral_tools.codes.HKWind2004 import get


# ===============================
# ====== General Functions ======
# ===============================

def get_storey_shear(s_forces: list[float], ascending:bool=True) -> list[float]:
    """Returns a list of storey shears

    Args:
        s_forces (list[float]): list of storey forces (ascending or descending)
        ascending (bool, optional): forces are ordered with rising storey levels. 
            Defaults to True.

    Returns:
        list[float]: List of storey shears
    """
    if ascending:
        return list(accumulate(s_forces[::-1]))[::-1]
    else:
        return list(accumulate(s_forces))


def get_storey_moment(storey_shears: list[float], 
                      storey_heights: list[float], 
                      ascending:bool=True) -> list[float]:
    """Returns a list of storey moments from storey shears and storey heights

    Args:
        storey_shears (list[float]): List of storey shears
        storey_heights (list[float]): List of storey heights
        ascending (bool, optional): forces are ordered with rising storey levels. 
            Defaults to True.

    Returns:
        list[float]: List of storey moments
    """
    if ascending:
        return list(accumulate(V * h 
                               for V, h 
                               in zip(storey_shears[::-1], 
                                      storey_heights[::-1])))[-2::-1] + [0]
    else:
        return [0] + (accumulate(V * h 
                                 for V, h 
                                 in zip(storey_shears, 
                                        storey_heights)))[:-1]


def gen_ff_heights(overall_height: float, h_list: list[float]) -> list[float]:
    """Generates storey heights to match the overall height. As a starter
    it takes any storey heights you provide and repeats the last value to
    extend the list.

    Args:
        overall_height (float): Target for total height (sum of storey heights)
        h_list (list[float]): bottom storey heights as a list - the last in the list
                will be repeated until the building height is reached.

    Returns:
        list[float]: Storey heights that sum to match the overall height.
    """
    # [5.5] + [3.5 for i in range(int((100-5.5)//3.5))]
    i_range = range(int((overall_height - sum(h_list[:-1])) // h_list[-1]))
    ff_heights = h_list[:-1] + [h_list[-1] 
                        for _ in i_range]
    return ff_heights


def gen_elevs(list_of_storey_heights: list[float], 
              ascending:bool=True) -> list[float]:
    """Returns list of storey elevations with same length as storey heights list
    if a base storey height of zero is provided, otherwise output list
    has one additional member.

    Args:
        list_of_storey_heights (list[float]): List of storey heights
        ascending (bool, optional): Whether storey heights are 
                in ascending order. Defaults to True.

    Returns:
        list[float]: storey elevations in the same order as the input storey heights
    """
    if not ascending:
        list_of_storey_heights = list(list_of_storey_heights)[::-1]
    if list_of_storey_heights[0] == 0:
        storey_elevs = list(accumulate(list_of_storey_heights))
    else:
        storey_elevs = [0] + list(accumulate(list_of_storey_heights))
    return storey_elevs if ascending else storey_elevs[::-1]


def gen_trib_heights(list_of_storey_heights: list[float], 
                     ascending:bool=True) -> list[float]:
    """Returns list of tributary heights with same length as storey heights list
    if a base storey height of zero is provided, otherwise output list
    has one additional member.

    Args:
        list_of_storey_heights (list[float]): List of storey heights
        ascending (bool, optional): Whether storey heights 
                are in ascending order. Defaults to True.

    Returns:
        list[float]: Tributary heights in the same order as the input storey heights
    """
    if not ascending:
        list_of_storey_heights = list(list_of_storey_heights)[::-1]
    if list_of_storey_heights[0] == 0:
        trib_heights = [round(0.5 * (a + b), 10) 
                        for b, a in zip(list_of_storey_heights,
                                        list_of_storey_heights[1:] + [0.0])]
    else:
        trib_heights = [round(0.5 * (a + b), 10) 
                        for b, a in zip([0.0] + list_of_storey_heights,
                                        list_of_storey_heights + [0.0])]
    return trib_heights if ascending else trib_heights[::-1]


def gen_trib_areas(list_of_storey_heights, list_of_storey_widths):
    """
    
    """
    if isinstance(list_of_storey_widths, list):
        h_len = len(list_of_storey_heights)
        w_len = len(list_of_storey_widths)
        if h_len != w_len:
            err_msg = 'number of storey heights ' + \
                f'({h_len}) does not match number of storey widths({w_len})'
            raise ValueError(err_msg)
        # TODO - fix 'list_of_zip'
        list_of_storey_areas = [h * w 
                                for h, w in zip(list_of_storey_heights, 
                                                list_of_storey_widths)]
        if list_of_storey_heights[0] == 0:
            return [0.5 * (a + b) for b, a in zip(list_of_storey_areas,
                                                list_of_storey_areas[1:] + [0.0])]
        else:
            return [0.5 * (a + b) for b, a in zip([0.0] + list_of_storey_areas,
                                            list_of_storey_areas + [0.0])]
    #elif isinstance(list_of_storey_widths,(float,int)):
    else:
        w = list_of_storey_widths
        if list_of_storey_heights[0] == 0:
            return [0.5 * w * (a + b) for b, a in zip(list_of_storey_heights,
                                                list_of_storey_heights[1:] + [0.0])]
        else:
            return [0.5 * w * (a + b) for b, a in zip([0.0] + list_of_storey_heights,
                                          list_of_storey_heights + [0.0])]



# =================================
# ====== Wind Code Functions ======
# =================================

def get_v_z_mean(z:float, z_g:float=500, 
                 alpha:float=0.11, v_g_mean:float=56.6
                 ) -> float:
    """hourly-mean wind velocity at height z (eq 1.5??)

    Args:
        z_g (float, optional): gauge height. Defaults to 500m.
        alpha (float, optional): power coefficient. Defaults to 0.11.
        v_g_mean (float, optional): hourly-mean velocity at gauge height. Defaults to 56.6.

    Returns:
        float: hourly-mean velocity at height z
    """
    z = min(max(z, 5), 500)
    return 1.05 * v_g_mean * (z / z_g)**alpha


def get_v_z(z, z_g=500, alpha=0.11, v_g_mean=56.6, g_v=3.7, I_g=0.087):
    """design 3 second gust wind speed v_z at height z (eq 1.2, 1.3, 1.5 combined)
    z - elevation in metres
    z_g - reference height in metres
    alpha - the power coefficient defining the wind gradient
    v_g_mean - hourly-mean wind velocity at reference height in m/s
    g_v - the peak factor, normally taken between 3.4 and 3.7
    I_g - 
    """
    z = min(max(z, 5), 500)
    if z > 0:
        return 1.05 * v_g_mean * (z / z_g)**alpha * (1 + g_v * I_g * (z / z_g)**(-alpha))
    else:
        return 1.05 * v_g_mean * g_v * I_g


def get_q_z(v_z:float, rho:float=1.2) -> float:
    """Returns design wind pressure based on wind velocity

    Args:
        v_z (float): velocity at elevation z in m/s
        rho (float, optional): air density in kg / m^3. Defaults to 1.2.

    Returns:
        float: design wind pressure in kPa
    """
    return 0.0005 * rho * v_z**2


def get_G(h:float, b:float, 
          n_a:float=None, zeta:float = 0.02, 
          g_v:float=3.7) -> float:
    """Dynamic Magnification Factor, G (Appendix F)
    
    Args:
        h (float): building height in metres
        b (float): building breadth (dimension perpendicular to the wind) in metres
        n_a (float, optional): building frequency in Hz. Defaults to None.
        zeta (float, optional): _description_needed_. Defaults to 0.02.
        g_v (float, optional): _description_needed_. Defaults to 3.7.

    Returns:
        float: Dynamic Magnification Factor, G
    """
    # h = min(max(h, 100), 200) # should h be limited?
    V_h_mean = get_v_z_mean(h)
    n_a = n_a if n_a else (46 / h)
    I_h = 0.1055 * (h / 90)**(-0.11)
    g_f = (2 * log(3600 * n_a))**0.5
    L_h = 1000 * (h / 10)**0.25
    B = 1 / (1 + (36 * h**2 + 64 * b**2)**0.5 / L_h)
    N = n_a * L_h / V_h_mean
    E = 0.47 * N / (2 + N**2)**(5/6)
    S = 1 / (1 + 3.5 * n_a * h / V_h_mean) / (1 + 4.0 * n_a * b / V_h_mean)
    return 1 + 2 * I_h * (g_v**2 * B + g_f**2 * S * E / zeta)**0.5


def get_Q_oz_table_3_1(Ze:float) -> float:
    """Returns Wind Reference Pressure Q_oz from table 3-1 (HK Wind 2019).
    Preferably use the function 'get_Q_oz' which uses equation 3-2

    Args:
        Ze (float): Ze is the effective height for open exposure,
            conservatively taken as the height above ground level

    Raises:
        ValueError: Seek specialist advice for Ze > 500m

    Returns:
        float: Q_oz from table 3-1 (HK Wind 2019)
    """
    # Table 3.1 data
    Zes = (2.5,5,10,20,30,50,75,100,150,200,250,300,400,500)
    Q_ozs = (1.59,1.77,1.98,2.21,2.36,2.56,2.73,2.86,3.05,3.20,3.31,3.41,3.57,3.70)
    
    if Ze <= 500: 
        err_msg = f'Effective height (Ze) is {Ze:4.2f} > 500m : Seek specialist advice'
        raise ValueError(err_msg)
    Q_oz = list_interp(Ze, Zes, Q_ozs, extend=True) 
    return Q_oz


def get_Q_oz(Ze:float) -> float:
    """Wind reference pressure  (Equation 3-2, HK Wind 2019).
    This is an alternative to the table lookup function 'get_table_3_1'

    Args:
        Ze (float): the effective height for open exposure,
            conservatively taken as the height above ground level.

    Raises:
        ValueError: If Ze is greater than 500m

    Returns:
        float: wind reference pressure Q_oz from Eq 3-2, HK Wind 2019
    """
    if Ze > 500: 
        err_msg = f'Effective height (Ze) is {Ze:4.2f} > 500m : Seek specialist advice'
        raise ValueError(err_msg)
    Q_oz =  3.7 * ((Ze / 500)**0.16) # Eq 3-2
    return Q_oz


def get_Q_z(Q_oz:float, S_t: float = 1.0, S_theta: float = 0.85) -> float:
    """design wind pressure, Q_z
    at an elevation Z defined by Eqn 3-1

    Args:
        Q_oz (float): defined in caluse 3.2
        S_t (float): topography factor (see section 3.4)
        S_theta (float): directionality factor (see Appendix A1)

    Returns:
        float: design wind pressure, Q_z
    """
    Q_z = Q_oz * S_t * S_theta
    return Q_z


def get_I_oz(Ze: float) -> float:
    """Turbulence intensity at effective height Ze (Eqn 3-4)

    Args:
        Ze (float): Ze is the effective height for open exposure,
            conservatively taken as the height above ground level

    Raises:
        ValueError: Invalid for Ze > 500m. Seek specialist advice.

    Returns:
        float: _description_
    """
    if Ze > 500: 
        raise ValueError('Seek specialist advice for Ze > 500m')
    return 0.087 * (Ze / 500)**-0.11 # Eq 3-3


def get_C_f(He_D: float, B_D:float) -> float:
    """Returns the force coefficient for a simple 
    rectangular prismatic building

    Args:
        He_D (float): building effective height divided by building depth
        B_D (float): building breadth to depth ratio

    Raises:
        ValueError: height to depth ratio exeeds limit of 12

    Returns:
        float: force coefficient for a simple rectangular prismatic building
    """
    
    # He_D = min(12, He_D)
    if He_D > 12:
        err_msg = f'H_e / D is {He_D:4.2f} - greater than 12'
        raise ValueError(err_msg)
    else:
        C_f = 1.1 + 0.055 * He_D / (
            exp(
                abs(
                    log(0.6 * B_D * (1 - 0.011 * He_D))
                    )**(1.7 - 0.0013 * (He_D**2))
                )
            )
        return C_f


table_B3_1_dict = {
    3: {'A': 2.3, 'B': 1.4, 'C':1.2, 'D': None},
    5: {'A': 2.9, 'B': 1.8, 'C':1.4, 'D': 1.2},
    10: {'A': 3.4, 'B': 2.1, 'C':1.7, 'D': 1.2},
    }


def get_table_4_2(phi):
    """Net pressure coefficient Cp for open frameworks (Table 4-2)."""
    phis = (0.01, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9, 1.0)
    C_ps = (2.0, 1.8, 1.7, 1.6, 1.5, 1.5, 1.6, 2.0)
    return list_interp(phi, phis, C_ps, extend=True)


# Table 4-1: net pressure coefficients for an enclosed building envelope
# without dominant openings.  Keyed by (condition, zone) -> Cp.
table_4_1_dict = {
    'wall': {
        ('negative_edge', 'A'):      -1.4,
        ('negative_elsewhere', 'B'): -1.0,
        ('positive', 'A'):           1.1,
        ('positive', 'B'):           1.1,
    },
    'flat_roof': {
        ('negative_corner', 'C'):    -2.2,
        ('negative_edge', 'D'):      -1.6,
        ('negative_elsewhere', 'E'): -1.0,
        ('positive', 'C'):           0.3,
        ('positive', 'D'):           0.3,
        ('positive', 'E'):           0.3,
    },
    'pitch_roof': {
        ('negative_corner', 'C'):    -1.4,
        ('negative_edge', 'D'):      -1.4,
        ('negative_elsewhere', 'E'): -1.0,
        ('positive', 'C'):           1.1,
        ('positive', 'D'):           1.1,
        ('positive', 'E'):           1.1,
    },
}


def get_table_4_1(case: str, condition: str, zone: str) -> float:
    """Net pressure coefficient Cp for an enclosed building envelope
    without dominant openings (Table 4-1).

    Args:
        case (str): 'wall', 'flat roof', or 'pitch roof' (case-insensitive).
        condition (str): 'negative edge', 'negative corner',
            'negative elsewhere', or 'positive' (case-insensitive).
        zone (str): zone letter 'A', 'B', 'C', 'D', or 'E' (case-insensitive).

    Returns:
        float: net pressure coefficient Cp.

    Notes:
        - Positive pressure applies to every zone of the case (Table 4-1);
          each zone is listed explicitly in ``table_4_1_dict``.
        - Linear interpolation is allowed for roof pitch between 30 and
          60 degrees (Code Table 4-1 note (e)); this function returns the
          two bounding cases ('flat roof' < 30 deg, 'pitch roof' > 60 deg)
          only.

    Raises:
        KeyError: if ``case``, ``condition``, or ``zone`` is not recognised.
    """
    case = case.strip().lower().replace(' ', '_')
    condition = condition.strip().lower().replace(' ', '_')
    zone = zone.strip().upper()
    return table_4_1_dict[case][(condition, zone)]


S_theta_dict = {
    'N': 0.82,
    'NE': 0.84,
    'E': 0.85,
    'SE': 0.85,
    'S': 0.85,
    'SW': 0.84,
    'W': 0.82,
    'NW': 0.80,    
}


def get_S_theta(wind_direction='S'):
    """Wind directionality factor S_theta
    HK Wind 2019 Appendix A1, Table A1-1
    It is used in Equation3-1
    Note that S_theta is 1.0 for circular buildings"""
    S_theta = S_theta_dict.get(wind_direction, 0.85)
    return S_theta


def get_S_r(return_period=10):
    """Wind reference pressure adjustment factor S_r
    for calculation of acelerations
    HK Wind 2019 Appendix A1.2, Table A1-2
    """
    S_r_dict = {1: 0.25, 10:0.55}
    S_r = S_r_dict.get(return_period, None)
    return S_r


def get_Ss(L05p: float, zone: str = None) -> float:
    """Returns Size factor Ss from Appendix C1, eq C1-1(a,b,c)
    
    Args:
        L05p: half perimeter length of the equivalent loaded area
        zone: 'Edge' or 'Corner' or None
    
    Returns:
        float: Size factor Ss based on the loaded area
    
    """
    Ss = exp(0.17 - 0.07 * L05p**0.32)
    if zone == 'Edge':
        return max(Ss, 1.3 - log(L05p) / 9.0)
    elif zone == 'Corner':
        return max(Ss, 1.5 - log(L05p) / 5.4)
    else:
        return Ss


def get_S_q(H: float, B:float):
    """Approximate "Size and Dynamic" Factor for buildings of 50m or less
    HK Wind 2019 Equation 5-3
    
    Args:
        H: building height
        B: building width
    
    """
    L05p = H / 1.5 + 2 * B
    S_q = 1.1 * get_Ss(L05p)
    if H > 50:
        err_msg = f'Building is too tall for approximate method ({H:4.2f} > 50m)'
        raise ValueError(err_msg)
    else:
        return S_q


def get_S_qh(H: float, B:float, N_x, ksi_x=0.02):
    """'Size and Dynamic' Factor for gust forces at top of buildings
    HK Wind 2019 Equation 5-1
    
    Args:
        H: building height
        B: building width
    
    """
    if H is None or B is None:
        print(f'H = {H}; B = {B}')
    L05p = B
    S_s = get_Ss(L05p)
    S_qh = 0.5 + ((S_s - 0.5)**2.0 + 0.25 / (B**0.5 * H * N_x**2 * ksi_x)) **0.5
    return S_qh


def get_S_qz(Z: float, H: float, B:float,
             N_x:float=None, ksi_x:float=0.02, debug:bool=False)->float:
    """'Size and Dynamic' Factor for gust forces at top of buildings
    HK Wind 2019 Equation 5-1
    
    Args:
        Z: height of measurement
        H: building height
        B: building width
        N_x: fundamental frequency for mode mainly aligned 
                with the along-wind direction
        ksi_x: ratio of critical damping in the relevant 
                direction
        debug: 
    
    Raises:
        ValueError: H exceeds height limit of 50m for simple method 
            so N_x must be specified

    Returns:
        float: 'Size and Dynamic' Factor for gust forces at 
                top of buildings (Eq 5-1)
    """
    if N_x is None and H <= 50:
        if debug:
            L_05p = {H/1.5 + 2*B}
            print(f'get_S_qz: N_x is None and H = {H:4.2f} (<= 50m) and L_0.5p = {L_05p}')
        S_qh =  get_S_q(H, B)  # independent of Z
        return S_qh
    elif N_x is None and H > 50:
        if debug:
            print(f'get_S_qz: N_x is None and H = {H:4.2f} (> 50m)')
        err_msg = f'H is {H} (> 50) but N_x is not specified'
        raise ValueError(err_msg)
    else:
        S_qh = get_S_qh(H, B, N_x, ksi_x)
        S_qz = S_qh - 1.2 * (S_qh - (10/H)**0.14) * (1 - Z / H)
        if debug:
            print(f'get_S_qz: Z = {Z:4.4g}; N_x = {N_x:6.6g} (not None); ksi_x = {ksi_x:6.6g} ' + \
            f'\n\tand H = {H:4.2f} (> 50m) and L_0.5p = {B}, s_qz = {S_qz:8.6g}')
        return S_qz


def get_W_z(Q_z: float, C_f:float, S_qz: float, B: float):
    """Returns wind load per unit height at a height Z

    Args:
        Q_z (float): wind reference pressure (adjusted for sheltering, 
            topography & direction - see clause 3.1) 
        C_f (float): force coefficient in accordance with clause 4.2
        S_qz (float): size and dynamic factor from clause 5.5
        B (float): breadth of building

    Returns:
        float: wind load per unit height at height Z
    """
    W_z = Q_z * C_f * S_qz * B
    W_z - get_q_z()
    return W_z


def get_W_z_func(H_e: float, 
                 N_x: float = None, ksi_x: float = 0.01, 
                 S_t: float = 1.0, S_theta: float = 0.85,
                 as_dict:bool = False, debug:bool=False):
    """_summary_

    Args:
        H_e (float): Building effective height
        N_x (float): _description_
        ksi_x (float): _description_
        S_t (float): _description_
        S_theta (float): _description_
        as_dict (bool):  False, 
        debug (bool): False

    Returns:
        function[float]: a function that returns the pressure 
                when provided elevation Z, breadth B and depth D
    """
    def W_z_func(Ze: float, B: float, D: float) -> float:
        """Returns along-wind load per unit height when
        provided with elevation, breadth and depth

        Args:
            Ze (float): effective elevation
            B (float): building breadth at height Ze
            D (float): building width at height Ze

        Returns:
            float: along-wind load per unit height
        """
        Q_z = get_Q_oz(Ze) * S_t * S_theta
        C_f = get_C_f(H_e/D, B/D)
        if debug:
            print(f'*** Ze = {Ze:4.4g}  N_x = {N_x:4.4g}  ksi_x = {ksi_x:4.4g}***')
        S_qz = get_S_qz(Ze, H_e, B, N_x, ksi_x)
        W_z = Q_z * C_f * S_qz * B        
        return W_z
    
    def W_z_dict_func(Ze: float, B: float, D: float) -> float:
        """Returns along-wind load per unit height when
        provided with elevation, breadth and depth

        Args:
            Ze (float): effective elevation
            B (float): building breadth at height Ze
            D (float): building width at height Ze

        Returns:
            float: along-wind load per unit height
        """
        Q_z = get_Q_oz(Ze) * S_t * S_theta
        C_f = get_C_f(H_e/D, B/D)
        if debug:
            print(f'*** Ze = {Ze:4.4g}  N_x = {N_x:4.4g}  ksi_x = {ksi_x:4.4g}***')
        S_qz = get_S_qz(Ze, H_e, B, N_x, ksi_x)
        W_z = Q_z * C_f * S_qz * B
        W_z_dict_func = {'W_z': W_z, 'Ze': Ze, 'B': B, 
                         'Q_z': Q_z, 'C_f': C_f, 'S_qz': S_qz,
                         'N_x': N_x, 'ksi_x': ksi_x, 
                         'S_t': S_t, 'S_theta': S_theta,}
        
        return W_z_dict_func
    return W_z_dict_func if as_dict else W_z_func


def typical_damping_ratio(aspect_ratio:float, 
                        structural_type:str='RC',
                        analysis:str = 'loads') -> float:
    """Damping ratio for typical RC and steel buildings (tables C2-1 & C2-2)

    Args:
        aspect_ratio (float): building aspect ratio equal to the 
                total tower height divided by the dimension of the 
                tower in the direction of vibration.
        structural_type (str, optional): Structural type (RC' or 'steel'). 
                Defaults to 'RC'.
        analysis (str, optional): 'loads' for structural response damping',
                'accel' for acceleration damping. Defaults to 'loads'.

    Returns:
        float: Damping ratio for typical RC and steel buildings 
                (tables C2-1 & C2-2)
    """
    d_dict = {
        (True, True):[0.020,0.016,0.013,0.011,0.010], 
        (True, False):[0.030,0.024,0.020,0.017,0.015], 
        (False, True):[0.010,0.008,0.007,0.006,0.005], 
        (False, False):[0.015,0.012,0.010,0.009,0.008], 
        }
    aspect_ratios = [4, 5, 6, 7, 8] 
    key = (structural_type == 'RC', analysis != 'loads')
    return list_interp(aspect_ratio, aspect_ratios, d_dict.get(key), extend=True)


def approx_fundamental_frequency(H_e:float) -> float:
    """Approximate building period N from section 1.2 (Symbols)

    Args:
        H_e (float): Height of building 

    Raises:
        ValueError: Specified height H_e is > 100m so dyanmic analysis is required

    Returns:
        float: Approximate building period N
    """
    if H_e > 100:
        err_msg = f'H is {H_e} (> 100)'
        raise ValueError(err_msg)
    else:
        return 46 / H_e


def freestanding_wall_pressure(phi: float, l_h: float, lc_h: float = 0.0) -> dict:
    """Net pressure coefficients Cp for free-standing walls (Table B3-1).

    Args:
        phi (float): solidity ratio (0.8 to 1.0; use Table 4-2 if < 0.8).
        l_h (float): wall length / height ratio.
        lc_h (float, optional): return-corner length / height ratio
            (0 = no return corners, >= 1 = full return corners). Defaults to 0.

    Returns:
        dict: zone pressure coefficients {'A', 'B', 'C', 'D'}. Zone 'D' is
        None for walls with l/h <= 3 (no zone D).

    Raises:
        ValueError: if phi < 0.8 (use Table 4-2 for open frameworks).
    """
    if phi < 0.8:
        raise ValueError('phi < 0.8: use Table 4-2 for open frameworks')

    lh_keys = (3, 5, 10)
    full_return = {'A': 2.1, 'B': 1.8, 'C': 1.4, 'D': 1.2}

    no_return = {z: list_interp(l_h, list(lh_keys),
                                [table_B3_1_dict[k][z] for k in lh_keys],
                                extend=True) for z in 'ABCD'}

    f_ret = max(0.0, min(1.0, lc_h))  # interpolate return-corner length
    cp1 = {z: (no_return[z] + (full_return[z] - no_return[z]) * f_ret
               if no_return[z] is not None else None) for z in 'ABCD'}

    if phi >= 1.0:
        return cp1
    cp08 = {'A': 1.2, 'B': 1.2, 'C': 1.2, 'D': 1.2}
    f = (phi - 0.8) / 0.2
    return {z: (cp08[z] + (cp1[z] - cp08[z]) * f
                if cp1[z] is not None else None) for z in 'ABCD'}


def cantilever_deflection(w, L, E, Ix):
    return 0.125 * w * L**4 / E / Ix


def cantilever_forces(p, B, L):
    V = p * B * L
    return V, 0.5 * V * L


# =============================
# ====== Plot Functions =======
# =============================


def plot_Ss(Lmax = 3, L05p = 6.1, n=30, debug=False):
    """
    Appendix C1, eq C1-1
    """
    
    L_range = [10**(i/n*Lmax) for i in range(n+1)]
    if debug:
        print(f'Q_oz(6.1) is {get_Q_oz(6.1):.3f}')
        print(f'S_s(6.1) is {get_Ss(6.1):.3f}')
        print(get_Ss(7.5), exp(0.17 - 0.07 * 7.5**0.32), 
            get_Ss(12), exp(0.17 - 0.07 * 12**0.32))
    S_range = [get_Ss(L,'Corner') for L in L_range]
    plt.semilogx(L_range, S_range)
    S_range = [get_Ss(L, 'Edge') for L in L_range]
    plt.semilogx(L_range, S_range)
    S_range = [get_Ss(L) for L in L_range]
    plt.semilogx(L_range, S_range)
    plt.legend(['Corner', 'Edge', 'Other'])
    plt.grid(True)
    plt.scatter(L05p, get_Ss(L05p), color = '#88c999')


def print_table_F1():
    """

    """
    print('Table F1 : Dynamic Magnification Factor G for ζ = 2.0%')
    [[round(get_G(h, b, zeta = 0.015),3) 
    for b in [20, 30, 40]] 
    for h in [200, 180, 160, 140, 120, 100]]

heights = [0, 5, 10, 20, 30, 
            50, 75, 100, 150, 
            200, 250, 300, 400, 
            500, 600]


def print_table_F2():
    """

    """
    print('Table F2 : Dynamic Magnification Factor G for ζ = 1.5%')
    [[round(get_G(h, b, zeta = 0.02),3) 
      for b in [20, 30, 40]] 
      for h in [200, 180, 160, 140, 120, 100]]


def print_table_F3(heights=heights):
    """

    """
    print('Table F3: Design hourly-mean wind velocity')
    v_z_mean = [get_v_z_mean(z) for z in heights]
    [print(f'{z:4.0f} m | {v:4.2f} m/s') for z, v in zip(heights, v_z_mean)]


def print_table_3sec_wind(heights=heights):
    """

    """
    print('Design 3 second gust wind speed at height z')
    v_z = [get_v_z(z) for z in heights]
    [print(f'{z:4.0f} m | {v:4.2f} m/s') for z, v in zip(heights, v_z)]


def print_table_1(heights=heights):
    """

    """
    v_z = [get_v_z(z) for z in heights]
    print('Table 1: Design Wind Pressure (kPa)')
    [print(f'{z:4.0f} m | {v:4.2f} m/s ' + \
        f'| {get_q_z(v):4.2f} kPa') 
        for z, v in zip(heights, v_z)]


# ===============================
# ====== Summary Functions ======
# ===============================


def full_wind_calcs(
        storey_heights, 
        storey_widths_x, storey_widths_y,
        N_x = None, N_y = None, ksi_x = 0.005, ksi_y = 0.005,  
        ascending = True,
        rigorous=True, print_data=False, debug=False,
        **kwargs):
    """Returns full set of calculations for HK Wind Code 2019"""
    S_e = 0.85
    if print_data:
        print('Initial settings')
        fx = '4.4g' if isinstance(N_x, float) else ''
        fy = '4.4g' if isinstance(N_y, float) else ''
        if debug:
            print(f'-- N_x = {N_x:{fx}}; ksi_x = {ksi_x:4.4g}')
            print(f'-- N_y = {N_y:{fy}}; ksi_y = {ksi_y:4.4g}')
    
    # This could be speeded up significantly by using numpy
    # print(f'Bottom storey height is {elevs[0]} m, 
    # top storey elevation is {elevs[-1]} m, sum of storeys is {sum(storey_heights)} m')
    dd = {}
    
    if not ascending:
        storey_heights = list(storey_heights)[::-1]
    if storey_heights[0] != 0:
        storey_heights = [0] + storey_heights
    
    dd['storey_heights'] = storey_heights

    if not isinstance(storey_widths_x, (list, tuple)):
        bx_ave = float(storey_widths_x)
        dd['storey_widths_x'] = [bx_ave for _ in storey_heights]
    else:
        if len(storey_widths_x) != len(storey_heights):
            sw_x, sh = len(storey_widths_x), len(storey_heights)
            raise ValueError(f'Number of storey heights, {sh}, ' + \
                f'does not match number of storey widths (x), {sw_x}')
        bx_ave = sum(storey_widths_x) / len(storey_widths_x)
        dd['storey_widths_x'] = storey_widths_x

    if not isinstance(storey_widths_y, (list, tuple)):
        by_ave = float(storey_widths_y)
        dd['storey_widths_y'] = [by_ave for n in storey_heights]
    else:
        if len(storey_widths_y) != len(storey_heights):
            raise ValueError('Number of storey heights ' + \
                'does not match number of storey widths (y)')
        by_ave = sum(storey_widths_y) / len(storey_widths_y)
        dd['storey_widths_y'] = storey_widths_y
    
    if debug:
        print(f'Storey Widths (x) [dim {len(dd["storey_widths_x"])}]: ', 
                    dd['storey_widths_x'])
    if debug:
        print(f'Storey Widths (y) [dim {len(dd["storey_widths_y"])}]: ', 
                    dd['storey_widths_y'])
    if debug:
        print(f'bx_ave = {bx_ave}; by_ave = {by_ave} ')
        print(f'Storey Heights [dim {len(dd["storey_heights"])}]: ', 
                    dd['storey_heights'])
    dd['storey_elevs'] = gen_elevs(dd['storey_heights'])
    if debug:
        print(f'Storey Elevations [dim {len(dd["storey_elevs"])}]: ', 
                    dd['storey_elevs'])
    dd['trib_heights'] = gen_trib_heights(dd['storey_heights'])
    if debug:
        print(f'Tributary heights [dim {len(dd["trib_heights"])}]: ', 
                    dd['trib_heights'])
    dd['storey_mean_vels'] = [get_v_z(z) for z in dd['storey_elevs']] 
    dd['storey_pressures'] = [get_q_z(v) for v in dd['storey_mean_vels']]
    
    H_e = max(dd['storey_elevs'])
    if debug:
        print('Building Height: ', H_e)
    
    dd['C_fx'] = [get_C_f(H_e / D, B / D) 
                    for B, D in zip(dd['storey_widths_x'], 
                                dd['storey_widths_y'])]
    dd['C_fy'] = [get_C_f(H_e / D, B / D) 
                    for B, D in zip(dd['storey_widths_y'], 
                                dd['storey_widths_x'])]
    
    dd['Q_oz'] = [get_Q_oz(Z) 
                    for Z in dd['storey_elevs']] 
     
    dd['Q_z'] = [S_e * Q_z 
                                for Q_z in dd['Q_oz']]    
    
    damp_x, damp_y = ksi_x, ksi_y
    if N_x is None and H_e <= 100:
        N_x = approx_fundamental_frequency(H_e)
        damp_x = typical_damping_ratio(H_e / bx_ave, 'RC', 'loads')
    if N_y is None and H_e <= 100:
        N_y = approx_fundamental_frequency(H_e)
        damp_y = typical_damping_ratio(H_e / by_ave, 'RC', 'loads')

    if print_data: 
        print(f'aspect = {H_e / bx_ave:4.4g}; N_x = {N_x:4.4g}; damp_x = {damp_x:4.4g}')
        print(f'aspect = {H_e / by_ave:4.4g}; N_y = {N_y:4.4g}; damp_y = {damp_y:4.4g}')
    
    S_t = kwargs.get('S_t', 1.0)
    S_theta = kwargs.get('S_theta', 0.85)
    if print_data:
        dd['S_qxz'] = [get_S_qz(Z, H_e, B,
                N_x, ksi_x, debug=debug) for Z, B in zip(dd['storey_elevs'], dd['storey_widths_x'])]
        dd['S_qyz'] = [get_S_qz(Z, H_e, B,
                N_y, ksi_y, debug=debug) for Z, B in zip(dd['storey_elevs'], dd['storey_widths_y'])]
        
    get_W_z_x = get_W_z_func(H_e, N_x, damp_x, S_t, S_theta)
    get_W_z_y = get_W_z_func(H_e, N_y, damp_y, S_t, S_theta)
    
    dd['W_x'] = [get_W_z_x(z, b, d) 
                 for z, b, d in zip(dd['storey_elevs'],
                                    dd['storey_widths_y'], 
                                    dd['storey_widths_x'],
                                    )]
    dd['W_y'] = [get_W_z_y(z, b, d) 
                 for z, b, d in zip(dd['storey_elevs'],
                                    dd['storey_widths_x'], 
                                    dd['storey_widths_y'],
                                                        )]
    if debug:
        print('W_x (dim), W_y (dim): ', len(dd['W_x']), ',', len(dd['W_y']))
    
    dd['storey_forces_x'] = [F * h 
                        for F, h in zip(dd['W_x'], dd['trib_heights'])]
    dd['storey_forces_y'] = [F * h 
                        for F, h in zip(dd['W_y'], dd['trib_heights'])]
    
    dd['storey_shears_x'] = get_storey_shear(dd['storey_forces_x'])
    dd['storey_shears_y'] = get_storey_shear(dd['storey_forces_y'])
    dd['storey_moments_yy'] = get_storey_moment(dd['storey_shears_x'], 
                                                dd['storey_heights'])
    dd['storey_moments_xx'] = get_storey_moment(dd['storey_shears_y'], 
                                                dd['storey_heights'])
    
    ## END OF CALCULATIONS #######

    # Check that maximum shear matches sum of forces and 
    # maximum moment is the sum of the products of the storey forces and storey elevations
    sum_of_storey_forces_x = sum(dd['storey_forces_x'])
    sum_of_storey_moments_yy = sum(F * h 
        for F, h in zip(dd['storey_forces_x'], 
                        dd['storey_elevs']))
    sum_of_storey_forces_y = sum(dd['storey_forces_y'])
    sum_of_storey_moments_xx = sum(F * h 
        for F, h in zip(dd['storey_forces_y'], 
                        dd['storey_elevs']))

    dd['Vx_base'] = sum_of_storey_forces_x
    dd['Vy_base'] = sum_of_storey_forces_y
    dd['Myy_base'] = sum_of_storey_moments_yy
    dd['Mxx_base'] = sum_of_storey_moments_xx
        
    
    err_list = []
    max_storey_shear_x = dd['storey_shears_x'][0]
    max_storey_shear_y = dd['storey_shears_y'][0]
    max_storey_moment_xx = dd['storey_moments_xx'][0]
    max_storey_moment_yy = dd['storey_moments_yy'][0]
    
    if round(sum_of_storey_forces_x) != round(max_storey_shear_x):
        err_msg = 'Error in sum of x-dir forces ' + \
                f'({sum_of_storey_forces_x:6,.6g} vs {max_storey_shear_x:6,.6g})'
        err_list.append(err_msg)
    if round(sum_of_storey_moments_yy) != round(max_storey_moment_yy):
        err_msg = 'Error in sum of yy-dir moments ' + \
                f'({sum_of_storey_moments_yy:6,.6g} vs {max_storey_moment_yy:6,.6g})'
        err_list.append(err_msg)
    if round(sum_of_storey_forces_y) != round(max_storey_shear_y):
        err_msg = 'Error in sum of y-dir forces ' + \
                f'({sum_of_storey_forces_y:6,.6g} vs {max_storey_shear_y:6,.6g})'
        err_list.append(err_msg)
    if round(sum_of_storey_moments_xx) != round(max_storey_moment_xx):
        err_msg = 'Error in sum of xx-dir moments ' + \
                f'({sum_of_storey_moments_xx:6,.6g} vs {max_storey_moment_xx:6,.6g})'
        err_list.append(err_msg)
    if len(err_list) > 0 and rigorous:
        err_msg = '\n'.join(err_list)
        raise ValueError(err_msg)
    elif len(err_list) > 0: 
        print('WARNING:', err_msg)
    else:
        if debug:
            print('No errors / warnings')

    ## PRINT RESULTS if print_data is True #######
    if print_data:
        print('print_data is True\n')
        S_qx = get_S_qz(H_e, H_e, by_ave, N_x, ksi_x, debug=debug)
        S_qy = get_S_qz(H_e, H_e, bx_ave, N_x, ksi_x, debug=debug)
        print(f'S_e : {S_e:4.3f}     S_qx : {S_qx:4.3f}     S_qy : {S_qy:4.3f}\n')

        # print(f'C_hx : {C_hx:4.3f}     C_hy : {C_hy:4.3f}\n')

        # print('Height | Storey Widths (x&y)| Trib. Ht. (x & y)| Mean Vels | Pressure')
        # [print(f'{a:4.2f} m | {bx:4.2f} m | {by:4.2f} m | ' + \
        #                f'{hx:4.2f} m | {hy:4.2f} m    | {v:4.2f} m/s | {p:4.2f} kPa') \
        #             for a, bx, by, hx, hy, v, p 
        #             in zip(dd['storey_heights'], 
        #                            dd['storey_widths_x'], dd['storey_widths_y'], 
        #                            dd['trib_heights_x'], dd['trib_heights_y'], 
        #                            dd['storey_mean_vels'], dd['storey_pressures'])]
        # print()
        
        keys = ['storey_heights', 'storey_elevs', 'Q_oz', 'Q_z', 
                        'storey_widths_x', 'storey_widths_y', 'C_fx', 'C_fy']
        if debug:
            print([len(dd[k]) for k in keys])
        print('Height | Elev    | Q_oz     | Q_z       ' + \
                '| B_x       | B_y          | C_fx         | C_fy     ')
        [print(f'{h:4.2f} m | {z:5.5g} m | {Q_oz:4.2f} kPa| {Q_z:4.2f} kPa | ' + \
                     f'{bx:5.5g} m | {by:5.5g} m    | {C_fx:5.5f} kPa | {C_fy:5.5f} kPa') \
                    for h, z, Q_oz, Q_z, bx, by, C_fx, C_fy 
                    in zip(dd['storey_heights'], dd['storey_elevs'],
                                 dd['Q_oz'], dd['Q_z'],
                                 dd['storey_widths_x'], dd['storey_widths_y'], 
                                 dd['C_fx'], dd['C_fy'])]
        print()

        # print('Height | Trib. Ht. (x & y)| Mean Vels | Press.     | Forces (x & y)')
        # [print(f'{a:4.2f} m | {ax:4.2f} m | {ay:4.2f} m    | ' + \
        #                f'{v:4.2f} m/s | {p:4.2f} kPa | {x:7,.0f} kN | {y:7,.0f} kN') \
        #         for a, ax, ay, v, p, x, y 
        #         in zip(dd['storey_heights'], dd['trib_heights'], 
        #                        dd['trib_heights'], dd['storey_mean_vels'], 
        #                        dd['storey_pressures'], dd['storey_forces_x'], 
        #                        dd['storey_forces_y'])]
        # print()

        keys = ['storey_elevs', 'trib_heights', 'S_qxz', 'S_qyz',
                        'W_x', 'W_y', 'storey_forces_x', 'storey_forces_y']
        if debug:
            print([len(dd[k]) for k in keys])
        print('Elev.        | Trib. H| S_q,z            | Wx        | Wy        | ' + \
            'Forces (x & y)')
        [print(f'{z:6.2f} m | {h:4.2f} m    | ' + \
                     f'{sqxz:6.4g} | {sqyz:6.4g} | ' + \
                     f'{wx:6.6g} kN/m | {wy:6.6g} kN/m | ' + \
                     f'{x:7,.0f} kN | {y:7,.0f} kN') \
                for z, h, sqxz, sqyz, wx, wy, x, y 
                in zip(dd['storey_elevs'], dd['trib_heights'], 
                             dd['S_qxz'], dd['S_qyz'],
                             dd['W_x'], dd['W_y'], 
                             dd['storey_forces_x'], dd['storey_forces_y'])]
        print()

        print('Vx,base = {:7,.0f} kN,\t\t Myy,base = {:7,.0f} kNm'.format(
            sum_of_storey_forces_x, sum_of_storey_moments_yy))
        print('Vy,base = {:7,.0f} kN,\t\t Mxx,base = {:7,.0f} kNm'.format(
            sum_of_storey_forces_y, sum_of_storey_moments_xx))
        print()

        print('Height    | Elev         | F_x                | ' + \
            'V_x                | M_yy                | ' + \
                'V_y                | M_xx     ')
        [print(f'{h:5.5g} m | {z:6.6g} m | {Fx:7,.0f} kN | ' + \
                     f'{Vx:7,.0f} kN | {Myy:7,.0f} kNm | ' + \
                     f'{Vy:7,.0f} kN | {Mxx:7,.0f} kNm') \
                for h, z, Fx, Vx, Myy, Vy, Mxx 
                in zip(dd['storey_heights'], dd['storey_elevs'], 
                             dd['storey_forces_x'], dd['storey_shears_x'], 
                             dd['storey_moments_yy'], dd['storey_shears_y'], 
                             dd['storey_moments_xx'])]
        print()

    return dd


def main():
    pass


if __name__ == "__main__":
    doctest.testmod()
    print('== Doctest Complete ==')
    main()

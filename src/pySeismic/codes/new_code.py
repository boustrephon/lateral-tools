"""

"""

from math import pi, exp, log, log10, isclose #, tan, sin, cos, atan2  # noqa: F401
from itertools import accumulate
from csv import writer 
import matplotlib.pyplot as plt
import doctest
from sys import path as sys_path
from os.path import abspath, join as path_join, dirname
# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

# from pySeismic.codes.new_code import T_list_gen, get_S_e, get_S_p
from pySeismic.utilities import shear_coords



def interp(x, x1, x2, y1, y2):
    """
    """
    if y1 is None or y2 is None:
        return None
    elif x1 == x2:
        return 0.5 * (y1 + y2)
    else:
        return y1 + (x - x1) / (x2 - x1) * (y2 - y1)


def list_interp(val, list_1, list_2):
    """
    Returns an interpolated value based on a list of reference values and a list of target values
    val is the lookup value
    list_1 is the lookup list
    list_2 is the result list with values corresponding to those on list_1
    """
    i_list = [i for i, n in enumerate(list_1) if n == val]
    if i_list:
        return list_2[i_list[0]]
    elif val <= list_1[0]:
        return list_2[0]
    elif val >= list_1[-1]:
        return list_2[-1]
    else:
        index_list = [i for i, (x1, x2) in enumerate(zip(list_1[:-1], list_1[1:])) if val >= x1 and val <= x2]
        if len(index_list) == 1:
            j = index_list[0]
            return interp(val, list_1[j], list_1[j+1], list_2[j], list_2[j+1])
        else:
            return None


def exceedance2return(exceedance_probability, period=50):
    """Returns the return period for ground motion with the given probability of
    being exceeded over the period considered. For example, the return period 
    for a ground motion with a 10% probability of being exceeded in 50 years
    is around 475 years. 
    >>> isclose(exceedance2return(exceedance_probability=0.1, period=50), 474.56107905149526, rel_tol=1e-9)
    True
    >>> isclose(exceedance2return(exceedance_probability=0.2, period=50), 224.07100588622754, rel_tol=1e-9)
    True
    """
    return - period / log(1 - exceedance_probability)


def return2exceedance(return_period, period=50):
    """Returns the probability of a ground motion being exceeded over a period  
    when given a ground motion return period. For example, the probability of a  
    ground motion with return period of 475 years being exceeded over a time  
    period of 50 years is around 10%. 
    >>> isclose(return2exceedance(return_period=475, period=50), 0.09991237374774076, rel_tol=1e-9)
    True
    >>> isclose(return2exceedance(return_period=225, period=50), 0.19926259708319194, rel_tol=1e-9)
    True
    """
    return 1 - exp(-period/return_period)


def PGA2S475(a_g_475):
    """"Returns the (S_1, S_s) tuple for a reference PGA.
    Reference: 15WCEE, 'Deriving S_s and S_1 parameters from PGA maps'
    """
    def S_s_calc(PGA): return 0.3386 * PGA + 2.1696
    def S_1_calc(PGA): return 0.5776 * PGA + 0.5967
    return {'S_s':S_s_calc(a_g_475), 'S_1':S_1_calc(a_g_475)}


def S_calc(S_s_475, S_1_475, T):
    """"Returns the (S_1, S_s) tuple for a given return period based on a 
    reference PGA and reference return period.
    Reference: 15WCEE, 'Deriving S_s and S_1 parameters from PGA maps'
    """
    def S_calc(S_475, k, T): return S_475 * (T / 475)**k
    S_s_475_S_keys = (0.20, 0.55, 0.90, 1.20, 1.50)
    S_s_475_k_vals = (0.46, 0.42, 0.39, 0.35, 0.30)
    S_1_475_S_keys = (0.01, 0.02, 0.05, 0.10, 0.20)
    S_1_475_k_vals = (0.40, 0.40, 0.40, 0.40, 0.40)
    k_s = list_interp(S_s_475, S_s_475_S_keys, S_s_475_k_vals)
    S_s = S_calc(S_s_475, k_s, T)
    k_1 = list_interp(S_1_475, S_1_475_S_keys, S_1_475_k_vals)
    S_1 = S_calc(S_1_475, k_1, T)
    return {'S_s':S_s, 'S_1':S_1, 'k_s':k_s, 'k_1':k_1}


# Note that these classes have the same names as EC8, but reversed importance
importance_classes = ['I', 'II', 'III', 'IV']


importance_class_dict_EC8 = {
  'I': 'Civil protection - buildings and structures essential for civil protection.',
  'II':	'High consequences in case of failure - buildings and structures of high importance to the society.',
  'III':	'Ordinary - structures of ordinary importance and not classified as I, II and IV.',
  'IV':	'Minor or temporary structures - structures of minor importance for public safety and temporary structures.'
}


importance_factor_dict_EC8 = {
  'I': 1.7,
  'II': 1.4,
  'III': 1.0,
  'IV': None
  }


importance_class_dict = {
  'A': '''Civil protection - buildings and structures essential for civil protection:
  Critical infrastructure, hazardous facilities, life rescue facilities, ATC towers, communication & hazardous''',
  'B':	'''High consequences in case of failure - buildings and structures of high importance to the society:
  Important infrastructure, dangerous good storage, public health care, transportation, communication, education, public assembly''',
  'C':	'Ordinary - structures of ordinary importance and not classified as I, II and IV.',
  'D':	'Minor or temporary structures - structures of minor importance for public safety and temporary structures.'
}


importance_factor_dict = {
  'A': 1.7,
  'B': 1.4,
  'C': 1.0,
  'D': 0.0
  }


def get_importance_factor(importance_class:str, 
                          importance_factor_dict:dict=importance_factor_dict,
                          ) -> float:
    """Returns the Importance factor I_E for occupancy classes A-C (D is not covered)
    Default values are provided, but can be replaced with EC8 NAD values
    by specifying the importance_factor_dict

    Args:
        importance_class (str): _description_
        importance_factor_dict (dict, optional): _description_. Defaults to importance_factor_dict.

    Returns:
        float: _description_
    """
    if importance_class in ('A', 'B', 'C', 'D'):
        if importance_class == 'D':
            print('No seismic required')
        return importance_factor_dict.get(importance_class)
    else:
        err_msg = 'Importance class should be one of A, B, C or D'
        err_msg += f' (currently {importance_class})'
        raise ValueError(err_msg)


def get_v_s(N, weathering=True):
    """"""
    return 113*N**0.27 if weathering else 121*N**0.19

def get_v_s_30(descending_list_of_tuples_of_thickness_and_v_s):
    """"""
    depth = 0
    d_v = 0
    #print('   d_i, depth, v_s_i, d/v_s,  d_v')
    for d_i, v_s_i in descending_list_of_tuples_of_thickness_and_v_s:
        if depth + d_i >= 30:
            d_i = 30 - depth
            d_v += d_i / v_s_i
            #print(f' {d_i:4.1f},  {(depth+d_i):4.1f},  {v_s_i:5.1f},  {(d_i / v_s_i):5.3f},  {d_v:5.3f},  (end)')
            return 30 / d_v
        else:
            depth += d_i
            d_v += d_i / v_s_i
            #print(f' {d_i:4.1f},  {(depth):4.1f},  {v_s_i:.1f},  {(d_i / v_s_i):5.3f},  {d_v:5.3f}')
    err_msg = f'Soil data needs to cover 30m. Current data only covers {depth}m'
    raise ValueError(err_msg)
    return None


def get_soil_class(v_s_30):
    """Table 3.1.1-1"""
    if v_s_30 < 100:
        return 'S1'
    elif v_s_30 < 180:
        return 'D'
    elif v_s_30 < 360:
        return 'C'
    elif v_s_30 < 800:
        return 'B'
    else:
        return 'A'

horizontal_response_parameter_dict = {
    'A': {'S': 1.0, 'T_B': 0.12, 'T_C':0.32, 'T_D':2.5},
    'B': {'S': 1.6, 'T_B': 0.12, 'T_C':0.34, 'T_D':1.5},
    'C': {'S': 1.5, 'T_B': 0.22, 'T_C':0.50, 'T_D':1.5},
    'D': {'S': 1.5, 'T_B': 0.22, 'T_C':0.65, 'T_D':1.5}
}


# Clause 4.8
gamma_ov_dict = {
    'Regular Concrete':{'q':3, 'gamma_ov':1.5}, 
    'High-strength Concrete':{'q':1, 'gamma_ov':1.5}, 
    'Composite':{'q':4, 'gamma_ov':1.5}, 
    'Steel_General':{'q':4, 'gamma_ov':1.25}, 
    'SCBF_DV':{'q':2, 'gamma_ov':1.25}, 
    'SCBF_K':{'q':1.5, 'gamma_ov':1.25},
    }


def get_gamma_ov(material: str) -> dict:
    """Returns dictionary of beahaviour factor 'q' 
    and 'overstrength factor , 'gamma_ov'
    

    Args:
        material (str): Concrete is high strength if the grade is 
        > C60 for normal and prestressed concrete, 
        or > C50 for precast concrete used in dissipative members

    Returns:
        dict: dictionary of q and gamma_ov
    """
    mat_dict = gamma_ov_dict.get(material, {'q':1.5, 'gamma_ov':1.5})
    return mat_dict


def get_seismic_weight(D, L_P, L_R):
    """Returns seismic weight (Section 4.6.1)"""
    return D + L_P + 0.3 * L_R


def get_T_code(h):
    return h / 46


def get_T_rayleigh(
    storey_weights, 
    F_b=None, 
    storey_deflections = [], 
    delta_R = None, 
    g = 9.81,
    ):
    """Returns building period based on the Rayleigh method (section 4.3.1)
    If only storey_weights and delta_R are provided it returns an approximate value
      storey_weights - seismic weight of each floor - should be in same units as base shear
      F_b - base shear
      storey_deflections - deflection of each floor (should be consistent units with g)
      delta_R - roof displacement (if this is not defined, the maximum storey deflection is used)
      g - constant of gravitational acceleration (should be in consistent units)
    Note that these calculations assume storey weights not masses
    """
    if delta_R is None and len(storey_deflections) > 0:
        delta_R = max(storey_deflections)

    if len(storey_deflections) == 0 and delta_R is None:
        err_msg = f'Either delta_R ({delta_R}) or storey deflections must be defined'
        raise ValueError(err_msg)
    elif F_b is None and delta_R is None:
        err_msg = f'delta_R ({delta_R}) and / or F_b ({F_b}) and storey_deflections must be provided'
        raise ValueError(err_msg)
    elif len(storey_deflections) == 0 or F_b is None:
        G_e = sum(storey_weights) if isinstance(storey_weights,list) else storey_weights
        w_eq = 0.5 * G_e # in units of weight, i.e. force
        Gamma = 1.5  # dimensionless
        F_b = 0.85 * G_e  # in units of force
    elif len(storey_deflections) != len(storey_weights):
        err_msg = (
            f'Length of storey_deflections ({len(storey_deflections)})' 
            + f' must equal that of storey_weights ({len(storey_weights)})'
            )
        raise ValueError(err_msg)
    else:
        sum_wd = sum(w * d for w, d in zip(storey_weights, storey_deflections))
        sum_wd2 = sum(w * d**2 for w, d in zip(storey_weights, storey_deflections))
        w_eq =  sum_wd / delta_R
        Gamma = w_eq * delta_R**2 / sum_wd2
    F_eq = F_b / Gamma
    d_eq = delta_R / Gamma
    return 2 * pi * (w_eq * d_eq / F_eq / g)**0.5


def T_list_gen(soil_class, n_pt=20):
    """Returns a list of periods suitable for plotting response spectra
        n_pt is the number of points between each corner period
    """
    params_dict = horizontal_response_parameter_dict[soil_class]
    T_vals = [0] + [params_dict.get(p) for p in ['T_B', 'T_C', 'T_D']] + [10]
    lists = [[(T1 * float(n_pt - i) + T2 * float(i)) / n_pt for i in range(n_pt)] 
           for T1, T2 in zip(T_vals[:-1], T_vals[1:])]
    return sum(lists, []) + [10]

def get_a_g(a_gR, importance_factor):
    return a_gR * importance_factor


def get_S_e(
    T:float, 
    soil_class:str, 
    importance_factor:float=1.0, 
    S_T:float=1.0, 
    eta:float=1.0, 
    a_gR:float = 0.12, 
    horizontal_response_parameter_dict:dict=horizontal_response_parameter_dict,
    debug = False,
    ) -> float:
    """Returns Horizontal Elastic Spectral Acceleration, S_e
    Section 3.3.1 "Horizontal Elastic Response Spectrum"
    
    Args:
        T (float): _description_
        soil_class (str): _description_
        importance_factor (float, optional): _description_. Defaults to 1.0.
        S_T (float, optional): topographic amplification factor 
            defined in clause 3.6. Defaults to 1.0.
        eta (float, optional): modification factor for hysteretic damping. Defaults to 1.0.
        a_gR (float, optional): reference ground spectral acceleration. Defaults to 0.12.
        horizontal_response_parameter_dict (dict, optional): Dictionary containing 
                soil response parameters from Table 3.3.1-1. 
                Defaults to horizontal_response_parameter_dict.
        debug (bool, optional): flag for turning on debug printing

    Raises:
        ValueError: Building period is too long (> 10 secs)

    Returns:
        float: Returns Elastic Spectral Acceleration, S_e
    """
    if debug:
        print(T, soil_class,importance_factor,S_T,eta,a_gR)
    params_dict = horizontal_response_parameter_dict[soil_class]
    S, T_B, T_C, T_D = [params_dict.get(p) for p in ['S', 'T_B', 'T_C', 'T_D']]
    a_g = importance_factor * a_gR  # units are 'g'
    if debug:
        print('\n### get_S_e function ###')
        print('S,    T_B,  T_C,  T_D')
        print(f'{S:4.2f}, {T_B:4.2f}, {T_C:4.2f}, {T_D:4.2f}\n')
        print(f'T = {T:5.3f}; a_gR = {a_gR}; I_E = {importance_factor}')
    if T > 10:
        err_msg = f'Period T is too high (T = {T:.2f} > 10'
        raise ValueError(err_msg)
    elif T > T_D:
        S_e = a_g * S * S_T * eta * 2.5 * T_C * T_D / T**2
        if debug:
            print(f'T > T_D ({T_D}) : S_e = {S_e:5.3f}\n')
        return S_e
    elif T > T_C:
        S_e = a_g * S * S_T * eta * 2.5 * T_C / T
        if debug:
            print(f'T > T_C ({T_C}) : S_e = {S_e:5.3f}\n')
        return S_e
    elif T > T_B:
        S_e = a_g * S * S_T * eta * 2.5
        if debug:
            print(f'T > T_B ({T_B}) : S_e = {S_e:5.3f}\n')
        return S_e
    else:
        return a_g * S * S_T * (1 + T / T_B * (eta * 2.5 - 1))


def get_S_p(
        T:float, 
        soil_class:str, 
        importance_factor:float=1.0, 
        S_T:float=1.0, 
        eta:float=1.0, 
        a_gR:float = 0.12, 
        q:float=1.5, 
        gamma_ov:float=1.5, 
        horizontal_response_parameter_dict:dict=horizontal_response_parameter_dict, 
        debug:bool=False,
    ) -> float:
    """Returns Horizontal Inlastic Spectral Acceleration, S_p
    From Section 4.2.1 "Behaviour Factor"
    
    Args:
        T (float): _description_
        soil_class (str): _description_
        importance_factor (float, optional): _description_. Defaults to 1.0.
        S_T (float, optional): topographic amplification factor 
            defined in clause 3.6. Defaults to 1.0.
        eta (float, optional): modification factor for hysteretic damping. Defaults to 1.0.
        a_gR (float, optional): reference ground spectral acceleration. Defaults to 0.12.
        q (float, optional): behaviour factor from Table 5.3-1 for concrete and 
                Table 6.4-1 for steel. Defaults to 1.5, but note that 
                high strength concretes may require q to be 1.0.
        gamma_ov (float, optional): material overstrength (Table 4.8-1). Defaults to 1.5.
        horizontal_response_parameter_dict (dict, optional): Dictionary containing 
                soil response parameters from Table 3.3.1-1. 
                Defaults to horizontal_response_parameter_dict.
        debug (bool, optional): flag for turning on debug printing

    Returns:
        float: Returns Inelastic Spectral Acceleration, S_p
    """
    if debug:
        print(T, soil_class,importance_factor,S_T,eta,a_gR,q,gamma_ov)
    
    T_C = horizontal_response_parameter_dict[soil_class]['T_C']
    S_e = get_S_e(
        T, 
        soil_class, 
        importance_factor=importance_factor, 
        S_T=S_T, 
        eta=eta, 
        a_gR = a_gR, 
        horizontal_response_parameter_dict=horizontal_response_parameter_dict,
        debug = debug
      )
    scaling_factor = 1 / q if T > T_C else 1 / (gamma_ov + (q - gamma_ov) * T/T_C)
    S_p = S_e * scaling_factor
    if debug:
        print('\n### get_S_p function ###')
        print(f'T ({T:4.2f} secs) ' + \
            ('> ' if T > T_C else '<= ') + \
            f'T_C  ({T_C:4.2f} secs)')
        print(f'q = {q}; gamma_ov = {gamma_ov}')
        print(f'Scaling factor is {scaling_factor:4.3f}')
        print(f'S_e = {S_e:4.3f}')
        print(f'S_p = {S_p:4.3f}\n')
    return S_p


#"""I think it would be helpful to make clear that ground acceleration should be 
#used in a dimensionless form. Otherwise people might use units of m/s^2"""


def get_base_shear(G_e:float, S_p:float) -> float :
  """Returns the building base shear (Section 4.3.1, Static Method).

  Extended description of function.

  Args:
      G_e: seismic weight
      S_p: inelastic spectral acceleration in units of g

  Returns:
      building base shear (Section 4.3.1, Static Method)

  """
  return 0.85 * G_e * S_p

def get_vertical_force_distribution(
        base_shear, 
        storey_elevs, 
        storey_weights, 
        k=2,
    ):
    """"""
    # total_weight = sum(storey_weights)
    weight_height_power_k = [w * z**k for z, w in zip(storey_elevs, storey_weights)]
    sum_wz_k = sum(weight_height_power_k)
    return [base_shear * wz_k / sum_wz_k for wz_k in weight_height_power_k]


def storey_shear_moment(force_list, height_list, ascending = True):
    """Returns a tuple of lists containing storey shear and storey moment
    (from ground up, increasing z if 'ascending' is True)"""
    shear_list_rev = list(accumulate(force_list[::-1] if ascending else force_list))
    v_h_rev = [v * h for v, h in zip(shear_list_rev, height_list[::-1] if ascending else height_list)]
    mom_list_rev = list(accumulate([0]+v_h_rev[:-1]))
    return shear_list_rev[::-1] if ascending else shear_list_rev, mom_list_rev[::-1] if ascending else mom_list_rev


# def shear_coords(elevs, shears):
#     """Returns a list of (z,V) tuples for plotting storey shear diagrams
#     Inputs are lists of elevations and shear values
#     NB. Interleaving is required because of the steps in the shear diagram"""
#     l_1 = list(zip(elevs[1:],shears[1:]))
#     l_2 = list(zip(elevs,list(shears[1:]) + [0.0]))
#     out = []
#     while l_2:
#         out.append(l_2.pop()) if l_2 else None
#         out.append(l_1.pop()) if l_1 else None
    # return out[::-1]


#### Plotting Utilities #####


def tick_steps(a, scales = (1, 2, 4, 5, 10)):
    """"""
    d, m = (1, 0) if a==0 else divmod(log10(a),1) 
    for i in scales:
        log_i = log10(i)
        if m < log_i:
            break
    step = i * 10**(d-1)
    number = (int(a/step)+1)
    #k = tick * n
    return number, step


def x_lims(num_list):
    """Returns limits for plotting that match the ticks"""
    res = [tick_steps(a) for a in num_list]
    n_max = max(r for r,_ in res)
    return tuple([(0, n_max * step) for _, step in res])


def diagram_colours(colour_map:str, n_c:int, clip:float=0):
    """Colour and linetype settings for the plot_S_p plotting function.

    returns colour and linetype settings for the plot_S_p plotting function
    for discrete maps, it picks sequential colours. For continuous maps (like    
    viridis) it returns out evenly-spaced colours over the whole colour scale.

    Args:
            colour_map: the name (string) of a standard colour map (e.g. 'viridis' or 
            'tab10')
            n_c: numbers of colours to be returned
            clip: a way to cut off colours at the end of the scale which may be 
            too bright

    Returns:
            colour and linetype settings

    """
    cmap = plt.get_cmap(colour_map)    # viridis
    if colour_map.startswith('Paired'):
        return ('-', range(1,2*n_c,2)), ('--', range(0,2*n_c,2))
    elif colour_map.startswith('tab20b') or colour_map.startswith('tab20c'):
        return ('-', range(0,4*n_c,2)), ('--', range(2,4*n_c,2))
    elif colour_map.startswith('tab20'):
        return ('-', range(0,2*n_c,2)), ('--', range(1,2*n_c,2))
    elif cmap.N < 19:
        return ('-', range(n_c)), ('--', range(n_c))
    elif cmap.N < 25:
        return ('-', range(0,2*n_c,2)), ('--', range(1,2*n_c,2))
    else:
        x = max(0, clip)
        c_params = [0] if n_c < 2 else [float(i) / (n_c - 1 + x) for i in range(n_c)]
        return ('-', c_params), ('--', c_params)


def plot_S_p(soil_classes=['A', 'B', 'C', 'D'], q=1.5, gamma_ov=1.5, n_pt=10, 
                         a_gR=0.12, colour_map='tab10', T_vals=None, I_E=1.0, csvsave=None, **kwargs):
    """Returns the plot tuple (fig, axes) for the design spectrum

    This function simplifies the process of creating plots. 
    Note that kwargs are passed through to the plt.subplots function, which allows
    for more detailed control of the plots    
    
    Args:
        soil classes (list, tuple): a list of characters designating soil types (default is ['A', 'B', 'C', 'D'])
        q (float): force reduction factor (set to 1 to see elastic plots only) (default = 1.5)
        gamma_ov (float): material overstrength factor, section 6.2(3)a (default = 1.5)
        n_pt (integer): number of points per arc on the plot
        a_gR (float): reference ground acceleration (default is 0.12)
        colour_map: colour map used for curves (default is 'viridis') 
        T_vals (list, tuple): periods at which values should be plotted (optional)
        I_E (float): important class (I, II, III or IV) which determines 
                    the importance factor, 2.1(3)P & 2.1(4) (default = 1.0), 
                    which multiplies a_gR to return a_g
        csvsave (str): filename for csv output
    
    Returns:

    """
    # Set up colour maps
    cmap = plt.get_cmap(colour_map)    # viridis
    n_c = len(soil_classes)
    
    # The following parameters are for setting the chart colours
    (lt1, c_params1), (lt2, c_params2) = diagram_colours(colour_map, n_c)


    T_series = {sc: T_list_gen(sc, n_pt) for sc in soil_classes}
    # S_e is the elastic spectrum, S_p is the plastic spectrum
    S_e_series = {sc:[get_S_e(T, sc, a_gR=a_gR, importance_factor=I_E) 
        for T in T_series[sc]] for sc in soil_classes}
    S_p_series = {sc:[get_S_p(T, sc, a_gR=a_gR, q=q, gamma_ov=gamma_ov, importance_factor=I_E) 
        for T in T_series[sc]] for sc in soil_classes}
    

    # Write data to CSV files
    if csvsave is not None:
        with open(csvsave, 'w') as csvfile:
            csvwriter = writer(csvfile, lineterminator='\n')
            headers = sum([['T', f'Se({sc})'] for sc in soil_classes], []) + \
                                    sum([['T', f'Sd({sc}, q={q})'] for sc in soil_classes], [])
            csvwriter.writerow(headers)
            S_lists = sum([[T_series.get(sc), S_e_series.get(sc)] for sc in soil_classes], []) + \
                                    sum([[T_series.get(sc), S_p_series.get(sc)] for sc in soil_classes], [])
            [csvwriter.writerow(r) for r in zip(*S_lists)]


    # Generate Plots
    x_lim = kwargs.pop('x_lim', None)
    y_lim = kwargs.pop('y_lim', None)
    
    fig, axes = plt.subplots(**kwargs)
    if q > 1:
        [axes.plot(T_series[sc], S_e_series[sc], lt1, c=cmap(c)) for sc, c in zip(soil_classes, c_params1)]
        [axes.plot(T_series[sc], S_p_series[sc], lt2, c=cmap(c)) for sc, c in zip(soil_classes, c_params2)]
        axes.set_title('Horizontal Response Spectrum, 5% Damping' + ('' if I_E==1 else f' (I_E={I_E})'))
        if isinstance(T_vals, (list, tuple)):
            [axes.plot(T_vals, [get_S_e(T_val, sc, a_gR=a_gR, importance_factor=I_E) 
                    for T_val in T_vals], 'o', c=cmap(c)) 
                    for sc, c in zip(soil_classes, c_params1)]
            [axes.plot(T_vals, [get_S_p(T_val, sc, a_gR=a_gR, q=q, gamma_ov=gamma_ov, importance_factor=I_E) 
                    for T_val in T_vals], 'o', c=cmap(c)) 
                    for sc, c in zip(soil_classes, c_params2)]
    else:
        [axes.plot(T_series[sc], S_e_series[sc], lt1, c=cmap(c)) for sc, c in zip(soil_classes, c_params1)]
        if isinstance(T_vals, (list, tuple)):
                [axes.plot(T_vals, [get_S_e(T_val, sc, a_gR=a_gR) for T_val in T_vals], 'o', c=cmap(c)) 
                    for sc, c in zip(soil_classes, c_params1)]
        axes.set_title('Unscaled Elastic Horizontal Response Spectrum, 5% Damping')
    axes.legend(labels = soil_classes)
    #axes.title(f'')
    if x_lim is not None:
        axes.set_xlim(x_lim)
    if y_lim is not None:
        axes.set_ylim(y_lim)

    axes.grid(True)
    return fig, axes


def get_S_e_plot_data(soil_class):
    """"""
    T_series = T_list_gen(soil_class)
    S_series = [get_S_e(T, soil_class) for T in T_series]
    return T_series, S_series


def plot_VM(storey_elevs, storey_shears, storey_moments, figsize=(5,8), debug=False): # HVM_list):
    """"""
    #storey_elevs, storey_shears, storey_moments = zip(*HVM_list)
    M_units = 'kNm'
    if max(storey_moments) >= 1E9:
        M_units = 'GNm'
        storey_moments = [0.000001 * m for m in storey_moments]
    elif max(storey_moments) >= 1E6:
        M_units = 'MNm'
        storey_moments = [0.001 * m for m in storey_moments]
    
    max_shear = max(storey_shears)
    max_moment = max(storey_moments)
    
    shear_x_lims, moment_x_lims = x_lims((max_shear, max_moment))
    if debug:
        print(f'Plot Limits: Moment = {shear_x_lims}, Shear = {moment_x_lims}')

    elev_shear_pairs = shear_coords(tuple(storey_elevs), tuple(storey_shears))
    elevs, shears = zip(*elev_shear_pairs)
    if debug:
      print('  z  ,   V')
      [print(f'{z:5.4g}, {V:6,.5g}') for z, V in zip(storey_elevs, storey_shears)]
      print('  z  ,   V')
      [print(f'{z:5.4g}, {V:6,.5g}') for z, V in zip(elevs, shears)]
    
    fig, ax1 = plt.subplots(figsize=figsize)

    
    scolor = (0.5, 0.2, 0.1, 0.75)
    scolor1 = (0.5, 0.2, 0.1, 0.9)
    ax1.set_xlabel('Shear (kN)', color=scolor)
    ax1.set_ylabel('Elevation (m)')
    ax1.plot(shears, elevs, color=scolor1)
    ax1.set_ylim(0, None)
    ax1.set_xlim(*shear_x_lims)
    ax1.tick_params(axis='x', labelcolor=scolor)
    ax1.set_xticks(ax1.get_xticks())
    # label_format = '{:,.0f}'
    ax1.set_xticklabels(['{:,.0f}'.format(x) for x in ax1.get_xticks()])
    ax1.grid(True, which='both', axis='y')
    ax1.grid(True, which='major', axis='x', color=scolor)
    
    ax2 = ax1.twiny()    # instantiate a second axes that shares the same y-axis

    mcolor = (0.1, 0.2, 0.5, 0.75)
    mcolor1 = (0.1, 0.2, 0.5, 0.9)
    ax2.set_xlabel(f'Moment ({M_units})', color=mcolor)    # we already handled the y-label with ax1
    ax2.plot(storey_moments, storey_elevs, color=mcolor1)
    ax2.set_xlim(*moment_x_lims)
    ax2.tick_params(axis='x', labelcolor=mcolor)
    ax2.set_xticks(ax2.get_xticks())
    ax2.set_xticklabels(['{:,.0f}'.format(x) for x in ax2.get_xticks()])
    ax2.grid(True, which='major', axis='x', color=mcolor)

    fig.tight_layout()    # otherwise the right y-label is slightly clipped
    
##### Code-specific Print Utilities ########

def print_low_rise_data(storey_heights_m, storey_weights_kN, S_p_g, drift_ratio = 500):
    storey_masses = [w / 9.81 for w in storey_weights_kN]
    storey_elevs = list(accumulate(storey_heights_m))
    height = sum(storey_heights_m)
    weight = sum(storey_weights_kN)
    base_shear = get_base_shear(weight, S_p_g)
    T_code = get_T_code(height)
    
    print(f'\nBased on a drift ratio of H/{drift_ratio}, ')
    print(f'T_rayleigh would be {get_T_rayleigh(storey_weights_kN, delta_R=height/drift_ratio):4.2f} s\n')

    print(f'Height is {height:4.2f} m, weight is {weight:4,.0f} kN, T_code is ' + \
        f'{T_code:4.2f} s, \n\tS_p is {S_p_g:4.2f} g, base shear is {base_shear:4,.0f} kN')
    print(f'Base shear coefficient is {100 * base_shear / weight:4.2f}%')

    storey_forces = get_vertical_force_distribution(
        base_shear, storey_elevs, storey_weights_kN)
    storey_shears, storey_moments = storey_shear_moment(
        storey_forces, storey_heights_m, ascending = True)
    # base_moment = max(storey_moments)
    
    M_units = 'kNm'
    if max(storey_moments) >= 1E9:
        M_units = 'GNm'
        storey_moments = [0.000001 * m for m in storey_moments]
    elif max(storey_moments) >= 1E6:
        M_units = 'MNm'
        storey_moments = [0.001 * m for m in storey_moments]    
    
    print('Storey heights  (m): ', ' | '.join([f'{f:7,.2f}' for f in storey_heights_m]))
    print('Storey elevs    (m): ', ' | '.join([f'{f:7,.2f}' for f in storey_elevs]))
    print('Storey masses  (te): ', ' | '.join([f'{f:7,.1f}' for f in storey_masses]))
    print('Storey weights (kN): ', ' | '.join([f'{f:7,.0f}' for f in storey_weights_kN]))
    print('Storey forces  (kN): ', ' | '.join([f'{f:7,.0f}' for f in storey_forces]))
    print('Storey shears  (kN): ', ' | '.join([f'{f:7,.0f}' for f in storey_shears]))
    print(f'Storey mom\'ts ({M_units}): ', ' | '.join([f'{m:7,.0f}' for m in storey_moments]))
    print()
    print(f'Sum of storey forces is {sum(storey_forces):6,.0f} kN')
    print(f'Base Moment is {max(storey_moments):6,.0f} {M_units}')
    

def main():
    pass


if __name__ == "__main__":
    doctest.testmod()
    print('== Doctest Complete ==')
    main()

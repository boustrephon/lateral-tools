"""

"""

from math import pi, exp, log, log10, tan, sin, cos, atan2, isclose
from itertools import accumulate


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
  >>> isclose(exceedance2return(exceedance_probability=0.1, period=50), 474.56107905149526)
  True
  >>> isclose(exceedance2return(exceedance_probability=0.2, period=50), 224.07100588622754)
  True
  """
  return - period / log(1 - exceedance_probability)


def return2exceedance(return_period, period=50):
  """Returns the probability of a ground motion being exceeded over a period  
  when given a ground motion return period. For example, the probability of a  
  ground motion with return period of 475 years being exceeded over a time  
  period of 50 years is around 10%. 
  >>> isclose(return2exceedance(return_period=475, period=50), 0.09991237374774076)
  True
  >>> isclose(return2exceedance(return_period=225, period=50), 0.19926259708319194)
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


importance_class_dict = {
  'I': 'Civil protection - buildings and structures essential for civil protection.',
  'II':	'High consequences in case of failure - buildings and structures of high importance to the society.',
  'III':	'Ordinary - structures of ordinary importance and not classified as A, B and D.',
  'IV':	'Minor or temporary structures - structures of minor importance for public safety and temporary structures.'
}

importance_factor_dict = {
  'I': 1.7,
  'II': 1.4,
  'III': 1.0,
  'IV': None
  }


def get_importance_factor(importance_class, importance_factor_dict=importance_factor_dict):
  """Returns the Importance factor I_E
  Default values are provided, but can be replaced with NAD values
  by specifying the importance_factor_dict"""
  return importance_factor_dict.get(importance_class)


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

response_param_dict = {
    'A': {'S': 1.0, 'T_B': 0.12, 'T_C':0.32, 'T_D':2.5},
    'B': {'S': 1.6, 'T_B': 0.12, 'T_C':0.34, 'T_D':1.5},
    'C': {'S': 1.5, 'T_B': 0.22, 'T_C':0.50, 'T_D':1.5},
    'D': {'S': 1.5, 'T_B': 0.22, 'T_C':0.65, 'T_D':1.5}
}


gamma_ov_dict = {
    'Concrete':{'q':3, 'gamma_ov':1.5}, 
    'Composite':{'q':4, 'gamma_ov':1.5}, 
    'Steel_General':{'q':4, 'gamma_ov':1.25}, 
    'SCBF_DV':{'q':2, 'gamma_ov':1.25}, 
    'SCBF_K':{'q':1.5, 'gamma_ov':1.25},
    }


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
    params_dict = response_param_dict[soil_class]
    T_vals = [0] + [params_dict.get(p) for p in ['T_B', 'T_C', 'T_D']] + [10]
    lists = [[(T1 * float(n_pt - i) + T2 * float(i)) / n_pt for i in range(n_pt)] 
           for T1, T2 in zip(T_vals[:-1], T_vals[1:])]
    return sum(lists, []) + [10]

def get_a_g(a_gR, importance_factor):
    return a_gR * importance_factor


def get_S_e(
    T, 
    soil_class, 
    importance_factor=1.0, 
    S_T=1.0, 
    eta=1.0, 
    a_gR = 0.12, 
    response_param_dict=response_param_dict,
    ):
    """Returns Elastic Spectral Acceleration
    S_T - topographic amplification factor defined in clause 3.6
    eta - damping factor
    a_gR - reference ground spectral acceleration
    """
    params_dict = response_param_dict[soil_class]
    S, T_B, T_C, T_D = [params_dict.get(p) for p in ['S', 'T_B', 'T_C', 'T_D']]
    a_g = importance_factor * a_gR  # units are 'g'
    if T > 10:
        err_msg = f'Period T is too high (T = {T:.2f} > 10'
        raise ValueError(err_msg)
    elif T > T_D:
        return a_g * S * S_T * eta * 2.5 * T_C * T_D / T**2
    elif T > T_C:
        return a_g * S * S_T * eta * 2.5 * T_C / T
    elif T > T_B:
        return a_g * S * S_T * eta * 2.5
    else:
        return a_g * S * S_T * (1 + T / T_B * (eta * 2.5 - 1))


def get_S_p(
        T, 
        soil_class, 
        importance_factor=1.0, 
        S_T=1.0, 
        eta=1.0, 
        a_gR = 0.12, 
        q=1.5, 
        gamma_ov=1.5, 
        response_param_dict=response_param_dict,
    ):
    """"""
    T_C = response_param_dict[soil_class]['T_C']
    S_e = get_S_e(
        T, 
        soil_class, 
        importance_factor=importance_factor, 
        S_T=S_T, 
        eta=eta, 
        a_gR = a_gR, 
        response_param_dict=response_param_dict,
      )
    S_p = S_e / (q if T > T_C else (gamma_ov + (q - gamma_ov) * T/T_C))
    return S_p


"""I think it would be helpful to make clear that ground acceleration should be 
used in a dimensionless form. Otherwise people might use units of m/s^2"""


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


def storey_shear_moment(force_list, height_list, rev = True):
    """Returns a tuple of lists containing storey shear and storey moment
    (from ground up, increasing z if 'rev' is True)"""
    shear_list_rev = list(accumulate(force_list[::-1] if rev else force_list))
    v_h_rev = [v * h for v, h in zip(shear_list_rev, height_list[::-1] if rev else height_list)]
    mom_list_rev = list(accumulate([0]+v_h_rev[:-1]))
    return shear_list_rev[::-1] if rev else shear_list_rev, mom_list_rev[::-1] if rev else mom_list_rev


def shear_coords(elevs, shears):
    """Returns a list of (z,V) tuples for plotting storey shear diagrams
    Inputs are lists of elevations and shear values
    NB. Interleaving is required because of the steps in the shear diagram"""
    l_1 = list(zip(elevs[1:],shears[1:]))
    l_2 = list(zip(elevs,list(shears[1:]) + [0.0]))
    out = []
    while l_2:
        out.append(l_2.pop()) if l_2 else None
        out.append(l_1.pop()) if l_1 else None
    return out[::-1]


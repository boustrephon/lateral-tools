"""A library with general seismic utilities
"""

# from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
from math import exp, log
from itertools import accumulate


##### Basic functions #####

# Functions:
# printt, interp, list_interp, 


def printt(table, h1 = 20, h2 = 12, h3 = 4):
    """
    Printing tables where first column is text and the remainder are floats
    h1 is the width of the text column, 
    h2 is the width of all other columns
    h3 is the number of decimal places
    """
    for line in table:
        print(f'{line[0]:{h1}s}', end = ' ')
        [print(f'{item:{h2-h3}.{h3}f}', end = ' ') for item in line[1:]]
    return None


def interp(x, x1, x2, y1, y2):
    """Returns an interpolated y-value for a line 
    between two points (x1, y1) and (x2, y2)
    for a given x-value - i.e. linear interpolation 
    along a line.
    """
    if y1 is None or y2 is None:
        return None
    elif x1 == x2:
        return 0.5 * (y1 + y2)
    else:
        return y1 + (x - x1) / (x2 - x1) * (y2 - y1)


def list_interp(val, list_1, list_2):
    """Returns interpolated values from list_2 based 
    on values related to list_1, returning the bookends 
    if the provided values are outside the range of 
    list_1. If any of the values are None, it will 
    return None.
    
    val is the lookup value
    list_1 is the lookup list
    list_2 is the result list with values corresponding 
        to those on list_1
    
    >>> list_interp(val, list_1, list_2)
    """
    i_list = [i for i, n in enumerate(list_1) if n == val]
    if i_list:  
        # if lookup value matches a value in list_1
        return list_2[i_list[0]]
    elif val <= list_1[0]:
        # if lookup value is lower than all in list_1
        return list_2[0]
    elif val >= list_1[-1]:
        # if lookup value is higher than all in list_1
        return list_2[-1]
    else:
        # carry out interpolation
        index_list = [i for i, (x1, x2) 
                      in enumerate(zip(list_1[:-1], list_1[1:])) 
                      if val >= x1 and val <= x2]
        if len(index_list) == 1:
            j = index_list[0]
            vals = (list_1[j], list_1[j+1], 
                    list_2[j], list_2[j+1])
            if any([v is None for v in vals]):
                return None
            else:
                return interp(val, *vals)
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
  """"Returns a dictionary containing S_s & S_1 values 
  from a reference 475-year return PGA.
  
  Reference: 15WCEE, 'Deriving S_s and S_1 parameters 
    from PGA maps'
  """
  def S_s_calc(PGA): return 0.3386 * PGA + 2.1696
  def S_1_calc(PGA): return 0.5776 * PGA + 0.5967
  return {'S_s':S_s_calc(a_g_475), 'S_1':S_1_calc(a_g_475)}


def S_calc(S_s_475, S_1_475, T):
  """"Returns a dictionary containing S_s, S_1, 
  k_s & k_1 values for a given return period based 
  on a reference PGA and reference return period.
  
  Reference: 15WCEE, 'Deriving S_s and S_1 parameters 
    from PGA maps'
  """
  def S_func(S_475, k, T): 
      return S_475 * (T / 475)**k
  
  S_s_475_S_keys = (0.20, 0.55, 0.90, 1.20, 1.50)
  S_s_475_k_vals = (0.46, 0.42, 0.39, 0.35, 0.30)
  S_1_475_S_keys = (0.01, 0.02, 0.05, 0.10, 0.20)
  S_1_475_k_vals = (0.40, 0.40, 0.40, 0.40, 0.40)
  
  k_s = list_interp(S_s_475, S_s_475_S_keys, S_s_475_k_vals)
  S_s = S_func(S_s_475, k_s, T)
  
  k_1 = list_interp(S_1_475, S_1_475_S_keys, S_1_475_k_vals)
  S_1 = S_func(S_1_475, k_1, T)
  
  return {'S_s':S_s, 'S_1':S_1, 'k_s':k_s, 'k_1':k_1}



#### Storey-based functions #####

# Functions:
# storey_force_coeffs, storey_force_distribution, storey_forces, storey_shear_moment, shear_coords
# storey_force_coeffs, storey_force_distribution, storey_shear_moment, shear_coords

def storey_force_coeffs(elev_list, weight_list, k = 1.0):
    """
        Inputs:
            elev_list - list of storey elevations (m)
          weight_list - list of storey weights (kN)
    """
    divisor = sum(w * h**k for w, h in zip(weight_list, elev_list))
    return [w * h**k / divisor for w, h in zip(weight_list, elev_list)]


def storey_force_distribution(V, elev_list, weight_list, k = 1.0):
    """
        Inputs:
            V  - base shear in kN
            elev_list - list of storey elevations (m)
          weight_list - list of storey weights (kN)
    """
    divisor = sum(w * h**k for w, h in zip(weight_list, elev_list))
    return [w * h**k * V / divisor for w, h in zip(weight_list, elev_list)]


def storey_forces(force_list, rev = True):
    """Returns a list of accumulated storey forces - i.e. storey shears
    NB lists have to be reversed before applying the accumulate function"""
    shear_list_rev = accumulate(force_list[::-1] if rev else force_list)
    return list(shear_list_rev)[::-1] if rev else list(shear_list_rev)


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
    l_2 = list(zip(elevs,shears[1:] + [0.0]))
    out = []
    while l_2:
        out.append(l_2.pop()) if l_2 else None
        out.append(l_1.pop()) if l_1 else None
    return out[::-1]


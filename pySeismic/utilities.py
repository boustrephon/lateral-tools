"""A library with general seismic utilities
"""

# from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
from math import exp, log
from itertools import accumulate
from typing import Any

##### Basic functions #####

# Functions:
# printt, interp, list_interp, 

def as_float(x):
    try:
        return float(x)        
    except:
        return x


units_SI_scaling_dict = {
    'mi':1604.344, 'yd':0.01, 'ft':0.3048, 'in':0.0254, 
    'km':1000.0, 'm':1.0, 'cm':0.01, 'mm':0.001,
    'GN':1.0E9, 'MN':1.0E6, 'kN':1000.0, 'N':1.0,
    'GNm':1.0E9, 'MNm':1.0E6, 'kNm':1000.0, 'Nm':1.0, 'Nmm':0.001,
    'kN/m':1000.0, 'N/m':1.0, 'N/mm':1000.0,
    'kN/m2':1000.0, 'N/m2':1.0, 'N/mm2':1.0E6,
    'kN/m^2':1000.0, 'N/m^2':1.0, 'N/mm^2':1.0E6,
    'ksf': 4, 'psf'
    'GPa':1.0E9, 'MPa':1.0E6, 'kPa':1000.0, 'Pa':1.0,
    'ft/s':0.3048, 'in/s':0.0254, 
    'km/s':1000.0, 'm/s':1.0, 'cm/s':0.01, 'mm/s':0.001,
    'ft/s2':0.3048, 'in/s2':0.0254, 'ft/s^2':0.3048, 'in/s^2':0.0254, 
    'km/s2':1000.0, 'm/s2':1.0, 'cm/s2':0.01, 'mm/s2':0.001,
    'km/s^2':1000.0, 'm/s^2':1.0, 'cm/s^2':0.01, 'mm/s^2':0.001,
    'g':9.80665
}


def scale_to_SI_units(units: str, default:Any=None) -> float:
    """Returns the factor by which a value with a certain unit
    should be scaled to get to SI units
    
    Args:
        unit: name of scale factor. Note that this 
        will have problems with 
        default: default result if lookup fails

    Returns: 
        scale: scaling factor   
    """
    scale = units_SI_scaling_dict.get(units, None)
    if scale is None:
        casefold_dict = {k.casefold():v 
                            for k, v in units_SI_scaling_dict.items()}
        scale = casefold_dict.get(units.casefold(), default)
    return scale


def printt(table, h1 = 20, h2 = 12):
    """
    Printing tables where first column is text and the remainder are floats
    h1 is the 
    """
    for line in table:
        print('{:20s}'.format(line[0]), end = ' ')
        [print('{:8.4f}'.format(item), end = ' ') for item in line[1:]]
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


def list_interp(val, list_1, list_2, extrapolate=False):
    """Returns interpolated values from list_2 based 
    on values related to list_1, returning the bookends 
    if the provided values are outside the range of 
    list_1 and if `extrapolate` is False. If any of 
    the values are None, it will return `None`.
    
    val is the lookup value
    list_1 is the lookup list
    list_2 is the result list with values corresponding 
        to those on list_1
    extrapolate is the option to extrapolate linearly
    
    >>> list_interp(val, list_1, list_2)
    """
    i_list = [i for i, n in enumerate(list_1) if n == val]
    if i_list:  
        # if lookup value matches a value in list_1
        return list_2[i_list[0]]
    elif val <= list_1[0]:
        # if lookup value is lower than all in list_1
        if extrapolate:
            vals = list_1[:2] + list_2[:2]
            return interp(val, *vals)
        else:
            return list_2[0]
    elif val >= list_1[-1]:
        # if lookup value is higher than all in list_1
        if extrapolate:
            vals = list_1[-2:] + list_2[-2:]
            return interp(val, *vals)
        else:
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


def storey_forces(force_list, ascending = True):
    """Returns a list of accumulated storey forces - i.e. storey shears
    NB lists have to be reversed before applying the accumulate function"""
    shear_list_desc = accumulate(force_list[::-1] if ascending else force_list)
    return list(shear_list_desc)[::-1] if ascending else list(shear_list_desc)


def storey_shear_moment(force_list, height_list, ascending = True, debug=False):
    """
    Returns a tuple of lists containing storey shear and storey moment
    (from ground up, increasing z if 'ascending' is True)

    Note that the shear and moment should be considered to be at
    the bottom of each storey
    """
    shear_list_desc = list(accumulate(force_list[::-1] if ascending else force_list))
    v_h_desc = [v * h for v, h in zip(shear_list_desc, height_list[::-1] if ascending else height_list)]
    mom_list_desc = list(accumulate(v_h_desc))
    if debug: print(mom_list_desc[-1], '=?', sum(f * z for f, z in zip(
            accumulate(height_list) if ascending else accumulate(height_list[::-1]), 
            force_list if ascending else force_list[::-1])))
    return shear_list_desc[::-1] if ascending else shear_list_desc, mom_list_desc[::-1] if ascending else mom_list_desc


def shear_coords(elevs, shears, base_level = 0):
    """Returns a list of (z,V) tuples for plotting storey shear diagrams
    Inputs are lists of elevations and shear values.
    The list of elevations should not include the base level
    NB. Interleaving is required because of the steps in the shear diagram"""
    # l_1 = list(zip(elevs[1:],shears[1:]))
    # l_2 = list(zip(elevs,shears[1:] + [0.0]))
    l_1 = list(zip(elevs,shears))
    l_2 = list(zip(list(elevs)[1:] + [base_level], shears))
    out = []
    while l_2:
        out.append(l_2.pop()) if l_2 else None
        out.append(l_1.pop()) if l_1 else None
    return out[::-1]


# =======================================
# ======= Printing Functions ============
# =======================================

from os.path import exists, abspath
from itertools import accumulate
from csv import reader, excel_tab
from tkinter import filedialog
from tkinter import *


def import_storey_data(file_path, delimiter=','):
    """
    Importing 4-column data from CSV file
    First column is assumed to be keys (storey names) and
    other three are required to be numeric. Any row containing
    non-numeric data in the last three columns will be ignored
    """
    if file_path is None:
        root = Tk()
        root.filename =  filedialog.askopenfilename(
            initialdir = "/", title = "Select file", 
            filetypes = (("csv files","*.csv"),("all files","*.*"))
            )
        file_path = root.filename
        root.destroy()

    print(f'Building data file path is "{file_path}"')
    print(f'File EXISTS' if exists(file_path) else f'File NOT found')

    with open(file_path, newline='') as csvfile:
        csv_reader = reader(csvfile, dialect = 'excel', delimiter=delimiter)
        data = [r for r in csv_reader]
    
    # _, b, c, d = data[0]
    # storey_headers = data[0] if not all([isinstance(as_float(x), float) 
    #                                         for x in (b, c, d)]) else []
    
    storey_data = [[a, float(b), float(c), float(d)] 
                for a,b,c,d in data if all([isinstance(as_float(x), float) 
                                            for x in (b, c, d)])]
    
    return storey_data


def print_table(table, h1 = 20, h2 = 12, h3 = 4):
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


def print_storey_data(
        storey_elevs_m, 
        storey_weights_kN, 
        storey_names = [], 
        base_level_m = 0,
        ascending = True):
    """Prints out storey data
    All data should be for ascending height. The base level should 
    not be included in the data, except as the `base_level_m` value.
    
    """
    if ascending:
        storey_elevs_m.reverse()
        storey_weights_kN.reverse()
        storey_names.reverse()
    
    storey_heights_m = [h2 - h1 for h2, h1 in zip(storey_elevs_m, list(storey_elevs_m[1:]) + [base_level_m])]
    storey_masses_kg = [1000 * w / 9.80665 for w in storey_weights_kN] # convert kN to kg

    if len(storey_weights_kN) == len(storey_elevs_m) == len(storey_names):
        pass
    else:
        storey_names = [f'L{i:04d}' for i, _ in enumerate(storey_elevs_m)]
        storey_names.reverse()
        
    basic_f = 'Height is {:.2f}m, total mass is {:8,.0f} tonne, total weight is {:8,.0f} kN\n'
    print(basic_f.format(max(storey_elevs_m), 0.001 * sum(storey_masses_kg), sum(storey_weights_kN)))

    tit_f1 = '{:20s} {:9s} {:10s} {:10s} {:10s}'
    tit_f2 = '{:20s} {:9s} {:10s} {:>10,.0f} {:8,.0f}'
    print(tit_f1.format('Storey', 'Ht. (m)', 'Elev (m)', 'Mass (kg)', 'Weight (kN)'))
    num_f = '{:20s} {:6.3f}   {:8.3f}   {:>10,.0f}   {:8,.0f}'
    [print(num_f.format(n, h, e, m, w)) for n, h, e, m, w in zip(storey_names, storey_heights_m, 
                                                    storey_elevs_m, storey_masses_kg, storey_weights_kN)]
    # base_level = f'{base_level_m:8.3f}'
    # print(base_level)
    print(tit_f2.format('BASE', '------', f'{base_level_m:7.3f}  ', sum(storey_masses_kg), sum(storey_weights_kN)))
    print()

    return None


def print_storey_force_data(
        V_kN, 
        storey_elevs_m, 
        storey_weights_kN, 
        storey_names = [], 
        k = 1,
        base_level_m = 0,
        ascending = True
        ):
    """Calculates and prints out storey force data
    All data should be for ascending height. The base level should 
    not be included in the data, except as the `base_level_m` value.
    
    """
    if ascending:
        list(storey_elevs_m).reverse()
        list(storey_weights_kN).reverse()
        list(storey_names).reverse()
    
    storey_heights_m = [h2 - h1 for h2, h1 in zip(storey_elevs_m, list(storey_elevs_m[1:]) + [base_level_m])]
    # storey_masses_kg = [1000 * w / 9.80665 for w in storey_weights_kN] # convert kN to kg

    if len(storey_weights_kN) == len(storey_elevs_m) == len(storey_names):
        pass
    else:
        storey_names = [f'L{i:04d}' for i, _ in enumerate(storey_elevs_m)]
        storey_names.reverse()

    k = k
    storey_force_coeff_list = storey_force_coeffs(storey_elevs_m, storey_weights_kN, k)
    storey_loads_kN = storey_force_distribution(V_kN, storey_elevs_m, storey_weights_kN, k)

    storey_shears_kN, storey_moments_kNm = storey_shear_moment(storey_loads_kN, storey_heights_m, False)

    tit_f1 = '{:20s} {:9s} {:10s} {:10s} {:10s} {:10s} {:10s}'
    print(tit_f1.format('Storey', 'Elevs (m)', 'Weights (kN)', ' Coeffs ', 'Forces (kN)', 'Shears (kN)', 'Moments (kNm)'))
    num_f = '{:20s} {:8.2f}   {:8,.0f}   {:8.4f}   {:8,.0f}   {:>10,.0f}   {:>10,.0f}'
    [print(num_f.format(n, h, wt, c, f, s, m)) for n, h, wt, c, f, s, m in zip(storey_names, storey_elevs_m, storey_weights_kN, 
                            storey_force_coeff_list, storey_loads_kN, storey_shears_kN, [0] + storey_moments_kNm[:-1])]
    tit_f2 = '{:20s} {:8.2f} {:10s} {:12s} {:10s} {:12s} {:>10,.0f}'
    print(tit_f2.format('BASE', base_level_m, '    ------', '    ------', ' ------   ', '  ------  ', storey_moments_kNm[-1]))
    tit_f3 = '{:20s} {:10s} {:8,.0f} {:10.4f} {:10,.0f} {:14s} {:12s}'
    print(tit_f3.format('SUM ', '  ------', sum(storey_weights_kN), sum(storey_force_coeff_list), sum(storey_loads_kN), '     ------', '     ------'))
    
    return None


def print_bar_areas_given_spacing(bar_spacing_mm, 
                       bar_list = [16, 20, 25, 32, 36, 40]):
    
    [print(f'T{D}@200 = {1000 * D**2 / 4 / bar_spacing_mm}mm2/m') 
        for D in bar_list]
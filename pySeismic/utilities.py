"""A library with general seismic utilities
"""

# from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
from itertools import accumulate


##### Basic functions #####

# Functions:
# printt, interp, list_interp, 


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
    """
    Linear interpolation based on a line between 
    two points (x1, y1) & (x2, y2)
    """
    if y1 is None or y2 is None:
        return None
    elif x1 == x2:
        return 0.5 * (y1 + y2)
    else:
        return y1 + (x - x1) / (x2 - x1) * (y2 - y1)


def list_interp(val, list_1, list_2):
    """Returns interpolated values from list_2 based on values related to list_1, returning the bookends if the 
    provided values are outside the range of list_1
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


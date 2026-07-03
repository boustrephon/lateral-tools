"""
analysis.py — Storey and member structural analysis functions.

Part of the lateral-tools package for seismic and lateral load engineering.

This module provides:
- Storey analysis: storey_elevations, storey_force_coeffs, storey_force_distribution,
  storey_forces, storey_shear_moment, shear_coords
- Member shear analysis: split_shear_list, shear_interp

These functions compute lateral force distribution and internal forces
for multi-storey buildings under seismic or wind loading.
"""
from __future__ import annotations
from itertools import accumulate
from typing import Any, Sequence
def storey_elevations(storey_height_list: list[float], 
                      ascending:bool=True, 
                      rounded = False,
                      rounding_value = 10) -> list[float]:
    """list of storey elevations from list of storey heights

    Args:
        storey_height_list (list[float]): ordered list of storey heights
        ascending (bool, optional): True if list is ascending 
        rounded (bool, optional): True if values should be rounded to clean numbers
        rounding_value (int, optional): number of decimal points to round to. Defaults to 10.
    Returns:
        list[float]: list of storey elevations (m)
    """
    if ascending is False:
        storey_height_list = storey_height_list[::-1]
    storey_height_list = storey_height_list if storey_height_list[0] == 0 else [0] + storey_height_list
    elev_list = list(accumulate(storey_height_list))
    rounding_value = 10
    if rounded:
        return [round(elev, rounding_value) for elev in elev_list] \
            if ascending else [round(elev, rounding_value) for elev in elev_list[::-1]]
    else:
        return elev_list if ascending else elev_list[::-1]




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




def storey_shear_moment(force_list, height_list, debug=False):
    """
    Returns a tuple of lists containing storey shear and storey moment
    (from ground up, increasing z if 'ascending' is True)

    Note that the shear and moment should be considered to be at
    the bottom of each storey

    Args:
        force_list (list[float]): list of storey forces (kN)
        height_list (list[float]): list of storey heights (m)
    
    Returns:
        tuple[list[float], list[float]]: shear_list, moment_list

    Example:
        >>> shear_list, moment_list = storey_shear_moment([1000, 1500, 2000], [10, 20, 30])

    """
    if ascending := all([h2 > h1 for h1, h2 in zip(height_list[:-1], height_list[1:])]):
        force_list = force_list[::-1]
        height_list = height_list[::-1]
    elif ascending := not all([h2 < h1 for h1, h2 in zip(height_list[:-1], height_list[1:])]):
        pass
    else:
        raise ValueError("height_list is not monotonic")
    
    shear_list_desc = list(accumulate(force_list))
    heights = [h2 - h1 for h1, h2 in zip(height_list[1:], height_list[:-1])]
    v_h_desc = [v * h for v, h in zip(shear_list_desc, heights)]
    mom_list_desc = list(accumulate(v_h_desc))
    if debug:
        print(mom_list_desc[-1], '=?', sum(f * z for f, z in zip(
            accumulate(height_list) if ascending else accumulate(height_list[::-1]), 
            force_list if ascending else force_list[::-1])))
    return shear_list_desc[::-1] if ascending else shear_list_desc, mom_list_desc[::-1] if ascending else mom_list_desc




def shear_coords(elevs: Sequence[float], shears:Sequence[float], 
                 base_level:float = 0, tolerance:float = 1e-6) -> Sequence[Sequence[float]]:
    """Returns a list of (z,V) tuples for plotting storey shear diagrams
    Inputs are lists of ascending elevations and matching shear values.
    If the list of elevations includes the base level,
    the first item of each list will be removed
    NB. Interleaving is required because of the steps in the shear diagram

    Args:
        elevs (list[float]): storey elevations (ascending)
        shears (list[float]): shear values at those elevations
        base_level (float, optional): Elevation at base. Defaults to 0.

    Raises:
        ValueError: Length mismatch between elevs and shears

    Returns:
        list[tuple[float, float]]: List of elevations and shears
    """
    # If the base level is included in the list of elevations
    if abs(elevs[0] - base_level) < tolerance:
        elevs = elevs[1:]
        # If there is then one more storey shear than elevation:
        if len(shears) == len(elevs):
            pass
        # elif len(shears) - len(elevs) == 1:
        #     err_msg = f'Length mismatch between elevs {len(elevs)} and shears {len(shears)}'
        #     # shears = shears[1:]
        #     raise ValueError(err_msg)
        else:
            err_msg = f'Length mismatch between elevs {len(elevs)} and shears {len(shears)}'
            err_msg += f'elevs = {elevs}\nshears = {shears}'
            raise ValueError(err_msg)
    out = [(base_level, 0.0)]
    z_old = base_level
    for z, V in zip(elevs, shears):
        out.append((z_old, V))
        out.append((z, V))
        z_old = z
    out.append((elevs[-1], 0.0))
    out.append((base_level, 0.0))
    return out


# =======================================
# ===== Member-based Functions ==========
# =======================================




def split_shear_list(xV_pairs: list[tuple[float, float]], tol:float=1E-8) -> list[list[tuple[float, float]]]:
    """Splits a list of x-V pairs wherever the x-value is repeated, so that the
    list becomes a list of lists.

    Args:
        xV_pairs (list[tuple[float, float]]): pairs of normalised position and shear 
            along a beam for ascending values of x between 0.0 and 1.0.
        tol (float, optional): tolerances. Defaults to 1E-8.

    Returns:
        list[list[tuple[float, float]]]: A list of lists of x-V pairs where no x-values are repeated.
    """
    xV_pairs_split_list = []
    xV_pairs_sublist = []
    x2 = None
    y2 = None
    for (x1, y1), (x2,y2) in zip(xV_pairs[:-1], xV_pairs[1:]):
        xV_pairs_sublist.append((x1, y1))
        if abs(x1 - x2) < tol:
            # add current sublist to output list
            xV_pairs_split_list.append(xV_pairs_sublist)
            # reset sublist to empty
            xV_pairs_sublist = []
    if x2 is not None and y2 is not None:
        xV_pairs_sublist.append((x2, y2))
        xV_pairs_split_list.append(xV_pairs_sublist) 
    return xV_pairs_split_list




def shear_interp(x_vals: list[float], xV_pairs: list[tuple[float, float]], tol: float = 1E-6) -> list[tuple[float, float | None]]:    
    """Interpolated values position and shear generated for specified normalised position values
    when provided with position and shear values.

    No checks are carried out regarding validity of data. Possible checks include:
    - Check that shear discontinuities occur at points in the list of x_vals
    - Check that x-values are ascending
    - Check that x-values are in the range of 0.0 to 1.0 inclusive
    
    Args:
        x_vals (list[float]): list of normalised position parameters, x (from 0.0 to 1.0)
        xV_pairs (list[tuple[float, float]]): list of (x, V) tuples where x ranges from 0.0 to 1.0 (ascending)
        tol (float, optional): how close . Defaults to 1E-6.

    Returns:
        list[tuple[float, float]]: New list of xV pairs based on x_vals
    
    Examples:
        >>> shear_interp([0.0, 0.25, 0.5, 0.75, 1.0], [(0.00, -5000.0), (0.1667, -4500.0), \
            (0.3333, -4000.0), (0.5000, -3500.0), (0.5000, 3000.0), (0.6667, 3500.0), (0.8333, 4000.0), (1.000, 4500.0)])
        [(0.0, -5000.0), (0.25, -4250.0), (0.5, -3500.0), (0.5, 3000.0), (0.75, 3750.0), (1.0, 4500.0)]
    """
    split_shear_list_pairs = split_shear_list(xV_pairs)
    x_and_y_lists = [[thing for thing in list(zip(*shear_list))] for shear_list in split_shear_list_pairs]
    interpolated_split_xV_pairs = [[(x, list_interp(x, x_list, y_list)) for x in x_vals if (x >= min(x_list) - tol and x <= max(x_list) + tol)] for x_list, y_list in x_and_y_lists]
    interpolated_xV_pairs = list(chain(*interpolated_split_xV_pairs))
    # This has problems if list_interp returns None
    return interpolated_xV_pairs


# =======================================
# ======= Printing Functions ============
# =======================================

# from os.path import exists, abspath
# from itertools import accumulate
# from csv import reader, excel_tab
# from tkinter import filedialog
# from tkinter import *





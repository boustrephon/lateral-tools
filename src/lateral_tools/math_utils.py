"""
math_utils.py — General mathematical utility functions for lateral-tools.

Part of the lateral-tools package for seismic and lateral load engineering.

This module provides:
- Type conversion helpers: as_float, as_numeric, scale_to_SI_units
- Interpolation functions: interp, list_interp

These are low-level utilities used across the package.
"""
from __future__ import annotations
from math import isclose
from typing import Any, Union, Sequence
def as_float(txt:str, default:Any=None) -> Union[float, str]:
    """Converts string to number (float), if possible
    
    Args:
        txt: text for conversion to float
        default: default value if conversion fails
    
    Returns:
        text as float or default value

    >>> as_float('6')
    6.0
    >>> as_float('3.14')
    3.14
    >>> as_float('3.14mm')
    '3.14mm'
    >>> as_float('3.14mm', 0)
    0
    """
    scale=1
    if txt.rstrip().endswith('%'):
        txt = txt.rstrip().replace('%', '')
        scale=0.01
    try:
        return scale*float(txt)       
    except ValueError:
        if default is None:
            if scale == 1:
                return txt
            else:
                return txt + '%'
        else:
            return default




def as_numeric(txt:str, default:Any=None) -> Union[int, float, str, None]:
    """Converts string to number (int or float), if possible
        
    Args:
        txt: text for conversion to integer or float
        default: default value if conversion fails
    
    Returns:
        text as int, float or default value

    >>> as_numeric('6')
    6
    >>> as_numeric('3.14')
    3.14
    >>> as_numeric('23.3%')
    0.233
    >>> as_numeric('3.14mm')
    '3.14mm'
    >>> as_numeric('3.14mm', 0)
    0
    """
    scale=1
    if txt.rstrip().endswith('%'):
        txt = txt.rstrip().replace('%', '')
        scale=0.01
    try:
        return scale*int(txt)
    except ValueError:
        try:
            return scale*float(txt)
        except ValueError:
            if default is None:
                if scale == 1:
                    return txt
                else:
                    return txt + '%'
            else:
                return default


units_SI_scaling_dict = {
    'mi':1609.344, 'yd':0.9144, 'ft':0.3048, 'in':0.0254, 
    'km':1000.0, 'm':1.0, 'cm':0.01, 'mm':0.001,
    'GN':1.0E9, 'MN':1.0E6, 'kN':1000.0, 'N':1.0,
    'GNm':1.0E9, 'MNm':1.0E6, 'kNm':1000.0, 'Nm':1.0, 'Nmm':0.001,
    'kN/m':1000.0, 'N/m':1.0, 'N/mm':1000.0,
    'kN/m2':1000.0, 'N/m2':1.0, 'N/mm2':1.0E6,
    'kN/m^2':1000.0, 'N/m^2':1.0, 'N/mm^2':1.0E6,
    'ksf': 47880.26, 'psf': 47.88026,
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




def interp(x:float, x1: float, x2: float, y1: float, y2: float) -> float:
    """Returns an interpolated y-value for a line 
    between two points (x1, y1) and (x2, y2)
    for a given x-value - i.e. linear interpolation 
    along a line.

    Args:
        x (float): point at which interpolation is to occur
        x1 (float): start value of x
        x2 (float): end value of x
        y1 (float): start value of y
        y2 (float): end value of y

    Returns:
        float: interpolated y-value for a line at point x
    
    Example:
        >>> interp(1.5, 0.7, 1.9, 1.4, 2.3)
        2.0
        >>> round(interp(2.3, 0.7, 1.9, 1.4, 2.3), 6)
        2.6
        >>> interp(-0.1, 0.7, 1.9, 1.4, 2.3)
        0.8
    """
    if y1 is None or y2 is None:
        return None
    elif x1 == x2:
        return 0.5 * (y1 + y2)
    else:
        return y1 + (x - x1) / (x2 - x1) * (y2 - y1)




def list_interp(val: float, 
                list_1: list[float]|tuple[float, ...], 
                list_2: list[float]|tuple[float, ...], 
                extend:bool = False, extrapolate:bool=False)-> float | None:
    """Returns interpolated values from list_2 based 
    on values related to list_1, returning the bookends 
    if the provided values are outside the range of 
    list_1 and if `extrapolate` is False. If any of 
    the values are None, it will return `None`.
    
    Args:
        val (float): the lookup value
        list_1 (list[float]): the lookup list
        list_2 (list[float]): is the result list with values corresponding 
                to those on list_1
        extend (Bool, optional): the option to extend the values 
            outside the limits of list_1 with the bookend values. Defaults to False.
        extrapolate (Bool, optional): the option to extrapolate linearly 
            outside the limits of list_1. Defaults to False.

    Returns:
        float | None: Linearly interpolated function
        
    Examples:
        >>> list_interp(0.5, [0.2, 0.8], [1.1, 1.35])
        1.225
        >>> list_interp(0.08, [0.2, 0.8], [1.1, 1.35], True, True)
        1.05
        >>> list_interp(0.08, [0.2, 0.8], [1.1, 1.35], True, False)
        1.1
        >>> list_interp(0.08, [0.2, 0.8], [1.1, 1.35], False)
        0
    """
    i_list = [i for i, n in enumerate(list_1) if n == val]
    if i_list:  
        # if lookup value matches a value in list_1
        return list_2[i_list[0]]
    elif val <= list_1[0]:
        # if lookup value is lower than all in list_1
        if not extend:
            return 0
        elif extrapolate:
            vals = tuple(list_1[:2]) + tuple(list_2[:2])
            return interp(val, *vals)
        else:
            return list_2[0]
    elif val >= list_1[-1]:
        # if lookup value is higher than all in list_1
        if not extend:
            return 0
        elif extrapolate:
            vals = tuple(list_1[-2:]) + tuple(list_2[-2:])
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





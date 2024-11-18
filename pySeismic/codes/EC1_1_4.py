"""
Calculations for Eurocode 1991-1-4 - Wind Loads
"""

import sys
from os.path import join as path_join, dirname, abspath
import doctest

# Add the parent directory of 'pySeismic' to sys.path
sys.path.append(abspath(path_join(dirname(__file__), '..', '..')))

# print('sys_path:', sys.path)

from math import log10
from pySeismic.utilities import list_interp  # noqa: E402
# from ..utilities import list_interp


def get_c_f0(d, b):
    """_summary_

    Args:
        d (float): width
        b (float): vertical dimension

    Returns:
        float: _description_
    
    Example:
        >>> get_c_f0(9.2, 10.0)
        2.18
        >>> get_c_f0(9.0, 1.0)
        0.92
        >>> [get_c_f0(d, 1) for d in (0.05, 0.1, 0.2, 0.7, 1.0, 2.0, 5.0, 10.0, 50.0, 60.0)]
        [2.0, 2.0, 2.0, 2.4, 2.1, 1.65, 1.0, 0.9, 0.9, 0.9]
    """
    d_b = (0.1, 0.2, 0.7, 1.0, 2.0, 5.0, 10.0, 50.0)
    c_f0 = (2.0, 2.0, 2.4, 2.1, 1.65, 1.0, 0.9, 0.9)
    return list_interp(d/b, d_b, c_f0, extend=True)


FLAT_ROOF_COEFF_DICT = {
    'sharp eaves': {'hp_h': 0.0, 'r_h': 0.0, 'alpha': 90, 'F_10': -1.8, 'G_10': -1.2, 'H_10': -0.7, 'Ip_10': +0.2, 'In_10': -0.2,
                                                               'F_1': -2.5, 'G_1': -2.0, 'H_1': -1.2, 'Ip_1': +0.2, 'In_1': -0.2,},
    'parapets 25': {'hp_h': 0.025, 'r_h': 0.0, 'alpha': 90, 'F_10': -1.6, 'G_10': -1.1, 'H_10': -0.7, 'Ip_10': +0.2, 'In_10': -0.2,
                                                               'F_1': -2.2, 'G_1': -1.8, 'H_1': -1.2, 'Ip_1': +0.2, 'In_1': -0.2,},
    'parapets 50': {'hp_h': 0.050, 'r_h': 0.0, 'alpha': 90, 'F_10': -1.4, 'G_10': -0.9, 'H_10': -0.7, 'Ip_10': +0.2, 'In_10': -0.2,
                                                               'F_1': -2.0, 'G_1': -1.6, 'H_1': -1.2, 'Ip_1': +0.2, 'In_1': -0.2,},
    'parapets 100': {'hp_h': 0.10, 'r_h': 0.0, 'alpha': 90, 'F_10': -1.2, 'G_10': -0.8, 'H_10': -0.7, 'Ip_10': +0.2, 'In_10': -0.2,
                                                               'F_1': -1.8, 'G_1': -1.4, 'H_1': -1.2, 'Ip_1': +0.2, 'In_1': -0.2,},
    'curved eaves 5': {'hp_h': 0.00, 'r_h': 0.05, 'alpha': 90, 'F_10': -1.0, 'G_10': -1.2, 'H_10': -0.4, 'Ip_10': +0.2, 'In_10': -0.2,
                                                                'F_1': -1.5, 'G_1': -1.8, 'H_1': -0.4, 'Ip_1': +0.2, 'In_1': -0.2,},
    'curved eaves 10': {'hp_h': 0.00, 'r_h': 0.10, 'alpha': 90, 'F_10': -0.7, 'G_10': -0.8, 'H_10': -0.3, 'Ip_10': +0.2, 'In_10': -0.2,
                                                                'F_1': -1.2, 'G_1': -1.4, 'H_1': -0.3, 'Ip_1': +0.2, 'In_1': -0.2,},
    'curved eaves 20': {'hp_h': 0.00, 'r_h': 0.20, 'alpha': 90, 'F_10': -0.5, 'G_10': -0.5, 'H_10': -0.3, 'Ip_10': +0.2, 'In_10': -0.2,
                                                                'F_1': -0.8, 'G_1': -0.8, 'H_1': -0.3, 'Ip_1': +0.2, 'In_1': -0.2,},
    'mansard eaves 30': {'hp_h': 0.00, 'r_h': 0.00, 'alpha': 30, 'F_10': -1.0, 'G_10': -1.0, 'H_10': -0.3, 'Ip_10': +0.2, 'In_10': -0.2,
                                                                'F_1': -1.5, 'G_1': -1.5, 'H_1': -0.3, 'Ip_1': +0.2, 'In_1': -0.2,},
    'mansard eaves 45': {'hp_h': 0.00, 'r_h': 0.00, 'alpha': 45, 'F_10': -1.2, 'G_10': -1.3, 'H_10': -0.4, 'Ip_10': +0.2, 'In_10': -0.2,
                                                                'F_1': -1.8, 'G_1': -1.9, 'H_1': -0.4, 'Ip_1': +0.2, 'In_1': -0.2,},
    'mansard eaves 60': {'hp_h': 0.00, 'r_h': 0.00, 'alpha': 60, 'F_10': -1.3, 'G_10': -1.3, 'H_10': -0.5, 'Ip_10': +0.2, 'In_10': -0.2,
                                                                'F_1': -1.9, 'G_1': -1.9, 'H_1': -0.5, 'Ip_1': +0.2, 'In_1': -0.2,},
    'mansard eaves 90': {'hp_h': 0.00, 'r_h': 0.00, 'alpha': 90, 'F_10': -1.8, 'G_10': -1.2, 'H_10': -0.7, 'Ip_10': +0.2, 'In_10': -0.2,
                                                               'F_1': -2.5, 'G_1': -2.0, 'H_1': -1.2, 'Ip_1': +0.2, 'In_1': -0.2,},
}


def get_flat_roof_c_pe(zone:str, area:float, roof_type:str = 'sharp eaves', 
                       hp_h:float = 0, r_h:float = 0, alpha:float = 90, 
                       limit:str='negative')->float:
    """Returns external pressure coefficient for flat roofs based on values
    in Table 7.2 (taken from dictionary: 'FLAT_ROOF_COEFF_DICT')
    (Table 7.2 of EC1-1-4)

    Args:
        roof_type (str): 'sharp eaves' , 'parapets', 'curved eaves', 'mansard eaves'
        zone (str): 'F' to 'I', or 'Ip' or 'In' (since positive or negative pressure 
            can be returned, default for 'I' is negative.)
        area (float): _area of element(s) under consideration (sq.m.)
        hp_h (float, optional): Required for parapet cases: height of parapet (h_p) 
            over building height (h). Defaults to 0.
        r_h (float, optional): Required for curved eaves cases: radius of eaves (r) 
            over building height (h). Defaults to 0.
        alpha (float, optional): Required for mansard cases: angle of mansard. 
            Defaults to 90.
        limit (str, optional): 'positive' or 'negative' value - only relevant to 
            zone I. Defaults to 'negative'.

    Raises:
        ValueError: if roof type is not contained in table

    Returns:
        float: external pressure coefficient for flat roofs (Table 7.2 of EC1-1-4)
    
    Examples:
    >>> get_flat_roof_c_pe('G', 10.0)
    -1.2
    >>> get_flat_roof_c_pe('G', 10.0, 'sharp eaves')
    -1.2
    >>> get_flat_roof_c_pe('In', 1.0, 'sharp eaves')
    -0.2
    >>> get_flat_roof_c_pe('I', 1.0, 'sharp eaves', limit = 'positive')
    0.2
    >>> get_flat_roof_c_pe('G', 1.0, 'curved eaves', r_h = 0.1)
    -1.4
    >>> get_flat_roof_c_pe('H', 1.0, 'mansard eaves', alpha = 45)
    -0.4
    >>> get_flat_roof_c_pe('F', 10.0, 'sharp eaves' )
    -1.8
    """
    if zone == 'I':
        if limit == 'positive':
            zone = 'Ip' # positive
        else:
            zone = 'In' # negative
    if roof_type == 'sharp eaves':
        c_pe_dict = FLAT_ROOF_COEFF_DICT.get('sharp eaves', None)
        c_pe_10, c_pe_1 = [c_pe_dict.get(zone + area) for area in ('_10', '_1')]
        c_pe = list_interp(area, (1.0, 10.0), (c_pe_1, c_pe_10), extend=True)
    else:
        pattern = ''
        if roof_type.startswith(pattern := 'parapets'):
            name, value, scale = 'hp_h', hp_h, 1000
        elif roof_type.startswith(pattern := 'curved eaves'):
            name, value, scale = 'r_h', r_h, 100
        elif roof_type.startswith(pattern := 'mansard eaves'):
            name, value, scale = 'alpha', alpha, 1
        else:
            raise ValueError(f'Unknown roof type: {roof_type}')
        
        # Identify relevant dictionaries based on category (pattern), creating a new scaled integer key for each.
        c_pe_dicts = {int(scale * v.get(name, 0)): v for k, v in FLAT_ROOF_COEFF_DICT.items() if k.startswith(pattern)}
        
        # Get the results for the area extremes (c_pe_10 and c_pe_1 for the range of the variable (h_p/h, r/h & alpha))
        c_pe_10_range, c_pe_1_range = [[(k, c_pe_dict.get(zone + area)) for k, c_pe_dict in c_pe_dicts.items()] for area in ('_10', '_1')]
        c_pe_10, c_pe_1 = [list_interp(value * scale, *zip(*c_pe_range), extend=True) for c_pe_range in (c_pe_10_range, c_pe_1_range)]
        c_pe_min, c_pe_max = sorted([c_pe_1, c_pe_10])
        c_pe = max(c_pe_min, min(c_pe_max, c_pe_1 - (c_pe_1 - c_pe_10) * log10(area)))
    
    return c_pe


def main():
    d, b = 10.91, 6.173
    print(f'd = {d}, b = {b}, c_f0 = {get_c_f0(d, b):.3f}')
    d, b = 6.611, 6.173
    print(f'd = {d}, b = {b}, c_f0 = {get_c_f0(d, b):.3f}')
    # print(f'd = 6.611, b = 6.173, c_f0 = {get_c_f0(6.611, 6.173):.3f}')


if __name__ == "__main__":
    main()
    doctest.testmod()
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
import matplotlib.pyplot as plt
from pySeismic.utilities import list_interp  # noqa: E402
# from ..utilities import list_interp


def get_c_f0(d, b):
    """Force coefficient c_f,0 of rectangular sections with sharp 
    corners and without free end flow (Figure 7.23).
    
    The force coefficient c_f is calculated by multiplying c_f,0 by
    psi_r and psi_lambda.

    Args:
        d (float): width (dimension in direction of wind flow)
        b (float): breadth (dimension perpendicular to wind flow)

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


def plot_c_f0(d_b_dict):
    """Force coefficient c_f,0 of rectangular sections with sharp 
    corners and without free end flow (Figure 7.23).
    
    Note: the force coefficient c_f is calculated by multiplying c_f,0 by
    psi_r and psi_lambda.

    Args:
        d_b_dict (dict[str]): width (dimension in direction of wind flow)
            divided by breadth (dimension perpendicular to wind flow)

    Returns:
        list[float]: list of c_f0 values
    
    Example:
        >>> get_c_f0(9.2, 10.0)
        2.18
        >>> get_c_f0(9.0, 1.0)
        0.92
        >>> [get_c_f0(d, 1) for d in (0.05, 0.1, 0.2, 0.7, 1.0, 2.0, 5.0, 10.0, 50.0, 60.0)]
        [2.0, 2.0, 2.0, 2.4, 2.1, 1.65, 1.0, 0.9, 0.9, 0.9]
    """
    c_f0_list = []
    
    d_b = (0.1, 0.2, 0.6, 0.7, 1.0, 2.0, 5.0, 10.0, 50.0)
    c_f0 = (2.0, 2.0, 2.35, 2.4, 2.1, 1.65, 1.0, 0.9, 0.9)
    
    # Plot the line
    plt.plot(d_b, c_f0, label='')

    # Plot the dots
    for k, d_b_val in d_b_dict.items():
        c_f0_val = list_interp(d_b_val, d_b, c_f0, extend=True)
        c_f0_list.append(c_f0_val)
        plt.plot(d_b_val, c_f0_val, 'o', label=k)

    # Set x and y limits
    plt.xlim(0.1, 60)
    plt.ylim(0.0, 2.8)

    # Set x-axis to log scale
    plt.xscale('log')

    # Set x and y labels
    plt.xlabel('d / b')
    plt.ylabel('$c_{f,0}$')

    # Add gridlines
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Add specific lines
    [plt.axhline(y=yval, color='grey', linestyle='-', linewidth=0.5, zorder=0) 
     for yval in c_f0[1:-1]]
    [plt.axvline(x=xval, color='grey', linestyle='-', linewidth=0.5, zorder=0)
     for xval in d_b[1:-2]]

    # Add text at specific coordinates
    [plt.text(0.1, y+0.02, f'{y}', fontsize=7, ha='left') for y in c_f0[1:-1]]

    # Add a title
    plt.title('Figure 7.23 Force coefficients $c_{f,0}$')

    # Add a legend
    plt.legend()

    # Show the plot
    plt.show()

    return c_f0_list


def get_c_fz(b, d_tot, theta=0):
    return list_interp(b / d_tot, (0,14,22), (0.75, 0.15, 0.15), extend=True)


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
    
    The effect of area upon the output is taken from figure 7.2

    Args:
        roof_type (str): 'sharp eaves' , 'parapets', 'curved eaves', 'mansard eaves'
        zone (str): 'F' to 'I', or 'Ip' or 'In' (since positive or negative pressure 
            can be returned, default for 'I' is negative.)
        area (float): area of element(s) under consideration (sq.m.)
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
        # c_pe = list_interp(area, (1.0, 10.0), (c_pe_1, c_pe_10), extend=True)
        c_pe = max(c_pe_10, min(c_pe_1, c_pe_1 - (c_pe_1 - c_pe_10) * log10(area))) # from Figure 7.2
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
        c_pe = max(c_pe_min, min(c_pe_max, c_pe_1 - (c_pe_1 - c_pe_10) * log10(area))) # from Figure 7.2
        
    return c_pe


def main():
    d, b = 10.91, 6.173
    print(f'd = {d}, b = {b}, c_f0 = {get_c_f0(d, b):.3f}')
    d, b = 6.611, 6.173
    print(f'd = {d}, b = {b}, c_f0 = {get_c_f0(d, b):.3f}')
    # print(f'd = 6.611, b = 6.173, c_f0 = {get_c_f0(6.611, 6.173):.3f}')


if __name__ == "__main__":
    doctest.testmod()
    print('== Doctest Complete ==')
    main()
    
"""
A library for carrying out calculations to ASCE7-16

"This Standard provides requirements for 
general structural design and includes 
means for determining dead, live, soil, 
flood, wind, snow, rain, atmospheric 
ice, and earthquake loads, as well as 
their combinations, which are suitable 
for inclusion in building codes and 
other documents."

Previous editions are:
* ASCE 7-10
* ASCE 7-05
* ASCE 7-02
* ASCE 7-98
* ASCE 7-95
* ASCE 7-93
* ASCE 7-88
"""

from math import pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
# from itertools import accumulate

from pySeismic.utilities import list_interp

#### Seismic Hazard ####

# Functions
# T_code_get

# F_a and F_v taken from ASCE7-16, tables 11.4-1 & 11.4-2
F_a_headers = (0.25,0.5,0.75,1.0,1.25,1.5)
F_a_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
           'B': (0.9,0.9,0.9,0.9,0.9,0.9),
           'C': (1.3,1.3,1.2,1.2,1.2,1.2),
           'D': (1.6,1.4,1.2,1.1,1.0,1.0),
           'E': (2.4,1.7,1.3,None,None,None),
           'F': (None,None,None,None,None,None)}

F_v_headers = (0.1,0.2,0.3,0.4,0.5,0.6)
F_v_dict = {'A': (0.8,0.8,0.8,0.8,0.8,0.8),
           'B': (0.8,0.8,0.8,0.8,0.8,0.8),
           'C': (1.5,1.5,1.5,1.5,1.5,1.4),
           'D': (2.4,2.2,2.0,1.9,1.8,1.7),
           'E': (4.2,None,None,None,None,None),
           'F': (None,None,None,None,None,None)}


def site_coeffs(site_class, S_s, S_1, debug=False):
    """Returns F_a, F_v when provided with 
    site class (A-F) and S_s and S_1 values (ASCE 7-16)
    """
    F_a_list = F_a_dict[site_class]
    F_v_list = F_v_dict[site_class]
    if debug:
        print(f'site_class = {site_class}; S_s = {S_s}; S_1 = {S_1}; F_a = {F_a_list}; F_v = {F_v_list}')
    return list_interp(S_s, F_a_headers, F_a_list), list_interp(S_1, F_v_headers, F_v_list)




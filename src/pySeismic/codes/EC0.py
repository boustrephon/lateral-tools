"""
Calculations for Eurocode 1991-1-4 - Wind Loads
"""

import sys
from os.path import join as path_join, dirname, abspath
import doctest

# Add the parent directory of 'pySeismic' to sys.path
sys.path.append(abspath(path_join(dirname(__file__), '..', '..')))

# print('sys_path:', sys.path)

# from math import log10
# import matplotlib.pyplot as plt
# from pySeismic.utilities import list_interp
# from ..utilities import list_interp

# Table A1.1 (EN 1990, 2005)
# Imposed loads in building by category
EN1990_TABLE_A1_1_DICT = {
    'A' : {'description': 'domestic, residential areas', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.3},
    'B' : {'description': 'office areas', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.3},
    'C' : {'description': 'congregation areas', 'psi_0': 0.7, 'psi_1': 0.7, 'psi_2': 0.6},
    'D' : {'description': 'shopping areas', 'psi_0': 0.7, 'psi_1': 0.7, 'psi_2': 0.6},
    'E' : {'description': 'storage areas', 'psi_0': 1.0, 'psi_1': 0.9, 'psi_2': 0.8},
    'F' : {'description': 'traffic area, vehicle weight ≤ 30kN', 'psi_0': 0.7, 'psi_1': 0.7, 'psi_2': 0.6},
    'G' : {'description': 'traffic area, 30kN < vehicle weight ≤ 160kN', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.3},
    'H' : {'description': 'roofs', 'psi_0': 0.0, 'psi_1': 0.0, 'psi_2': 0.0},
    'FINS Snow': {'description': 'Finland, Iceland, Norway, Sweden', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.2},
    'High Snow': {'description': 'Other CEN sites > 1,000 m', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.2},
    'Low Snow': {'description': 'Other CEN sites <= 1,000 m', 'psi_0': 0.5, 'psi_1': 0.2, 'psi_2': 0.0},
    'Wind': {'description': 'shopping areas', 'psi_0': 0.6, 'psi_1': 0.2, 'psi_2': 0.0},
    'Temperature': {'description': 'shopping areas', 'psi_0': 0.6, 'psi_1': 0.5, 'psi_2': 0.0},
}

EN1990_UK_TABLE_A1_1_DICT = {
    'A' : {'description': 'domestic, residential areas', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.3},
    'B' : {'description': 'office areas', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.3},
    'C' : {'description': 'congregation areas', 'psi_0': 0.7, 'psi_1': 0.7, 'psi_2': 0.6},
    'D' : {'description': 'shopping areas', 'psi_0': 0.7, 'psi_1': 0.7, 'psi_2': 0.6},
    'E' : {'description': 'storage areas', 'psi_0': 1.0, 'psi_1': 0.9, 'psi_2': 0.8},
    'F' : {'description': 'traffic area, vehicle weight ≤ 30kN', 'psi_0': 0.7, 'psi_1': 0.7, 'psi_2': 0.6},
    'G' : {'description': 'traffic area, 30kN < vehicle weight ≤ 160kN', 'psi_0': 0.7, 'psi_1': 0.5, 'psi_2': 0.3},
    'H' : {'description': 'roofs', 'psi_0': 0.7, 'psi_1': 0.0, 'psi_2': 0.0},
    'High Snow': {'description': 'Other CEN sites > 1,000 m', 'psi_0': 0.5, 'psi_1': 0.2, 'psi_2': 0.0},
    'Low Snow': {'description': 'Other CEN sites <= 1,000 m', 'psi_0': 0.5, 'psi_1': 0.2, 'psi_2': 0.0},
    'Wind': {'description': 'shopping areas', 'psi_0': 0.5, 'psi_1': 0.2, 'psi_2': 0.0},
    'Temperature': {'description': 'shopping areas', 'psi_0': 0.6, 'psi_1': 0.5, 'psi_2': 0.0},
}


def get_psi(category, region=None):
    if region is None:
        A1_1_dict = EN1990_TABLE_A1_1_DICT
    elif region == 'UK':
        A1_1_dict = EN1990_UK_TABLE_A1_1_DICT
    else:
        if isinstance(region, dict):
            A1_1_dict = region
        else:
            A1_1_dict = {}
    return A1_1_dict.get(category)
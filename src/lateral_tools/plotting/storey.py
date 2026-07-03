"""
storey.py — Storey plotting functions.

Part of the lateral-tools package for seismic and lateral load engineering.
"""
from __future__ import annotations

"""
These are plotting functions for seismic engineering. The data needs to 
be generated elsewhere, for example from eqsigs.

"""

from collections.abc import Sequence, Iterator
from typing import TypeVar, Any, overload, Literal


from typing import Union, List, Tuple, Dict, Any
from math import log, pi, tau, sin, cos, tan, acos, asin, atan, atan2, log, log10  # noqa: F401, F811
import math
from sys import path as sys_path
from os.path import abspath, join as path_join, dirname

# from itertools import accumulate

import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.figure import Figure

import numpy as np

from scipy.interpolate import make_interp_spline, Akima1DInterpolator, PchipInterpolator

from matplotlib.ticker import EngFormatter

sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

from lateral_tools.analysis import storey_force_distribution, \
    storey_shear_moment, shear_coords


# Type aliases to describe the nested data structures
Curve = Sequence[float]                     # e.g. a list of periods or velocities
CurveCollection = Sequence[Curve]            # multiple curves (e.g. several period lists)
Point = tuple[float, float] | Sequence[float]  # a (T, Sv) pair
PointCollection = Sequence[Sequence[Point]]   # curves given as (T, Sv) pairs

# For the generic transpose helper
T = TypeVar('T')

List_of_lists = Union[
    List[List[float]], 
    List[Tuple[float]], 
    Tuple[List[float]],
                            ]
List_of_list_of_lists = Union[
    List[List_of_lists],
    Tuple[List_of_lists],
                            ]
List_of_list_of_list_of_lists = Union[
    List[List_of_list_of_lists],
    Tuple[List_of_list_of_lists],
                            ]

# Union[List[List[List[float]]], List[List[List[float]]], List[List[Tuple[float, float]]], List[Tuple[Tuple[float, float], ...]], Tuple[List[Tuple[float, float]]], Tuple[Tuple[Tuple[float, float], ...]]]
# Union[List[List[List[float]]], List[List[List[float]]], List[List[Tuple[float, float]]], List[Tuple[Tuple[float, float], ...]], Tuple[List[Tuple[float, float]]], Tuple[Tuple[Tuple[float, float], ...]]]



# ========================================================
# =================== Storey Plots =======================
# ========================================================

# from lateral_tools.analysis import storey_force_distribution, \
#     storey_shear_moment, shear_coords
# from matplotlib.ticker import EngFormatter


def plot_storey_VM(
        base_shears: Union[float, list, tuple], 
        storey_elevs: Union[list[float], tuple[float]], 
        storey_wts: Union[list[float], tuple[float]], 
        k:float = 1, 
        base_level: float = 0, 
        **kwargs: dict):
    """
    
    Args:
        base_shears: base shear(s)
        storey_elevs: storey elevations
        storey_wts:  
        k: (default = 1)
        base_level: (default = 0), 
        **kwargs
    
    Returns:
        fig: Storey shear and moment plot
    """

    # x_limits = kwargs.pop('x_limits', (0, None))
    # y_limits = kwargs.pop('y_limits', (0, None))
    plot_title = kwargs.pop('plot_title', 'Building Shears & Moments')
    figsize = kwargs.pop('figsize', (6,8))
    force_units = kwargs.pop('force_units', 'kN')
    length_units = kwargs.pop('length_units', 'm')
    # debug = kwargs.pop('debug', False)

    force_scale_dict = {'N':1.0, 'kN':1000.0, 'MN':1.0E6, 'GN':1.0E9}
    length_scale_dict = {'mm':0.001, 'cm':0.01, 'm':1.0, 'km':1000.0}
    V_scale = force_scale_dict.get(force_units)
    M_scale = length_scale_dict.get(length_units) * V_scale

    if list(storey_elevs) == sorted(storey_elevs):
        ascending = True
    elif list(storey_elevs) == sorted(storey_elevs, reverse=True):
        ascending = False
    else:
        raise ValueError('Storey elevations are not monotonic')
    
    if ascending:
        storey_elevs = storey_elevs[::-1]
        storey_wts = storey_wts[::-1]

    storey_heights = [z2 - z1 for z2, z1 in zip(
        storey_elevs, 
        list(storey_elevs[1:]) + [0]
        )]

    V = base_shears if isinstance(base_shears, (list, tuple)) else [base_shears]

    storey_loads = [storey_force_distribution(
        vv, storey_elevs, storey_wts, k) for vv in V]
    
    storey_shears_and_moments = [storey_shear_moment(
        forces, storey_heights, ascending) for forces in storey_loads]
    storey_shears, storey_moments = zip(*storey_shears_and_moments)

    z_V = [shear_coords(
        storey_elevs, V, base_level) for V in storey_shears]

    fig, [ax_V, ax_M] = plt.subplots(nrows=1, ncols=2, figsize=figsize)
    fig.suptitle(plot_title)

    ax_V.invert_xaxis()
    for z_VV in z_V:
        z_range, V_range = zip(*z_VV)
        V_range = [V_scale * V for V in V_range]
        ax_V.plot([0] + list(V_range), list(z_range[:1]) + list(z_range))
    ax_V.set_xlim((None, 0))
    ax_V.set_title('Storey Shear')
    ax_V.set_xlabel('Shear')
    formatterV = EngFormatter(unit='N')
    ax_V.xaxis.set_major_formatter(formatterV)

    for M_range in storey_moments:
        M_range = [M_scale * M for M in M_range]
        ax_M.plot([0] + M_range, list(storey_elevs) + [0])
    ax_M.set_xlim((0, None))
    ax_M.set_title('Storey Moment')
    ax_M.set_xlabel('Moment')
    formatterM = EngFormatter(unit='Nm')
    ax_M.xaxis.set_major_formatter(formatterM)
    
    for ax in (ax_V, ax_M):
        ax.set_ylim((0, None))
        ax.minorticks_on()
        ax.grid(True, which='major', linestyle='-', linewidth=0.5)
        ax.grid(True, which='minor', linestyle='--', linewidth=0.25)
        ax.tick_params(axis='x', rotation=45)


    fig.set_tight_layout(True)
    return fig


def plot_storey_VM2(
        storey_elevs: Sequence[float], 
        storey_shears: Sequence[Sequence[float]], 
        storey_moments: Sequence[Sequence[float]], 
        base_level: float = 0, 
        **kwargs: Any):
    """
    
    Args:
        storey_elevs: storey elevations
        storey_shears:  list of list of storey shears
        storey_moments: list of list of storey moments
        base_level: (default = 0), 
        **kwargs
    
        Note that the storey shears and moments should be in the 
        same order as the storey elevations, and the moments should 
        include the base moment but not the top moment, which
        is assumed to be zero. In this way the length of the lists 
        of shears should match the length of the lists of moments
        
    Returns:
        fig: Storey shear and moment plot
    """

    # x_limits = kwargs.pop('x_limits', (0, None))
    # y_limits = kwargs.pop('y_limits', (0, None))
    plot_title = kwargs.pop('plot_title', 'Building Shears & Moments')
    figsize = kwargs.pop('figsize', (6,8))
    force_units = kwargs.pop('force_units', 'kN')
    length_units = kwargs.pop('length_units', 'm')
    # moment_units = (force_units, length_units)
    debug = kwargs.pop('debug', False)

    force_scale_dict = {'N':1.0, 'kN':1000.0, 'MN':1.0E6, 'GN':1.0E9}
    length_scale_dict = {'mm':0.001, 'cm':0.01, 'm':1.0, 'km':1000.0}

    V_scale = force_scale_dict[force_units]
    M_scale = length_scale_dict[length_units] * V_scale

    # Identify whether the storey elevations are ascending or descending
    if list(storey_elevs) == sorted(storey_elevs):
        ascending = True
    elif list(storey_elevs) == sorted(storey_elevs, reverse=True):
        ascending = False
    else:
        raise ValueError('Storey elevations are not monotonic')

    if debug:
        print(f'Storey levels ascending: {ascending}')

    # if elevations are descending, reorganise as ascending
    if not ascending:
        storey_elevs = storey_elevs[::-1]
        storey_shears = storey_shears[::-1]
        storey_moments = storey_moments[::-1]
        # storey_wts = storey_wts[::-1]

    # Calculate the shear and moment for each storey from top down
    z_V = [shear_coords(
        storey_elevs, V, base_level) for V in storey_shears]
    
    # z_M = [[(base_level, 0.0)] + [(elev, M) for elev, M in zip(storey_elevs, M_range)] # + [()] 
    #     for M_range in storey_moments]
    z_M = [[(elev, M) for elev, M in zip(storey_elevs, M_range)] # + [()] 
        for M_range in storey_moments]

    fig, [ax_V, ax_M] = plt.subplots(nrows=1, ncols=2, figsize=figsize)
    fig.suptitle(plot_title)

    ax_V.invert_xaxis()
    for z_VV in z_V:
        z_range, V_range = zip(*z_VV)
        V_range = [V_scale * V for V in V_range]
        if debug:
            print(f'V_range: {V_range}')
            print(f'z_range: {z_range}')
        ax_V.plot([0] + list(V_range), list(z_range[:1]) + list(z_range))
    ax_V.set_xlim((None, 0))
    ax_V.set_title('Storey Shear')
    ax_V.set_xlabel('Shear')
    formatterV = EngFormatter(unit='N')
    ax_V.xaxis.set_major_formatter(formatterV)

    if debug:
        print(f'length of z: {len(storey_elevs)}')
        print('length of z_M:')
        _ = [print(len(z_MM), [len(z) for z in z_MM]) for z_MM in z_M]
    for z_MM in z_M:
        z_range, M_range = zip(*z_MM)
        M_range = [M_scale * M for M in M_range]
        # M_range = M_range[1:] + M_range[:1]
        if debug:
            print(f'z_range: {z_range}')
            print(f'M_range: {M_range}')
        ax_M.plot([0] + list(M_range) + [0], [base_level]*2 + list(z_range))
    ax_M.set_xlim((0, None))
    ax_M.set_title('Storey Moment')
    ax_M.set_xlabel('Moment')
    formatterM = EngFormatter(unit='Nm')
    ax_M.xaxis.set_major_formatter(formatterM)

    for ax in (ax_V, ax_M):
        ax.set_ylim((0, None))
        ax.minorticks_on()
        ax.grid(True, which='major', linestyle='-', linewidth=0.5)
        ax.grid(True, which='minor', linestyle='--', linewidth=0.25)
        ax.tick_params(axis='x', rotation=45)


    fig.set_tight_layout(True)
    return fig


# ===============================================================
# EQ Plotting functions
# ===============================================================
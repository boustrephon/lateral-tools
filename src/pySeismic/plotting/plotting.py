"""
These are plotting functions for seismic engineering. The data needs to 
be generated elsewhere, for example from eqsigs.

"""
from __future__ import annotations

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

# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

from pySeismic.utilities import storey_force_distribution, \
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

# from pySeismic.utilities import storey_force_distribution, \
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

def reorganise_plot_data(
    list_of_list_of_lists: Sequence[Sequence[Sequence[T]]]
) -> Iterator[tuple[Sequence[T], ...]]:
    """
    Transpose a three‑level nested sequence.

    Example:
        input:  [ [(a1,b1), (a2,b2)], [(c1,d1), (c2,d2)] ]
        output: ( (a1,a2), (c1,c2) ), ( (b1,b2), (d1,d2) )
    """
    # Inner transpose: turn each curve's points into (first_coords, second_coords, ...)
    transposed_curves = [zip(*curve) for curve in list_of_list_of_lists]
    # Outer transpose: group all first_coords, all second_coords, etc.
    return zip(*transposed_curves)


@overload
def tripartite_plot(
    T_ranges: Curve,
    S_v_ranges: CurveCollection,
    series_names: Sequence[str] = ...,
    plot_pts: Sequence[Point] = ...,
    debug: bool = ...,
    figsize: tuple[float, float] = ...,
    suptitle: str = ...,
    spec_type: str = ...,
    units: str | None = ...,
    minor_lines: Sequence[int] = ...,
    **kwargs: Any,
) -> Figure:
    """Called when T_ranges is a single period list (same periods for all curves)."""
    ...


@overload
def tripartite_plot(
    T_ranges: CurveCollection,
    S_v_ranges: CurveCollection,
    series_names: Sequence[str] = ...,
    plot_pts: Sequence[Point] = ...,
    debug: bool = ...,
    figsize: tuple[float, float] = ...,
    suptitle: str = ...,
    spec_type: str = ...,
    units: str | None = ...,
    minor_lines: Sequence[int] = ...,
    **kwargs: Any,
) -> Figure:
    """Called when T_ranges is already a collection of period lists."""
    ...


@overload
def tripartite_plot(
    T_ranges: Literal[[]],  # empty list triggers the reorganise branch
    S_v_ranges: PointCollection,
    series_names: Sequence[str] = ...,
    plot_pts: Sequence[Point] = ...,
    debug: bool = ...,
    figsize: tuple[float, float] = ...,
    suptitle: str = ...,
    spec_type: str = ...,
    units: str | None = ...,
    minor_lines: Sequence[int] = ...,
    **kwargs: Any,
) -> Figure:
    """Called when periods are embedded in S_v_ranges as (T, Sv) pairs."""
    ...


def tripartite_plot(
    T_ranges: Curve | CurveCollection | list,
    S_v_ranges: CurveCollection | PointCollection,
    series_names: Sequence[str] = (),
    plot_pts: Sequence[Point] = (),
    debug: bool = False,
    figsize: tuple[float, float] = (10, 8),
    suptitle: str = '',
    spec_type: str = 'V',
    units: str | None = None,
    minor_lines: Sequence[int] = (1, 2, 3, 4, 5, 6, 7, 8, 9),
    **kwargs: Any,
) -> Figure:
    """
    Plots tripartite (log‑log) spectral data.

    Two main calling patterns:
    1. Provide periods in `T_ranges` (either a single list or a list of lists)
       and velocities in `S_v_ranges` (a list of lists).
    2. Provide `T_ranges = []` and give (T, Sv) pairs inside `S_v_ranges`
       (a list of lists of (T, Sv) points). `reorganise_plot_data` will split them.

    Parameters
    ----------
    T_ranges : single curve or collection of curves (periods)
        If empty list, periods are taken from `S_v_ranges` (pattern 2).
    S_v_ranges : collection of curves (velocities) or collection of (T, Sv) pairs
        Spectral values to plot.
    series_names : names for each curve (used in legend)
    plot_pts : extra (T, Sv) points to mark
    debug : print intermediate values
    figsize : figure size in inches
    suptitle : overall figure title
    spec_type : 'V' for velocity (default), 'A' to convert acceleration to velocity
    units : not used (kept for compatibility)
    minor_lines : which minor log‑cycle lines to draw (1…9)
    **kwargs : passed to `ax.loglog`
    """
    # ===== 1. Handle the embedded case =====
    if T_ranges == []:
        # S_v_ranges is PointCollection
        transposed = reorganise_plot_data(S_v_ranges)      # type: ignore[arg-type]
        T_tuples, Sv_tuples = transposed
        # Convert tuples of sequences to lists of lists
        T_ranges = [list(ts) for ts in T_tuples]
        S_v_ranges = [list(vs) for vs in Sv_tuples]

    # ===== 2. Normalise T_ranges to a list of lists =====
    # At this point T_ranges is either Curve or CurveCollection.
    # Convert a single curve into a list containing that curve repeated.
    if not T_ranges:
        raise ValueError("T_ranges cannot be empty after embedded case")
    first = T_ranges[0]
    if isinstance(first, (int, float)):
        # Single curve -> repeat for each velocity curve
        T_ranges = [T_ranges] * len(S_v_ranges)          # type: ignore[list-item]
    elif len(T_ranges) == 1 and isinstance(first, (tuple, list)):
        # One curve given as a list inside a list – repeat it
        T_ranges = list(T_ranges) * len(S_v_ranges)      # type: ignore[list-item]
    # Now T_ranges is guaranteed to be a list of lists (CurveCollection)

    # ===== 3. Acceleration → velocity conversion =====
    if spec_type == 'A':
        tau = 2.0 * math.pi
        def rspec_ag_to_v(T: float, a: float) -> float:
            return a * 9.81 * T / tau
        # Both T_ranges and S_v_ranges are now CurveCollection
        S_v_ranges = [
            [rspec_ag_to_v(t, sa) for t, sa in zip(T_range, Sa_range)]
            for T_range, Sa_range in zip(T_ranges, S_v_ranges)
        ]

    # ===== 4. Validate lengths =====
    if not all(len(sv) == len(tr) for sv, tr in zip(S_v_ranges, T_ranges)):
        raise ValueError(
            f"Length mismatch: T lengths {[len(tr) for tr in T_ranges]}, "
            f"Sv lengths {[len(sv) for sv in S_v_ranges]}"
        )

    # ===== Create the plot =====
    if debug:
        print(f"Figure size: {figsize}")
        print(f"Super‑title: {suptitle}")

    fig, ax = plt.subplots(figsize=figsize)
    if suptitle:
        fig.suptitle(suptitle)
    fig.tight_layout()
    ax.grid(True, "both")

    # Plot the main curves
    for T_range, S_v_range in zip(T_ranges, S_v_ranges):
        ax.loglog(T_range, S_v_range, **kwargs)

    # ===== Tripartite grid lines =====
    # (the rest of the function is unchanged from the previous refactoring)
    T_min, T_max = ax.get_xlim()
    v_min, v_max = ax.get_ylim()

    a_min = 2.0 * math.pi / T_max * v_min
    a_max = 2.0 * math.pi / T_min * v_max
    d_min = 0.5 / math.pi * T_min * v_min
    d_max = 0.5 / math.pi * T_max * v_max

    if debug:
        print(f"T: {T_min:.3g} - {T_max:.3g}, Sv: {v_min:.3g} - {v_max:.3g}")
        print(f"Sa: {a_min:.3g} - {a_max:.3g}, Sd: {d_min:.3g} - {d_max:.3g}")

    log_T_min, log_T_max = math.log10(T_min), math.log10(T_max)
    log_v_min, log_v_max = math.log10(v_min), math.log10(v_max)
    log_a_min, log_a_max = math.log10(a_min), math.log10(a_max)
    log_d_min, log_d_max = math.log10(d_min), math.log10(d_max)

    a_fac, a_unit = 1.0, "m/s²"
    a_list = [
        m * 10.0**n
        for n in range(int(log_a_min) - 2, int(log_a_max) + 1)
        for m in minor_lines
        if a_min < m * 10.0**n * a_fac < a_max
    ]
    a_str = [f"{a * a_fac:5.1g},{a_unit}" for a in a_list]
    log_a_vals = [math.log10(a) for a in a_list]

    d_fac, d_unit = 1.0, "m"
    d_list = [
        m * 10.0**n
        for n in range(int(log_d_min), int(log_d_max) + 1)
        for m in minor_lines
        if d_min < m * 10.0**n * d_fac < d_max
    ]
    d_str = [f"{d * d_fac:5.1g},{d_unit}" for d in d_list]
    log_d_vals = [math.log10(d) for d in d_list]

    a_linestyle = ["-" if abs(math.log10(a) - int(math.log10(a))) < 1e-12 else ":" for a in a_list]
    d_linestyle = ["-" if abs(math.log10(d) - int(math.log10(d))) < 1e-12 else ":" for d in d_list]

    log_ta = [(log_a - log_a_min) / (log_a_max - log_a_min) for log_a in log_a_vals]
    log_td = [(log_d - log_d_min) / (log_d_max - log_d_min) for log_d in log_d_vals]

    log_va = [log_v_min + t * (log_v_max - log_v_min) for t in log_ta]
    log_Ta = [log_T_max + t * (log_T_min - log_T_max) for t in log_ta]

    log_vd = [log_v_min + t * (log_v_max - log_v_min) for t in log_td]
    log_Td = [log_T_min + t * (log_T_max - log_T_min) for t in log_td]

    s_a1 = [max(log_v_min - v, log_T_min - T) for v, T in zip(log_va, log_Ta)]
    s_a2 = [min(log_v_max - v, log_T_max - T) for v, T in zip(log_va, log_Ta)]

    s_d1 = [min(v - log_v_min, log_T_max - T) for v, T in zip(log_vd, log_Td)]
    s_d2 = [max(v - log_v_max, log_T_min - T) for v, T in zip(log_vd, log_Td)]

    a_lines = [
        lines.Line2D(
            (10 ** (T + s1), 10 ** (T + s2)),
            (10 ** (v + s1), 10 ** (v + s2)),
            transform=ax.transData,
            color="grey",
            ls=ls,
            lw=0.5,
            zorder=0,
        )
        for v, T, s1, s2, ls in zip(log_va, log_Ta, s_a1, s_a2, a_linestyle)
    ]

    d_lines = [
        lines.Line2D(
            (10 ** (T + s1), 10 ** (T + s2)),
            (10 ** (v - s1), 10 ** (v - s2)),
            transform=ax.transData,
            color="grey",
            ls=ls,
            lw=0.5,
            zorder=0,
        )
        for v, T, s1, s2, ls in zip(log_vd, log_Td, s_d1, s_d2, d_linestyle)
    ]

    # Find the longest solid line for placing the axis labels
    _, i_a = max(
        ((s2 - s1), i)
        for i, (s1, s2, ls) in enumerate(zip(s_a1, s_a2, a_linestyle))
        if ls == "-"
    )
    _, i_d = max(
        ((s1 - s2), i)
        for i, (s1, s2, ls) in enumerate(zip(s_d1, s_d2, d_linestyle))
        if ls == "-"
    )

    # Origin lines (magenta)
    ax_lines = [
        lines.Line2D(
            (10 ** (log_Ta[i] + s_a1[i]), 10 ** (log_Ta[i] + s_a2[i])),
            (10 ** (log_va[i] + s_a1[i]), 10 ** (log_va[i] + s_a2[i])),
            transform=ax.transData,
            color="magenta",
            ls="-",
            lw=1.0,
            zorder=0,
        )
        for i in (i_a, i_d)
    ]

    for ln in a_lines + d_lines + ax_lines:
        ax.add_line(ln)

    def line_angle(line: lines.Line2D) -> float:
        xd, yd = zip(*[ax.transData.transform(p) for p in line.get_xydata()])
        dx, dy = xd[1] - xd[0], yd[1] - yd[0]
        return math.degrees(math.atan2(dy, dx))

    angle_a = line_angle(ax_lines[0])
    angle_d = line_angle(ax_lines[1])

    for T_log, v_log, txt in zip(log_Ta, log_va, a_str):
        ax.text(10 ** T_log, 10 ** v_log, f"${txt}$",
                transform=ax.transData, rotation=angle_d, rotation_mode="anchor",
                alpha=0.5, zorder=0.1)

    for T_log, v_log, txt in zip(log_Td, log_vd, d_str):
        ax.text(10 ** T_log, 10 ** v_log, f"${txt}$",
                transform=ax.transData, rotation=angle_a, rotation_mode="anchor",
                alpha=0.5, zorder=0.1)

    if series_names:
        if len(series_names) == len(S_v_ranges):
            ax.legend(series_names)
        else:
            ax.legend([f"Spectrum {i+1:02d}" for i in range(len(S_v_ranges))])

    for T, v in plot_pts:
        ax.plot(T, v, marker="o")

    return fig


def reorganise_plot_data_old(
        list_of_list_of_lists: Union[
            List[List[List[float]]], 
            List[List[Tuple[float, float]]], 
            List[Tuple[Tuple[float, float], ...]], 
            Tuple[List[Tuple[float, float]]], 
            Tuple[Tuple[Tuple[float, float], ...]]
            ]
        ) -> Any:
    """Transposes lists of lists of vectors into lists
    of tuples of like coordinates.
    
    This is a process of transposition. For example, the
        ((a1,b1),(a2,b2),(a3,b3), ...),
                    ((c1,d1),(c2,d2),(c3,d3), ...)
        ->
        ((a1,a2,a3, ...),(c1,c2,c3, ...)),
                    ((b1,b2,b3, ...),(d1,d2,d3, ...))
    This is particularly useful where data is in the form
    of (x,y) pairs and you have multiple sets of these
    that you want to plot, but the plot function is looking
    for lists of x and lists of y and you want to plot multiple
    sets of these.
    """
    transposed_data = zip(
    *[new_list 
      for new_list in [zip(*each_list) 
                       for each_list 
                       in list_of_list_of_lists]]
    )
    return transposed_data



def tripartite_plot_old(
    T_ranges: Union[List[float], 
                    List[List[float]], 
                    Tuple[List[float]], 
                    Tuple[List[float], ...]], 
    S_v_ranges: Union[List_of_list_of_lists,
                      List_of_list_of_list_of_lists], 
    series_names: List[str] = [],
    plot_pts: List[float] = [], 
    debug: bool = False, 
    figsize: Tuple[float,float] = (10,8), 
    suptitle: str ='', 
    spec_type: str = 'V',
    units: str|None = None,
    minor_lines: List[int] = [1, 2, 3, 4, 5, 6, 7, 8, 9],
    **kwargs,
    ):
    """Plots tripartite data
        minorlines is a list of the minor log lines to plot (default [1,2,3,4,5,6,7,8,9])
    
    The following documentation will be helpful which provides information on the 
    "Artist" elements that are drawn on the plot (namely the diagonal axis lines) and the
    transformations between data and display coordinates and between plain and log coords:
    * https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.loglog.html
    * https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.loglog.html
    * https://matplotlib.org/stable/tutorials/advanced/transforms_tutorial.html
    * https://matplotlib.org/stable/api/artist_api.html#matplotlib.artist.Artist
    * https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html

    Args:
        T_ranges: collection of time values (each as a collection of 
            floats) if there is only one item, it is assumed that the 
            range is the same for all plots. If an empty list is passed
            to the function, it assumes that all the data is in the
            PS_a_ranges function as T,S pairs.
        S_v_ranges: collection of spectral velocities (each as a collection of floats)
        series_names (list[str]): one name for each S_v_range
        plot_pts (list[tuple,list]): (T, S_v) points for plotting 
        debug: if True then some intermediate results are printed to terminal
        figsize: size of the plot (width, height)  in inches
        suptitle: main title for the figure
        spec_type: 
        kwargs: parameters for passing through to `plot`
    Returns:
        Figure
    """
    # One of the challenges is the need to switch between plain and log10 values
    # 
    if T_ranges == []:
        T_ranges, S_v_ranges = reorganise_plot_data_old(S_v_ranges)
    if spec_type == 'A':
        def rspec_ag_to_v(T, a):
            return a * 9.81 * T / tau
        S_v_ranges = [[rspec_ag_to_v(t, sa) 
            for t, sa in zip(T_range, Sa_range)] 
            for T_range, Sa_range in zip(T_ranges, S_v_ranges)]



    # Standardise period and check that numbers align
    if isinstance(T_ranges[0],(int, float)):
        T_ranges = len(S_v_ranges) * [T_ranges]
    elif len(T_ranges) == 1 and isinstance(T_ranges[0],(tuple, list)):
        T_ranges = len(S_v_ranges) * T_ranges
    elif len(T_ranges) == len(S_v_ranges):
        pass
    else:
        pass

    if not (all(len(S_v_range) == len(T_range) for S_v_range, T_range in zip(S_v_ranges, T_ranges))):
        T_ranges_len = [len(T_range) for T_range in T_ranges]
        S_v_ranges_len = [len(S_v_range) for S_v_range in S_v_ranges]
        err_msg = (f'Lengths of curves do not match: \n\tlen(T_ranges) = {T_ranges_len} \n\tlen(S_v_ranges) = {S_v_ranges_len}')
        raise ValueError(err_msg)

    # calculate the spectral pseudo velocities (plain) - not needed
    # S_v_ranges = [[9.81 * 0.5 / pi * T * S_a for T, S_a in zip(T_ranges, S_a_range)] for S_a_range in S_a_ranges]
    # S_a_ranges = [[2 * pi / T * S_v / 9.81 for T, S_v in zip(T_ranges, S_v_range)] for S_v_range in S_v_ranges]
    

    # === Set up plot ==========
    if debug:
        print(f'Figure size is {figsize}')
        print(f'Super-title is {suptitle}')
    
    # Set up the basic plot
    fig0, ax0 = plt.subplots(figsize=figsize)
    if len(suptitle) > 0:
        fig0.suptitle(suptitle)
    fig0.tight_layout()
    ax0.grid(True, 'both')
    # fig = plt.gcf(); fig.set_size_inches(16,8)
    

    # Plot spectral velocity on log-log plot  (log-log space)
    [ax0.loglog(T_range, S_v_range, **kwargs) 
        for T_range, S_v_range in zip(T_ranges, S_v_ranges)]
    
    # ax0.loglog(T_code_s, S_v_g, marker = 'o') # marking the results for the building period
    # === End set up plot ===


    # rspec_v_to_a = lambda v, T: v * (2.0 * pi) / T
    # rspec_a_to_v = lambda a, T: a * T / (2.0 * pi)
    # rspec_v_to_d = lambda v, T: v * T / (2.0 * pi)
    # rspec_d_to_v = lambda d, T: d * (2.0 * pi) / T


    # === Tripartite Lines ========
    # Get limits (max & min) of Sv & T  (** no log space **)
    # - These are the four corners that are used for setting extent of 
    #   lines required for acceleration and displacement lines
    T_min, T_max = fig0.axes[0].get_xlim() #; print(x0)
    v_min, v_max = fig0.axes[0].get_ylim() #; print(y0)
    a_min, a_max = [2.0 * pi / T * v for v, T in ((v_min, T_max),(v_max, T_min))]
    d_min, d_max = [0.5 / pi * T * v for v, T in ((v_min, T_min),(v_max, T_max))]
    if debug:
        print(f'T & Sv limits\n  T_min = {T_min}; T_max = {T_max}; v_min = {v_min}; v_max = {v_max}')
        print(f'Sa & Sd limits\n  a_min = {a_min}; a_max = {a_max}; d_min = {d_min}; d_max = {d_max}')
    
    # Generate log values of T, Sv, Sa & Sd limits  (** log-log space **)
    log_T_min, log_T_max = log10(T_min), log10(T_max)
    log_v_min, log_v_max = log10(v_min), log10(v_max)
    log_a_min, log_a_max = log10(a_min), log10(a_max)
    log_d_min, log_d_max = log10(d_min), log10(d_max)
    

    #== Acceleration grid lines =========
    # Calculate a-values and d-values on the log-log line between a_min & a_max and d_min & d_max respectively
    if debug:
        print(f'Which ticks shown - {minor_lines}')
    
    # Scaling factors (currently set to 1.0)
    a_fac, a_unit = (1.0, 'm/s^2')
    # a_fac, a_unit = (9.81, 'g') # does not work!
    d_fac, d_unit = (1.0, 'm')
    
    # generate the list of acceleration line values from list of values provided they are within the limits
    a_list = [(m * 10**float(n)) for n in range(int(log10(a_min))-2, int(log10(a_max)+1)) 
            for m in minor_lines 
            if (m * 10**float(n) * a_fac > a_min) and (m * 10**float(n) * a_fac < a_max)]
    # precision
    a_str = [f'{a:5.1g},{a_unit}' for a in a_list]
    a_values = [a * a_fac for a in a_list]
    log_a_values = [log10(a) for a in a_values]
    # [print(a,b,s) for a,b,s in zip(a_list, a_values, a_str)]
    

    #== Displacement grid lines =========
    # generate the list of displacement line values provided they are within the limits
    d_list = [(m * 10**float(n)) for n in range(int(log10(d_min)), int(log10(d_max))+1) 
            for m in minor_lines 
            if (m * 10**float(n) * d_fac > d_min) and (m * 10**float(n) / d_fac < d_max)]
    d_str = [f'{d:5.1g},{d_unit}' for d in d_list]
    d_values = [d * d_fac for d in d_list]
    # print('d_values (min & max)', d_values[0], d_values[-1])
    log_d_values = [log10(d) for d in d_values]

    # Line-style definitions
    # - '-' means solid lines
    #   ':' means dotted lines
    # change the following lines
    # a_linestyle = ['-' if abs(a + log10(a_fac) - int(a + log10(a_fac))) < 0.001 else ':' for a in log_a_values]
    a_linestyle = ['-' if abs(log10(a) - int(log10(a))) < 0.001 else ':' for a in a_list]
    # d_linestyle = ['-' if abs(d - int(d)) < 0.001 else ':' for d in log_d_values]
    d_linestyle = ['-' if abs(log10(d) - int(log10(d))) < 0.001 else ':' for d in d_list]

    # T-values corresponding to a-lines
    log_ta_values = [(log_a - log_a_min)/(log_a_max - log_a_min) for log_a in log_a_values]
    # print('ta_min, ta_max', ta_values[0], ta_values[-1])
    
    # T-values corresponding to d-lines
    log_td_values = [(log_d - log_d_min)/(log_d_max - log_d_min) for log_d in log_d_values]
    # print('td_min, td_max', td_values[0], td_values[-1])

    # v-values & T-values corresponding to a-lines
    # acceleration
    log_va_values = [log_v_min + log_t * (log_v_max - log_v_min) for log_t in log_ta_values]
    # print('va_min, va_max', 10**va_values[0], 10**va_values[-1])
    log_Ta_values = [log_T_max + log_t * (log_T_min - log_T_max) for log_t in log_ta_values]
    # print('Ta_min, Ta_max', 10**Ta_values[0], 10**Ta_values[-1])
    
    # 
    log_vd_values = [log_v_min + log_t * (log_v_max - log_v_min) for log_t in log_td_values]
    # print('vd_min, vd_max', 10**vd_values[0], 10**vd_values[-1])
    log_Td_values = [log_T_min + log_t * (log_T_max - log_T_min) for log_t in log_td_values]
    # print('Td_min, Td_max', 10**Td_values[0], 10**Td_values[-1])

    # s1
    s_a2 = [min(log_v_max - v, log_T_max - T) 
            for v, T in zip(log_va_values, log_Ta_values)]
    s_a1 = [max(log_v_min - v, log_T_min - T) 
            for v, T in zip(log_va_values, log_Ta_values)]
    # print('s_a1[0], s_a2[0]', s_a1[0], s_a2[0])

    # identifying the index of the longest acceleration line for placing the text
    _, i_a = max((v2-v1,i) for i, (v1, v2, als) in enumerate(zip(s_a1, s_a2, a_linestyle)) if als == '-')
    # print('m_a, i_a', m_a, i_a)

    # s1
    s_d2 = [max(v - log_v_max, log_T_min - T) 
            for v, T in zip(log_vd_values, log_Td_values)]
    s_d1 = [min(v - log_v_min, log_T_max - T) 
            for v, T in zip(log_vd_values, log_Td_values)]
    # print('s_d1[0], s_d2[0]', s_d1[0], s_d2[0])

    # identifying the index of the longest displacement line
    _, i_d = max((v1 - v2, i) for i, (v1, v2, dls) in enumerate(zip(s_d1, s_d2, d_linestyle)) if dls == '-')
    # print('m_d, i_d', m_d, i_d)

    # calculate the coordinates of the ends of the acceleration lines (** non-log space **)
    a_coords = [((10**(T+s1), 10**(T+s2)), (10**(v+s1), 10**(v+s2))) for v, T, s1, s2 in zip(log_va_values, log_Ta_values, s_a1, s_a2)]
    # define the acceleration lines (** non-log space **)
    a_lines = [lines.Line2D(P1, P2, transform=ax0.transData, axes=ax0, color='grey', ls = ls, lw = 0.5, zorder=0) for (P1, P2), ls in zip(a_coords, a_linestyle)]

    # calculate the coordinates of the ends of the displacement lines (** non-log space **)
    d_coords = [((10**(T+s1), 10**(T+s2)), (10**(v-s1), 10**(v-s2))) for v, T, s1, s2 in zip(log_vd_values, log_Td_values, s_d1, s_d2)]
    # define the displacement lines (** non-log space **)
    d_lines = [lines.Line2D(P1, P2, transform=ax0.transData, axes=ax0, color='grey', ls = ls, lw = 0.5, zorder=0) for (P1, P2), ls in zip(d_coords, d_linestyle)]

    # print('a_coords', a_coords[i_a])
    # print('d_coords', d_coords[i_d])

    # Define 'origins' lines in magenta
    ax_lines = [lines.Line2D(P1, P2, transform=ax0.transData, axes=ax0, color='magenta', ls = '-', lw = 1.0, zorder=0) for (P1, P2) in (a_coords[i_a], d_coords[i_d])]

    # print('\nOrigin Lines')
    # [print(P1, P2) for P1, P2 in (a_coords[i_a], d_coords[i_d])]
    # print([(P1, P2) for P1, P2 in (a_coords[i_a], d_coords[i_d])])
    # print('\n')
    
    # Calculate gradient of origin lines in display units
    (((xd1, xd2), (yd1, yd2)), ((xa2, xa1), (ya2, ya1))) = [(P1, P2) for P1, P2 in (a_coords[i_a], d_coords[i_d])]
    xd1_display, yd1_display = ax0.transData.transform((xd1, yd1))
    xd2_display, yd2_display = ax0.transData.transform((xd2, yd2))
    xa1_display, ya1_display = ax0.transData.transform((xa1, ya1))
    xa2_display, ya2_display = ax0.transData.transform((xa2, ya2))

    slope_D = atan2((yd2_display - yd1_display), (xd2_display - xd1_display))
    gradient_D = slope_D * 180 / pi
    slope_A = atan2((ya2_display - ya1_display), (xa2_display - xa1_display))
    gradient_A = slope_A * 180 / pi
    
    if debug:
        print(f'GradD = {slope_D:4.3f} radians ' +
                    f'({gradient_D:4.2f} degrees), ' +
                    f'GradA = {slope_A:4.3f} radians ' +
                    f'({gradient_A:4.2f} degrees)')
    
    
    # Add acceleration, displacement and origin lines to ('extend') existing lines in the plot
    # ax0.lines.extend(a_lines + d_lines + ax_lines) #, l3, l4])
    for ln in (a_lines + d_lines + ax_lines):
        ax0.add_line(ln)

    def format_txt(txt):
        # txt.replace('e', 'E')
        # txt.replace('-', '\textendash')
        # txt.replace('-', '\textendash')
        return '$'+txt+'$'

    # Add acceleration and displacement text to plot
    [ax0.text(10**T, 10**v, format_txt(txt), 
              transform=ax0.transData, axes=ax0, 
              rotation=gradient_D, rotation_mode='anchor', 
              alpha=0.5, zorder=0.1) 
                for T, v, txt in zip(log_Ta_values, log_va_values, a_str)]
    [ax0.text(10**T, 10**v, format_txt(txt), 
              transform=ax0.transData, axes=ax0, 
              rotation=gradient_A, rotation_mode='anchor',
              alpha=0.5, zorder=0.1) 
                for T, v, txt in zip(log_Td_values, log_vd_values, d_str)]
    
    # Legend
    if debug:
        print(f'Data length: {len(S_v_ranges)}, names:{len(series_names)}')
    if len(series_names) == 0:
        pass
    elif len(series_names) == len(S_v_ranges):
        ax0.legend(series_names)
    else:
        ax0.legend([f'Spectrum {i+1:02d}' for i in range(len(S_v_ranges))])

    # plot any additionally provided points
    if len(plot_pts) > 0:
        [ax0.plot(T, v, marker = 'o') for T,v in plot_pts]
    
    return fig0


def response_spectrum_plots(
        T_ranges: Union[List,
            List_of_lists,
            List_of_list_of_lists
        ], 
        PS_a_ranges: Union[List_of_list_of_lists,
                           List_of_list_of_list_of_lists], 
        series_names: List[str], 
        plot_pts: List_of_lists=[], 
        debug: bool =False, 
        **kwargs: Dict[str, Any]
                            ):
    """Generates a set of pseudo response spectrum plots based on 
    spectral acceleration (units of g)

    Args:
        T_ranges (list[float]: time values - it is assumed that the 
            range is the same for all plots. If an empty list is passed
            to the function, it assumes that all the data is in the
            PS_a_ranges function as T,S pairs.
        PS_a_ranges: pseudo spectral acceleration (assumed to be 
            in units of g). If Period data is in T_ranges, then data 
            should be list of list of lists. Otherwise each data point 
            should be a (T,S) pair.
        series_names (list[str]): one title for each S_a_range
        plot_pts (list[tuple,list]): (T, S_v) points for plotting 
        kwargs: parameters for passing through to? Nowhere at present
    """
    figsize = kwargs.get('figsize', (10,15))
    suptitle = kwargs.get('suptitle', '')
    T_max = kwargs.get('T_max', None)


    if debug:
        print()
        print(f'figsize = {figsize}')
        print(f'suptitle = {suptitle}')
        print('Number of series_names = len(series_names)\n')
    
    # Standardise period and check that numbers align
    if T_ranges == []:
        if debug:
            print('No period data provided')
        if isinstance(
            PS_a_ranges, 
            List_of_list_of_list_of_lists
            ):
            T_ranges, PS_a_ranges = reorganise_plot_data(PS_a_ranges)
        else:
            if debug:
                print('No period data provided as part of PS_a_ranges')
            return None
    elif isinstance(T_ranges[0],(int, float)):
        if debug:
            print('(1) T_ranges is a single collection')
        T_ranges = len(PS_a_ranges) * [T_ranges]
    elif len(T_ranges) == 1 and isinstance(T_ranges[0],(tuple, list)):
        if debug:
            print('(2) T_ranges is a collection containing a single collection')
        T_ranges = len(PS_a_ranges) * T_ranges
    elif len(T_ranges) == len(PS_a_ranges):
        if debug:
            print('(3) T_ranges is a collection of collections')
        pass
    else:
        pass
    
    # DEBUG: print shapes of input collections
    if debug:
        T_ranges_lens = [len(T_range) for T_range in T_ranges]
        print(f'T_ranges shape: {len(T_ranges)} x {T_ranges_lens}')
        PS_a_ranges_lens = [len(PS_a_range) for PS_a_range in PS_a_ranges]
        print(f'PS_a_ranges shape: {len(PS_a_ranges)} x {PS_a_ranges_lens}')

    # CHECK: check that shapes match
    if not (all(len(PS_a_range) == len(T_range) 
                for PS_a_range, T_range in zip(PS_a_ranges, T_ranges))):
        T_ranges_len = [len(T_range) for T_range in T_ranges]
        PS_a_ranges_len = [len(PS_a_range) for PS_a_range in PS_a_ranges]
        err_msg = ('Lengths of curves do not match: ' +
                   f'\n\tlen(T_ranges) = {T_ranges_len} ' + 
                   f'\n\tlen(PS_a_ranges) = {PS_a_ranges_len}\n')
        raise ValueError(err_msg)

    # T_max: truncate collections as necessary
    if T_max is not None:
        if debug:
            print(f'\nTruncation based on T_max = {T_max}:')
        T_list = []
        A_list = []
        for T_range, PS_a_range in zip(T_ranges, PS_a_ranges):
            ta = [(t, a) for t, a in zip(T_range, PS_a_range) if t <= T_max]
            tt, aa = zip(*ta)
            T_list.append(tt)
            A_list.append(aa)
        T_ranges = T_list
        PS_a_ranges = A_list
        # T_ranges = [T for T in T_ranges if T <= T_max]
        # T_num = len(T_ranges)
        # PS_a_ranges = [[2 * pi / T * S_v / 9.81 
        #                 for T, S_v in zip(T_ranges, S_v_range)] 
        #                 for S_v_range in S_v_ranges]
        # PS_a_ranges = [PS_a_range[:T_num] 
        #                for PS_a_range in PS_a_ranges]

        if debug:
            T_ranges_lens = [len(T_range) for T_range in T_ranges]
            print(f'  T_ranges shape: {len(T_ranges)} x {T_ranges_lens}')
            PS_a_ranges_lens = [len(PS_a_range) for PS_a_range in PS_a_ranges]
            print(f'  PS_a_ranges shape: {len(PS_a_ranges)} x {PS_a_ranges_lens}')
            print()
    
    # 
    S_v_ranges = [[9.81 * 0.5 / pi * T * PS_a 
                   for T, PS_a in zip(T_range, PS_a_range)] 
                   for T_range, PS_a_range in zip(T_ranges, PS_a_ranges)]
    # PS_a_ranges = [[2 * pi / T * S_v / 9.81 
    #                 for T, S_v in zip(T_ranges, S_v_range)] 
    #                 for S_v_range in S_v_ranges]
    PS_d_ranges = [[0.5 / pi * T * S_v 
                    for T, S_v in zip(T_range, S_v_range)] 
                    for T_range, S_v_range in zip(T_ranges, S_v_ranges)]

    # DEBUG
    if debug:
        print('\n -- preparing plot --')

    # Preparing plot
    figs, (ax1, ax2, ax3, ax4) = plt.subplots(4,figsize=figsize)
    if len(suptitle) > 0:
        figs.suptitle(suptitle)
    figs.tight_layout()

    # DEBUG
    if debug:
        print('\n -- starting to plot --')

    # Plotting
    [ax1.plot(T_range, PS_a_range) 
     for T_range, PS_a_range in zip(T_ranges, PS_a_ranges)]
    ax1.grid(True)
    ax1.legend(series_names)
    ax1.set_ylabel('Acceleration (g)')
    ax1.set_xlim((0,T_max))
    ax1.set_ylim((0,None))
        
    [ax2.plot(T_range, S_v_range) 
     for T_range, S_v_range in zip(T_ranges, S_v_ranges)]
    ax2.grid(True)
    ax2.legend(series_names)
    ax2.set_ylabel('Velocity (m/s)')
    ax2.set_xlim((0,T_max))
    ax2.set_ylim((0,None))
    
    [ax3.plot(T_range, PS_d_range) 
     for T_range, PS_d_range in zip(T_ranges, PS_d_ranges)]
    ax3.grid(True)
    ax3.legend(series_names)
    ax3.set_xlabel('Period (s)')
    ax3.set_ylabel('Displacement (m)')
    ax3.set_xlim((0,T_max))
    ax3.set_ylim((0,None))
    
    [ax4.plot(PS_d_range, PS_a_range) 
     for PS_d_range, PS_a_range in zip(PS_d_ranges, PS_a_ranges)]
    ax4.grid(True)
    ax4.legend(series_names)
    ax4.set_ylabel('Acceleration (g)')
    ax4.set_xlabel('Displacement (m)')
    ax4.set_xlim((0,None))
    ax4.set_ylim((0,None))

    if len(plot_pts) > 0:
        for pt in plot_pts:
            T, S_v = pt
            PS_a = S_v / 0.5 * pi / T / 9.81  # in m/s
            ax1.plot(T, PS_a, marker = 'o')
            # S_v = 0.5 / pi * T *  PS_a * 9.81  # in m/s
            ax2.plot(T, S_v, marker = 'o')
            PS_d = 0.5 / pi * T *  S_v  # in m
            ax3.plot(T, PS_d, marker = 'o')
            ax4.plot(PS_d, PS_a, marker = 'o')
    
    return figs


def smooth_interpolate(x_vector, y_vector, new_x_values=None, interpolators=None, 
                      plot=False, log_x=True, log_y=True, npoints=200, plot_func=None, **kwargs):
    """
    Perform multiple types of spline interpolation with flexible plotting options.
    
    This function provides a unified interface to compare different interpolation
    methods (Akima, PCHIP, cubic, quadratic) on datasets, with special support
    for response spectra applications. It handles both linear and logarithmic
    interpolation spaces and allows comparison of multiple interpolation techniques.
    
    Parameters
    ----------
    x_vector : array-like
        Original x-coordinates of the data points (e.g., periods for response spectra).
    y_vector : array-like
        Original y-coordinates of the data points (e.g., spectral accelerations).
    new_x_values : array-like, optional
        Specific x-values at which to interpolate. If None, generates a dense grid
        of points between the min and max of x_vector. Default is None.
    interpolators : list of str, optional
        List of interpolator names to use. If None, uses all available interpolators.
        Available options: 
        - 'Akima': Akima1DInterpolator - good balance of smoothness and conservation
        - 'PCHIP': PchipInterpolator - shape-preserving, monotonic between points  
        - 'Cubic (k=3)': Cubic spline - smooth but may overshoot
        - 'Quadratic (k=2)': Quadratic spline - smoother than linear, less oscillatory than cubic
        - 'Linear (k=1)': Linear spline - pievewise linear with sharp corners
        Default is None (uses all interpolators).
    plot : bool, optional
        If True, generates a comparison plot showing original data and all 
        interpolated curves. Default is False.
    log_x : bool, optional
        If True, performs interpolation in log-space for x-values. Recommended for
        response spectra where periods span multiple decades. Default is True.
    log_y : bool, optional
        If True, performs interpolation in log-space for y-values. Recommended for
        response spectra where spectral accelerations span multiple decades. 
        Default is True.
    npoints : int, optional
        Number of points to generate for interpolation when new_x_values is None.
        Default is 200.
    plot_func : function, optional
        Plotting function, e.g. plt.plot or plt.loglog - which will override
        log_x and log_y settings, so that your plot can be in linear space when
        the interpolation is done in log-space (or vice versa).
        Default is None.
    **kwargs
        Additional keyword arguments for plot customization:
        - 'title': Plot title string
        - 'xlim': Tuple for x-axis limits (e.g., (0.1, 10))
        - 'ylim': Tuple for y-axis limits (e.g., (0.01, 1))
        - Any other matplotlib.pyplot plotting parameters
    
    Returns
    -------
    dict
        Dictionary with interpolator names as keys and tuples (x_interp, y_interp) as values.
        Example: {'Akima': (x_akima, y_akima), 'PCHIP': (x_pchip, y_pchip), ...}
        where:
        - x_interp: Array of x-values used for interpolation (same for all methods)
        - y_interp: Array of interpolated y-values for each method
    
    Raises
    ------
    ValueError
        If input arrays have incompatible shapes, contain invalid values for
        the chosen interpolation space, or if unknown interpolator names are specified.
    
    Notes
    -----
    - For response spectra applications, use log_x=True and log_y=True (default)
    - Non-positive values are automatically filtered when using log-space interpolation
    - All interpolators operate in the same transformed space (log or linear) for fair comparison
    - Recommended interpolators for response spectra:
        - 'Akima': General purpose, good balance
        - 'PCHIP': Conservative, shape-preserving
        - 'Quadratic (k=2)': Smooth but less oscillatory than cubic
    
    Examples
    --------
    >>> # Compare all interpolation methods
    >>> results = smooth_interpolate(periods, spectral_data, plot=True)
    >>> akima_x, akima_y = results['Akima']
    >>> pchip_x, pchip_y = results['PCHIP']
    
    >>> # Compare only specific interpolators at target periods
    >>> target_periods = [0.5, 1.0, 1.5, 2.0, 2.5]
    >>> results = smooth_interpolate(
    ...     periods, spectral_data,
    ...     new_x_values=target_periods,
    ...     interpolators=['Akima', 'PCHIP'],
    ...     plot=True,
    ...     title='Response Spectrum Interpolation Comparison'
    ... )
    
    >>> # Linear space interpolation comparison
    >>> results = smooth_interpolate(
    ...     x_data, y_data,
    ...     log_x=False,
    ...     log_y=False,
    ...     interpolators=['Cubic (k=3)', 'Quadratic (k=2)'],
    ...     plot=True,
    ...     xlim=(0, 5),
    ...     ylim=(0, 10)
    ... )
    
    >>> # Get interpolated values without plotting
    >>> results = smooth_interpolate(periods, spectral_data)
    >>> for method, (x_vals, y_vals) in results.items():
    ...     print(f"{method}: {len(y_vals)} interpolated values")
    
    >>> # Semi-log interpolation (log x, linear y)
    >>> results = smooth_interpolate(
    ...     periods, spectral_data,
    ...     log_x=True,
    ...     log_y=False,
    ...     interpolators=['Akima'],
    ...     plot=True
    ... )
    
    See Also
    --------
    akima_interpolate : For single-method Akima interpolation with annotation options
    scipy.interpolate.Akima1DInterpolator : Underlying Akima interpolator
    scipy.interpolate.PchipInterpolator : Underlying PCHIP interpolator
    scipy.interpolate.make_interp_spline : Underlying spline interpolators
    """
    
    # Your data
    x_vector = np.array(x_vector)
    y_vector = np.array(y_vector)
    
    # Filter out non-positive values for log plotting
    if log_x:
        positive_mask = x_vector > 0
        x_vector = x_vector[positive_mask]
        y_vector = y_vector[positive_mask]
    
    if log_y:
        positive_mask = y_vector > 0
        x_vector = x_vector[positive_mask]
        y_vector = y_vector[positive_mask]

    # Generate interpolation x values if not provided
    if new_x_values is None:
        if log_x:
            x_min = np.min(x_vector)
            x_max = np.max(x_vector)
            new_x_values = np.logspace(np.log10(x_min), np.log10(x_max), npoints)
        else:
            new_x_values = np.linspace(x_vector.min(), x_vector.max(), npoints)

    # Convert to appropriate space for interpolation
    if log_x:
        interp_x_vector = np.log10(x_vector)
        interp_new_x = np.log10(new_x_values)
    else:
        interp_x_vector = x_vector
        interp_new_x = new_x_values
        
    if log_y:
        interp_y_vector = np.log10(y_vector)
    else:
        interp_y_vector = y_vector

    # Create interpolators
    interpolators_dict = {
        'Cubic (k=3)': make_interp_spline(interp_x_vector, interp_y_vector, k=3),
        'Quadratic (k=2)': make_interp_spline(interp_x_vector, interp_y_vector, k=2),
        'Linear (k=1)': make_interp_spline(interp_x_vector, interp_y_vector, k=1),
        'Akima': Akima1DInterpolator(interp_x_vector, interp_y_vector),
        'PCHIP': PchipInterpolator(interp_x_vector, interp_y_vector)
    }

    if interpolators is not None:
        interpolators_dict = {name: interpolators_dict[name] for name in interpolators}

    # Perform interpolation
    results = {}
    for name, interp_func in interpolators_dict.items():
        if log_y:
            # Interpolate in log space, then convert back
            log_new_y = interp_func(interp_new_x)
            new_y_values = 10 ** log_new_y
        else:
            # Interpolate directly in linear space
            new_y_values = interp_func(interp_new_x)
        
        results[name] = (new_x_values, new_y_values)

    # Optional plotting
    if plot:
        figsize = kwargs.get('figsize', (12, 8))
        plt.figure(figsize=figsize)
        
        # Plot data points
        if plot_func is not None: # override plot function choice
            plot_func(x_vector, y_vector, 'ko', markersize=8, label='Data points', zorder=5)
        elif log_x and log_y:
            plt.loglog(x_vector, y_vector, 'ko', markersize=8, label='Data points', zorder=5)
        elif log_x:
            plt.semilogx(x_vector, y_vector, 'ko', markersize=8, label='Data points', zorder=5)
        elif log_y:
            plt.semilogy(x_vector, y_vector, 'ko', markersize=8, label='Data points', zorder=5)
        else:
            plt.plot(x_vector, y_vector, 'ko', markersize=8, label='Data points', zorder=5)

        # Plot interpolated curves
        for name, (x_plot, y_plot) in results.items():
            if plot_func is not None: # override plot function choice
                plot_func(x_plot, y_plot, label=name, linewidth=2)
            elif log_x and log_y:
                plt.loglog(x_plot, y_plot, label=name, linewidth=2)
            elif log_x:
                plt.semilogx(x_plot, y_plot, label=name, linewidth=2)
            elif log_y:
                plt.semilogy(x_plot, y_plot, label=name, linewidth=2)
            else:
                plt.plot(x_plot, y_plot, label=name, linewidth=2)

        plt.xlabel(kwargs.get('xlabel', 'Period (s)'))
        plt.ylabel(kwargs.get('ylabel', 'Spectral Acceleration (g)'))
        plt.title(kwargs.get('title', 'Response Spectrum Interpolation Methods Comparison'))
        plt.xlim(kwargs.get('xlim', (0.5, 3)))
        plt.ylim(kwargs.get('ylim', (0.1, 0.3)))
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.3)
        plt.show()

    return results



def akima_interpolate(x_vector, y_vector, new_x_values=None, plot=False, 
                     log_x=True, log_y=True, npoints=500, plot_func=None,
                     annotate_points=None, annotation_style=None, **plot_kwargs):
    """
    Perform Akima interpolation on data points with optional plotting and annotations.
    
    This function is particularly useful for response spectra interpolation where
    smooth curves in log-space are desired. It handles both linear and logarithmic
    interpolation spaces and provides flexible annotation options.
    
    Parameters
    ----------
    x_vector : array-like
        Original x-coordinates of the data points (e.g., periods for response spectra).
    y_vector : array-like
        Original y-coordinates of the data points (e.g., spectral accelerations).
    new_x_values : array-like, optional
        Specific x-values at which to interpolate. If None, generates a dense grid
        of points between the min and max of x_vector. Default is None.
    plot : bool, optional
        If True, generates a plot showing the original data and interpolated curve.
        Default is False.
    log_x : bool, optional
        If True, performs interpolation in log-space for x-values. Recommended for
        response spectra. Default is True.
    log_y : bool, optional
        If True, performs interpolation in log-space for y-values. Recommended for
        response spectra. Default is True.
    npoints : int, optional
        Number of points to generate for interpolation when new_x_values is None.
        Default is 500.
    plot_func : function, optional
        Plotting function, e.g. plt.plot or plt.loglog - which will override
        log_x and log_y settings, so that your plot can be in linear space when
        the interpolation is done in log-space (or vice versa).
        Default is None.
    annotate_points : {None, True, array-like}, optional
        Controls which points to annotate on the plot:
        - None: No annotations (default)
        - True: Automatically selects a subset of interpolated points for annotation
        - array-like: Specific x-values to annotate (e.g., [0.5, 1.0, 1.5])
    annotation_style : dict, optional
        Dictionary controlling the appearance of annotations. Available keys:
        - 'color': Color of the 'x' marker (default: 'red')
        - 'marker_size': Size of the 'x' marker (default: 10)
        - 'text_offset': (x,y) offset for text in points (default: (15, 15))
        - 'bbox_style': Dictionary for text box style (default: yellow rounded box)
        - 'arrow_style': Dictionary for arrow properties (default: red arrow)
        - 'fontsize': Font size for annotations (default: 9)
        - 'precision_x': Decimal precision for x-values (default: 2)
        - 'precision_y': Decimal precision for y-values (default: 3)
    **plot_kwargs
        Additional keyword arguments for plot customization:
        - 'title': Plot title string
        - 'xlim': Tuple for x-axis limits (e.g., (0.1, 10))
        - 'ylim': Tuple for y-axis limits (e.g., (0.01, 1))
    
    Returns
    -------
    tuple
        (x_interpolated, y_interpolated) where:
        - x_interpolated: Array of x-values used for interpolation
        - y_interpolated: Array of interpolated y-values at x_interpolated points
    
    Raises
    ------
    ValueError
        If input arrays have incompatible shapes or contain invalid values for
        the chosen interpolation space.
    
    Notes
    -----
    - For response spectra applications, use log_x=True and log_y=True (default)
    - Non-positive values are automatically filtered when using log-space interpolation
    - The Akima interpolator is chosen for its balance of smoothness and conservation,
      avoiding excessive overshoots common in cubic splines for engineering data
    
    Examples
    --------
    >>> # Basic interpolation without plotting
    >>> periods, accelerations = akima_interpolate(headers, spectral_data)
    
    >>> # Interpolate at specific periods and plot with annotations
    >>> target_periods = [0.5, 1.0, 1.5, 2.0]
    >>> x, y = akima_interpolate(
    ...     headers, spectral_data,
    ...     new_x_values=target_periods,
    ...     plot=True,
    ...     annotate_points=target_periods,
    ...     annotation_style={'color': 'blue', 'precision_x': 2, 'precision_y': 4}
    ... )
    
    >>> # Linear space interpolation for non-spectral data
    >>> x_lin, y_lin = akima_interpolate(
    ...     x_data, y_data,
    ...     log_x=False,
    ...     log_y=False,
    ...     plot=True,
    ...     title="Linear Space Interpolation"
    ... )
    
    >>> # Automatic annotation of interpolated curve
    >>> x_auto, y_auto = akima_interpolate(
    ...     headers, spectral_data,
    ...     plot=True,
    ...     annotate_points=True,
    ...     npoints=100  # Fewer points for cleaner auto-annotation
    ... )
    """
    
    x_vector = np.array(x_vector)
    y_vector = np.array(y_vector)
    
    # Only filter positive values if using log space
    if log_x:
        mask = x_vector > 0
        x_vector = x_vector[mask]
        y_vector = y_vector[mask]
    
    if log_y:
        mask = y_vector > 0
        x_vector = x_vector[mask]
        y_vector = y_vector[mask]
    
    # Convert to appropriate space
    if log_x:
        interp_x = np.log10(x_vector)
    else:
        interp_x = x_vector
        
    if log_y:
        interp_y = np.log10(y_vector)
    else:
        interp_y = y_vector
    
    # Create interpolator
    akima = Akima1DInterpolator(interp_x, interp_y)
    
    # Generate or use provided x values
    if new_x_values is None:
        if log_x:
            new_x_values = np.logspace(np.log10(x_vector.min()), np.log10(x_vector.max()), npoints)
        else:
            new_x_values = np.linspace(x_vector.min(), x_vector.max(), npoints)
    else:
        new_x_values = np.array(new_x_values)
    
    # Convert new x values to interpolation space
    if log_x:
        interp_new_x = np.log10(new_x_values)
    else:
        interp_new_x = new_x_values
    
    # Interpolate
    interp_new_y = akima(interp_new_x)
    
    # Convert back from interpolation space if needed
    if log_y:
        new_y_values = 10 ** interp_new_y
    else:
        new_y_values = interp_new_y
    
    # Optional plotting with annotations
    if plot:
        plt.figure(figsize=(10, 6))
        
        # Plot style based on log/linear settings
        if plot_func is not None:
            plot_func(x_vector, y_vector, 'ko', markersize=6, label='Data points', alpha=0.7)
            plot_func(new_x_values, new_y_values, 'b-', label='Akima interpolation', linewidth=2, alpha=0.8)
        elif log_x and log_y:
            plt.loglog(x_vector, y_vector, 'ko', markersize=6, label='Data points', alpha=0.7)
            plt.loglog(new_x_values, new_y_values, 'b-', label='Akima interpolation', linewidth=2, alpha=0.8)
        elif log_x:
            plt.semilogx(x_vector, y_vector, 'ko', markersize=6, label='Data points', alpha=0.7)
            plt.semilogx(new_x_values, new_y_values, 'b-', label='Akima interpolation', linewidth=2, alpha=0.8)
        elif log_y:
            plt.semilogy(x_vector, y_vector, 'ko', markersize=6, label='Data points', alpha=0.7)
            plt.semilogy(new_x_values, new_y_values, 'b-', label='Akima interpolation', linewidth=2, alpha=0.8)
        else:
            plt.plot(x_vector, y_vector, 'ko', markersize=6, label='Data points', alpha=0.7)
            plt.plot(new_x_values, new_y_values, 'b-', label='Akima interpolation', linewidth=2, alpha=0.8)
        
        # Annotation settings
        default_annotation_style = {
            'color': 'red',
            'marker_size': 10,
            'text_offset': (15, 15),
            'bbox_style': dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8, edgecolor='red'),
            'arrow_style': dict(arrowstyle='->', color='red', lw=1.5),
            'fontsize': 9,
            'precision_x': 2,
            'precision_y': 3
        }
        
        annotation_style = {**default_annotation_style, **(annotation_style or {})}
        
        # Determine which points to annotate
        points_to_annotate = []
        if annotate_points is True:
            # Annotate all interpolated points (use subset for clarity)
            if len(new_x_values) > 20:
                # Select a reasonable number of points for annotation
                indices = np.linspace(0, len(new_x_values)-1, min(8, len(new_x_values)), dtype=int)
                points_to_annotate = list(zip(new_x_values[indices], new_y_values[indices]))
            else:
                points_to_annotate = list(zip(new_x_values, new_y_values))
        elif isinstance(annotate_points, (list, np.ndarray)):
            # Annotate specific x-values
            for x_val in annotate_points:
                if log_x:
                    interp_x_val = np.log10(x_val)
                else:
                    interp_x_val = x_val
                y_val = akima(interp_x_val)
                if log_y:
                    y_val = 10 ** y_val
                points_to_annotate.append((x_val, y_val))
        elif annotate_points is not None:
            # Custom selection logic
            points_to_annotate = annotate_points
        
        # Add annotations
        # if points_to_annotate is not None:
        for i, (x_point, y_point) in enumerate(points_to_annotate):
            plt.plot(x_point, y_point, 'x', color=annotation_style['color'], 
                    markersize=annotation_style['marker_size'], markeredgewidth=2)
            
            label = (f'T={x_point:.{annotation_style["precision_x"]}f} s\n'
                    f'Sa={y_point:.{annotation_style["precision_y"]}f} g')
            
            # Smart offset adjustment for log scales
            text_offset = annotation_style['text_offset']
            if log_x or log_y:
                # Use fixed points offset instead of data coordinates for better behavior
                arrowprops = {**annotation_style['arrow_style'], 'connectionstyle': 'arc3,rad=0.1'}
            else:
                arrowprops = annotation_style['arrow_style']
            
            plt.annotate(label, 
                        xy=(x_point, y_point), 
                        xytext=text_offset, 
                        textcoords='offset points',  # This ensures consistent positioning
                        ha='left', 
                        va='bottom',
                        fontsize=annotation_style['fontsize'],
                        bbox=annotation_style['bbox_style'],
                        arrowprops=arrowprops)

        
        plt.xlabel('Period (s)')
        plt.ylabel('Spectral Acceleration (g)')
        plt.title(plot_kwargs.get('title', 'Response Spectrum with Akima Interpolation'))
        plt.xlim(plot_kwargs.get('xlim', None))
        plt.ylim(plot_kwargs.get('ylim', None))
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.3)
        plt.show()
    
    return new_x_values, new_y_values


# # Fully customized annotations
# custom_style = {
#     'color': 'darkgreen',
#     'marker_size': 12,
#     'text_offset': (20, 20),
#     'bbox_style': dict(boxstyle='round,pad=0.6', facecolor='lightgreen', alpha=0.8, edgecolor='darkgreen'),
#     'arrow_style': dict(arrowstyle='-|>', color='green', lw=2),
#     'fontsize': 10,
#     'precision_x': 3,
#     'precision_y': 5
# }


def EQ_plots(EQ_dict, 
             event_key = 'EVENT_NAME',  # or EVENT_ID
             station_key = 'STATION_CODE', 
             eqdate_key = 'EVENT_DATE_YYYYMMDD',
             eqtime_key = 'EVENT_TIME_HHMMSS',
             eqdir_key = 'STREAM',
             datatype_key = 'DATA_TYPE', 
             data_key = 'Data',  
             events_list = [],
             scaling = 1.0 / 9.81,
             display_units = 'g'
             ):
    """Returns a list of plots based on data in a pair of matching 
    dictionaries of time history and response spectrum

    These dictionaries should be generated by the ASC_import function. The data 
    for 'date', 'station', 'direction' etc assume a certain key format:
        key = "20230206011732_4628" -> Date:2023-02-06 01:17 Station:4628

    It is further assumed that the directions of the motion are 
    
    Args:
        EQ_dict: a dictionary containing earthquake 
            ground acceleration data under the following keys
            - ground motion as list of (t, a) under keys such as 'E', 'N', 'U'.
            - xlabel and ylabel
            and response spectra as list of (T, Sa) under keys such as 'E', 'N', 'U'.
            
    Returns:
        a list of figure objects
    """

    # eq_dates = [f'{eq[:4]}-{eq[4:6]}-{eq[6:8]} {eq[8:10]}:{eq[10:12]}' for eq in TH_dict.keys()]
    # eq_stations = [f'{eq[-4:]}' for eq in TH_dict.keys()]

    keys = (event_key, station_key, eqdate_key, eqtime_key, eqdir_key, datatype_key)
    meta_data = [[v.get(key,None) for key in keys] 
              for k, v in EQ_dict.items()
              if v.get(datatype_key, ['None']).startswith('Acceleration')]

    # events, stations, eqdates, eqtimes, eqdirs, data_types = zip(*meta_data)
    events, stations, *_ = zip(*meta_data)
    if len(events_list) == 0:
        ev_st_list = [(ev, st) for ev in set(events) for st in set(stations)]
    else:
        ev_st_list = [(ev, st) for ev in set(events) for st in set(stations) if ev in events_list]

    meta_dict = {}
    for ev, st in ev_st_list:
        meta_list = meta_dict.get((ev,st), [])
        meta_list.extend([m for m in meta_data if (m[0] == ev and m[1] == st)])
        meta_dict[(ev,st)] = meta_list
    
    # Acquire EQ directions (we need three) and data types (we want TH & RS)
    dates_dict =  {k:list(set(vv[2] for vv in v)) for k,v in meta_dict.items()}
    times_dict =  {k:list(set(vv[3] for vv in v)) for k,v in meta_dict.items()}
    eq_dirs_dict = {k:list(set(vv[4] for vv in v)) for k,v in meta_dict.items()}
    d_types_dict = {k:list(set(vv[5] for vv in v)) for k,v in meta_dict.items()}
    # print('this', eq_dirs, d_types)

    # data_key = event_key
    figs_list = []
    for ev, st in ev_st_list:
        d_types = list(d_types_dict.get((ev, st), []))
        print(ev, st, d_types)
        
        if len(d_types) == 2:
            eq_dirs = eq_dirs_dict.get((ev, st), [])
            eqdate = dates_dict.get((ev, st),'')
            eqtime = times_dict.get((ev, st),'')
            TH_type, RS_type = d_types if 'spectrum' in d_types[1].casefold() else d_types[::-1]
            th_dict = {v.get(eqdir_key,'None'):v.get(data_key,[]) for _, v in EQ_dict.items() if 
                       all([v.get(event_key) == ev,
                           v.get(station_key) == st,
                           v.get(datatype_key) == TH_type,
                           ])}
            rs_dict = {v.get(eqdir_key,'None'):v.get(data_key, []) for _, v in EQ_dict.items() if 
                       all([v.get(event_key) == ev,
                           v.get(station_key) == st,
                           v.get(datatype_key) == RS_type,
                           ])}
            th_data = [th_dict.get(dirn, []) for dirn in eq_dirs]
            rs_data = [rs_dict.get(dirn, []) for dirn in eq_dirs]
            
            # print(th_dict)
            # print(rs_dict)
        
            suptitle = f'Station: {st},    Date: {eqdate[0]+"-"+eqtime[0]}'
            
            # print('Calling plot function')
            # print('lengths', len(th_data), len(rs_data))
            # print('lengths', [len(th) for th in th_data], [len(rs) for rs in rs_data])
            # print('eq_dirs:', eq_dirs)
            # print('suptitle:', suptitle)
            figs_list.append(
                EQ_plot(
                    th_data, rs_data, 
                    eq_name='', 
                    eq_dirs=eq_dirs, 
                    rs_xlim = (0,None), 
                    eq_date = eqdate[0],
                    suptitle = suptitle,
                    scaling=scaling,
                    display_units = display_units
                ))
    
    return figs_list


def EQ_plot(TH_data, RS_data, 
            eq_name='EQ', 
            eq_dirs = ['X', 'Y', 'Z'], 
            eq_date = '', 
            suptitle = '', 
            rs_xlim = (0, None), 
            figsize = (18,18),
            scaling = 1.0 / 9.81,
            display_units = 'g'):
    """Plots earthquake motion data for ground accelerations and 
    response spectra in a 3x2 grid and prints to file (SVG & PNG).
    
    Args:
        TH_data: 
        RS_data: 
        eq_name: is the key to the relevant earthquake in the dictionaries 
        generated by the ASC_import2 function
    'eq_dirs' must be keys in the dictionary referring to direction 
    for each record, 
        e.g. 'E', 'N', 'U'
    """
    
    # print('lengths', len(TH_data), len(RS_data))
    # print('lengths', [len(th) for th in TH_data], [len(rs) for rs in RS_data])
        
    if len(TH_data) == len(RS_data) == len(eq_dirs):
        nrows = len(eq_dirs)
    else:
        raise ValueError('TH_data, RS_data & eq_dirs should have same number of parts')
    
    series_names = []

    fig, axes = plt.subplots(
        nrows=nrows, ncols=2, figsize=figsize
        ) # this is the overall Figure and the Axes object
    fig.suptitle(suptitle)
    
    lts = ['-' for _ in eq_dirs]
    # print('\nlts', lts)
    
    for i, (th_data, rs_data, eqd, ax, lt) in enumerate(zip(TH_data, RS_data, eq_dirs, axes, lts)):
        color = f'C{i}'
        # Plotting the data...
        tE, aE = zip(*th_data) # assumes these come in as SI as default
        aE = [a * scaling for a in aE] # assumes plotting g units as default
        ax[0].hlines(0, 0, tE[-1], colors='k', lw=2)
        ax[0].plot(tE, aE, lt, color=color)
        ax[0].set_title(f'Acceleration - {eqd}')
        ax[0].set_xlabel('Time (s)')
        ax[0].set_ylabel(f'Acceleration ({display_units})')
        ax[0].set_xlim((0, tE[-1]))
        tE, aE = zip(*rs_data)
        aE = [a * scaling for a in aE] # assumes plotting g units as default
        ax[1].plot(tE, aE, lt, color=color)
        ax[1].set_title(f'Response Spectrum - {eqd}')
        series_names.append(f'{eq_date} ({eqd})')
        ax[1].set_xlabel('Period (s)')
        ax[1].set_ylabel(f'Acceleration ({display_units})')
        ax[1].set_xlim(rs_xlim)
        ax[1].set_ylim((0, None))
        
        for axx in ax:
            # axx.grid(True)
            axx.minorticks_on()
            axx.grid(True, which='major', linestyle='-', linewidth=0.5)
            axx.grid(True, which='minor', linestyle='--', linewidth=0.25)
        
        
    fig.tight_layout()
    # fig.savefig(eq_name + '.svg')
    # fig.savefig(eq_name + '.png')

    return fig


def RS_plot(RS_sets:List_of_list_of_lists, 
            data_type:str = 'A', 
            titles:List[str] = [], 
            plot_type:str = 'plot',
            g:float = 9.81, 
            figsize:Tuple[int,int]=(12,16),
            **kwargs) -> Any:
    """
    Generates three plots of A, V, D Response Spectra
    
    Args:
        RS_sets: these are the response spectra data
            which should be as lists of lists of (t, s) pairs
            i.e. if you only want to plot one spectrum you
            should place it as one item in a list, e.g. [spectrum_data]
        data_type: Specifies the type of data in the RS_sets:
            'A', 'V' or 'D' (default is 'A')
        titles: the number of titles should match the 
            number of spectra, otherwise generic titles will
            be substituted
        plot_type: plot type, either 'plot', 'loglog', 'semilogx'
            or 'semilogy' (default is 'plot', i.e. linear)
        g: optional scaling factor to apply to acceleration 
            spectum data (default 9.81), 
        figsize: figure size (default is (12,16))
        

    """
    # xlim = kwargs.get('xlim', (None, None))
    suptitle = kwargs.get('suptitle', None)
    
    if len(titles) != len(RS_sets):
        titles = [f'Spectrum {i:00d}' for i, _ in enumerate(RS_sets)]

    # Functions, scaling factor for acceleration
    # note: tau = 2.0 * pi
    def rspec_ag_to_a(T, a):
        return (T, a * g)
    def rspec_ag_to_v(T, a):
        return (T, a * g * T / tau)
    def rspec_v_to_a(T, v):
        return (T, v * tau / T)
    def rspec_a_to_v(T, a):
        return (T, a * T / tau)
    def rspec_v_to_d(T, v):
        return (T, v * T / tau)
    def rspec_d_to_v(T, d):
        return (T, d * tau / T)

    SA_max = []
    SV_max = []
    SD_max = []

    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=figsize)

    if plot_type == 'loglog':
        ax_plot = [ax.loglog for ax in axs]
    elif plot_type == 'semilogx':
        ax_plot = [ax.semilogx for ax in axs]
    elif plot_type == 'semilogy':
        ax_plot = [ax.semilogy for ax in axs]
    else:
        ax_plot = [ax.plot for ax in axs]
    
    for T_S in RS_sets:
        if data_type == 'A':
            T_SA = [rspec_ag_to_a(t, sa) for t, sa in T_S]
            T_SV = [rspec_a_to_v(t, sa) for t, sa in T_SA]
            T_SD = [rspec_v_to_d(t, sv) for t, sv in T_SV]
        elif data_type == 'V':
            T_SV = list(T_S)
            T_SA = [rspec_v_to_a(t, sv) for t, sv in T_SV]
            T_SD = [rspec_v_to_d(t, sv) for t, sv in T_SV]
        elif data_type == 'D':
            T_SD = T_S
            T_SV = [rspec_d_to_v(t, sd) for t, sd in T_SD]
            T_SA = [rspec_v_to_a(t, sv) for t, sv in T_SV]
        
        SA_max.append(max(T_SA, key=lambda x:x[1]))
        SV_max.append(max(T_SV, key=lambda x:x[1]))
        SD_max.append(max(T_SD, key=lambda x:x[1]))
        
        # legends = ['2% in 50yrs','10% in 50yrs','50% in 50yrs', 'Peak values']
        axs[0].grid(True, 'both')
        ax_plot[0](*list(zip(*T_SA)))
        axs[0].set_ylabel('Spectral Acceleration (m/s^2)')
        axs[1].grid(True, 'both')
        axs[1].set_ylabel('Spectral Velocity (m/s)')
        ax_plot[1](*list(zip(*T_SV)))
        axs[2].grid(True, 'both')
        axs[2].set_xlabel('Period (sec)')
        axs[2].set_ylabel('Spectral Displacement (m)')
        ax_plot[2](*list(zip(*T_SD)))

    titles.append('Peak values')
    ax_plot[0](*zip(*[(t, a) for t, a in SA_max]),'o')
    axs[0].legend(titles)
    ax_plot[1](*zip(*[(t, v) for t, v in SV_max]),'o')
    axs[1].legend(titles)
    ax_plot[2](*zip(*[(t, d) for t, d in SD_max]),'o')
    axs[2].legend(titles)
    # print(*zip(*[(t, a) for t, a in SA_max]))
    # print(SA_max)
    # print(SV_max)
    # print(SD_max)
    if suptitle is not None:
        fig.suptitle(suptitle)
    
    return fig


def EQ_plot_ss(eq, TH_dict, RS_dict, eq_dirs, eq_date = '', 
            suptitle = '', rs_xlim = (0, None)):
    """Plots earthquake motion data for ground accelerations and 
    response spectra in a 3x2 grid and prints to file (SVG & PNG).
    
    'eq' is the key to the relevant earthquake in the dictionaries 
    generated by the ASC_import2 function
    'eq_dirs' must be keys in the dictionary referring to direction 
    for each record, 
        e.g. 'E', 'N', 'U'
    """
    
    S_a_ranges = []
    series_names = []

    th_dict = TH_dict[eq]
    rs_dict = RS_dict[eq]
    # print(th_dict.keys())
    # print(rs_dict.keys())

    fig, [[ax11, ax21], [ax12, ax22], [ax13, ax23]] = plt.subplots(nrows=3, ncols=2, figsize=(18,18)) # this is the overall Figure and the Axes object
    # fig.set_label(eq)
    fig.suptitle(suptitle)
    
    for ax in [ax11, ax21, ax12, ax22, ax13, ax23]:
        ax.grid(True)

    
    axes = [[ax11, ax21], [ax12, ax22], [ax13, ax23]]
    lts = ['-', '-', '-']
    for i, (eqd, ax, lt) in enumerate(zip(eq_dirs, axes, lts)):
        
        color = f'C{i}'
        # Plotting the data...
        tE, aE = zip(*th_dict[eqd])
        ax[0].hlines(0, 0, tE[-1], colors='k', lw=2)
        ax[0].plot(tE, aE, lt, color=color)
        ax[0].set_title(f'Acceleration - {eqd}')
        ax[0].set_xlabel(th_dict.get('xlabel', 'Time (s)'))
        ax[0].set_ylabel(th_dict.get('ylabel', 'Acceleration'))
        ax[0].set_xlim((0, tE[-1]))
        tE, aE = zip(*rs_dict[eqd])
        S_a_ranges.append([a / 980.665 for a in aE]) # assumes cm/s^2
        ax[1].plot(tE, aE, lt, color=color)
        ax[1].set_title(f'Response Spectrum - {eqd}')
        series_names.append(f'{eq_date} ({eqd})')
        ax[1].set_xlabel(rs_dict.get('xlabel', 'Period (s)'))
        ax[1].set_ylabel(rs_dict.get('ylabel', 'Acceleration'))
        ax[1].set_xlim(rs_xlim)
        ax[1].set_ylim((0, None))

        
        
    fig.tight_layout()
    fig.savefig(eq + '.svg')
    fig.savefig(eq + '.png')

    return fig


# def plot_EC8_spectra(a_gR, site_class, source_type=1, 
#                     importance=1.0, q=None, beta=0.2, damp=0.05, T=None, plot_pts=[], 
#                     figsize=(12,9), debug=False, **kwargs):
#     """Spectral horizontal elastic acceleration (EC8, section 3.2.2.2, 
#     Table 3.2 & 3.3, eq 3.2 - 3.5).

#     Args:
#         T: period or list of periods (self-generated range if T=None)
#         a_gR: reference peak ground acceleration on type A ground (in units of g)
#         soil_type: A to E
#         source_type: type is 1 if M_s >= 5.5 else 2 (ref: 3.2.2.2(2P))
#         importance: default is 1.0 (ref:2.1(2-3,P))
#         damp: should be as a decimal not as a percentage (e.g. 0.05 for 5 % damping)
    
#     Returns:
#         T_s_list: tuple of tuple of periods, T & tuple of spectral accelerations (g) 

#     """
#     x_limits = kwargs.pop('x_limits', (None, None))
#     y_limits = kwargs.pop('y_limits', (None, None))

#     Te_range, Se_range = S_e(T, a_gR, site_class, source_type, 
#                                 importance=importance, damp=damp)
#     legend = ['Elastic Spectrum']
#     fig0, ax0 = plt.subplots(figsize=figsize)
#     ax0.set_title(f'EC8 Soil Type {site_class}')
#     ax0.set_xlabel('Building Period (s)')
#     ax0.set_ylabel('Base Shear Coefficient (g)')
#     ax0.plot(Te_range, Se_range)
#     ax0.grid(True)
#     ax0.set_xlim(*x_limits)
#     ax0.set_ylim(*y_limits)
#     if q is not None:
#         Td_range, Sd_range = S_d(T, a_gR, site_class, source_type, 
#                                     importance=importance, beta=beta, q=q)
#         ax0.plot(Td_range, Sd_range) 
#         legend.append(f'Design Spectrum (q={q})')
#     if len(plot_pts) > 0:
#         s_e = S_e(plot_pts, a_gR, site_class, source_type, 
#                     importance=importance, damp=damp)
#         ax0.plot(plot_pts, s_e)
#         if q is not None:
#             s_d = S_d(plot_pts, a_gR, site_class, source_type, 
#                         importance=importance, beta=beta, q=q)
#             ax0.plot(plot_pts, s_d)
#             if debug:
#               print(f's_d is {s_d}')
    
    
#     ax0.legend(legend)

#     return fig0


def _interp_sa(T_query, T_spec, Sa_spec):
    """Interpolate spectral acceleration values onto *T_query*."""
    return np.interp(np.asarray(T_query), np.asarray(T_spec), np.asarray(Sa_spec))


def _resolve_spectrum(entry: dict, T_plot: np.ndarray, g: float = 9.81) -> np.ndarray:
    """Return Sa (m/s\u00b2) for a single spectrum entry.

    Entries with pre-computed data (``'T'`` / ``'Sa'`` keys) are interpolated
    onto *T_plot*.  Otherwise *entry* must have a ``'code'`` key that
    dispatches to the appropriate module-level builder.
    """
    # --- Pre-computed data ------------------------------------------------
    if "T" in entry and "Sa" in entry:
        T_data = np.asarray(entry["T"], dtype=float)
        Sa_data = np.asarray(entry["Sa"], dtype=float)
        if entry.get("in_g", False):
            Sa_data = Sa_data * g
        return np.interp(T_plot, T_data, Sa_data)

    code = entry.get("code", "")
    params = entry.get("params", entry)

    # --- GB 50011 ---------------------------------------------------------
    if code == "GB50011":
        from pySeismic.codes.GB50011 import GB_get_alpha_max, GB_design_spectrum

        alpha_max = GB_get_alpha_max(
            intensity=params.get("intensity", 7),
            acceleration=params.get("acceleration", 0.10),
            level=params.get("level", "rare"),
        )
        return GB_design_spectrum(
            T=T_plot,
            alpha_max=alpha_max,
            site_class=params.get("site_class", "II"),
            group=params.get("group", 1),
            zeta=params.get("damping", 0.05),
        )

    # --- EC8 --------------------------------------------------------------
    if code == "EC8":
        from pySeismic.codes.EC8 import EC8_design_spectrum

        sa_func = EC8_design_spectrum(
            a_gR=params.get("a_gR", 0.12),
            site_class=params.get("site_class", "C"),
            source_type=params.get("source_type", 1),
            importance=params.get("importance", 1.0),
            damp=params.get("damp", 0.05),
            q=params.get("q", None),
            beta=params.get("beta", 0.2),
        )
        return sa_func(T_plot) * g  # convert from g to m/s\u00b2

    # --- ASCE 7 -----------------------------------------------------------
    if code == "ASCE7":
        from pySeismic.codes.ASCE7_16 import design_response_spectrum

        return (
            design_response_spectrum(
                S_DS=params.get("S_DS", 0.5),
                S_D1=params.get("S_D1", 0.3),
                T=T_plot,
                T_L=params.get("T_L", 12.0),
            )
            * g
        )  # convert from g to m/s\u00b2

    # --- ASCE 41 (same shape as ASCE 7) -----------------------------------
    if code == "ASCE41":
        from pySeismic.codes.ASCE7_16 import design_response_spectrum

        return (
            design_response_spectrum(
                S_DS=params.get("S_DS", 0.5),
                S_D1=params.get("S_D1", 0.3),
                T=T_plot,
                T_L=params.get("T_L", 12.0),
            )
            * g
        )

    raise ValueError(
        f"Unknown spectrum code {code!r}. Supported: GB50011, EC8, ASCE7, ASCE41"
    )


def plot_seismic_spectrum(
    spectra: Sequence[Dict[str, Any]],
    g: float = 9.81,
    figsize: Tuple[float, float] = (9, 5),
    T_max: float = 6.0,
    xlim: Tuple[float, float] | None = None,
    n_pts: int = 300,
    xlabel: str = "Period T (s)",
    ylabel: str | None = None,
    title: str = "Response Spectrum",
    unit: Literal["m/s\u00b2", "g"] = "m/s\u00b2",
    log_x: bool = False,
    log_y: bool = False,
    grid: bool = True,
    legend_fontsize: float | None = None,
    modal: dict | None = None,
    return_data: bool = False,
    **kwargs: Any,
) -> Figure | Tuple[Figure, dict]:
    """Plot one or more response spectra from any seismic code.

    Parameters
    ----------
    spectra : list of dict
        Each dict defines one spectrum to plot.  Two formats are supported:

        **Pre-computed data**
        ::
            {"T": [...], "Sa": [...], "label": "...", "color": "...", ...}

        **Code-specific parameters**
        ::
            {"code": "GB50011",
             "params": {"intensity": 7, "acceleration": 0.10, ...},
             "label": "...", "color": "...", ...}

        Entries may also include ``"in_g": True`` if *Sa* values are in g.

        Supported ``code`` values: ``GB50011``, ``EC8``, ``ASCE7``, ``ASCE41``.

        Style keys recognised per entry: ``color``, ``linestyle``,
        ``linewidth``, ``alpha``, ``marker``, ``zorder``.

    g : float
        Acceleration of gravity (m/s\\u00b2).
    figsize : tuple of float
        Figure size (width, height) in inches.
    T_max : float
        Maximum period on the x-axis (used when *xlim* is ``None``).
    xlim : tuple of float or None
        Explicit (xmin, xmax) period limits.  Overrides *T_max* when set.
    n_pts : int
        Number of points for the period axis.
    xlabel, ylabel : str
        Axis labels.  If *ylabel* is ``None``, auto-set from *unit*.
    title : str
        Figure suptitle.
    unit : ``'m/s\\u00b2'`` or ``'g'``
        Display unit for spectral acceleration.
    log_x, log_y : bool
        If ``True``, use logarithmic axes.
    grid : bool
        Whether to show a grid.
    legend_fontsize : float or None
        Font size for the legend.
    modal : dict or None
        If provided, vertical dashed lines mark dominant periods.
    return_data : bool
        If ``True``, return ``(fig, data_dict)`` where *data_dict* maps each
        label to ``(T_plot, Sa)`` arrays (Sa in m/s\\u00b2).

    Returns
    -------
    matplotlib.figure.Figure or (Figure, dict)
    """
    # --- Resolve y-label from unit -------------------------------------------
    if ylabel is None:
        ylabel = f"Sa ({unit})"

    scale = 1.0 if unit == "m/s\\u00b2" else (1.0 / g)

    T_plot = np.linspace(0.01, T_max, n_pts)

    fig, ax = plt.subplots(figsize=figsize)

    # Bind the correct plot method
    if log_x and log_y:
        ax_plot = ax.loglog
    elif log_x:
        ax_plot = ax.semilogx
    elif log_y:
        ax_plot = ax.semilogy
    else:
        ax_plot = ax.plot

    data_out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    for entry in spectra:
        # Resolve Sa values (m/s\\u00b2), then scale to display unit
        Sa = _resolve_spectrum(entry, T_plot, g) * scale

        label = entry.get("label", "")
        data_out[label] = (T_plot.copy(), Sa.copy() / scale)

        style = {
            k: entry[k]
            for k in ("color", "linestyle", "linewidth", "alpha", "marker", "zorder")
            if k in entry
        }
        ax_plot(T_plot, Sa, label=label, **style)

    # Vertical lines for fundamental periods (modal analysis)
    if modal is not None:
        periods = modal.get("periods", [])
        mp = modal.get("modal_props", {})
        for dir_label, ratio_key in [
            ("X", "partiMassRatiosMX"),
            ("Y", "partiMassRatiosMY"),
        ]:
            ratios = mp.get(ratio_key, [])
            if not ratios:
                continue
            best = max(range(len(ratios)), key=lambda i: abs(ratios[i]))
            if best < len(periods):
                T_dom = periods[best]
                ax.axvline(T_dom, color="grey", linestyle="--",
                           linewidth=0.8, alpha=0.7)
                ax.text(T_dom, ax.get_ylim()[1] * 0.95,
                        f"T\\u2081({dir_label})={T_dom:.3f}s",
                        fontsize=8, rotation=90, va="top", ha="right",
                        alpha=0.7)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if not log_x:
        ax.set_xlim(*(xlim if xlim is not None else (0, T_max)))
    if grid:
        ax.grid(True, alpha=0.3)
    if spectra:
        ax.legend(fontsize=legend_fontsize or 8)
    fig.tight_layout()

    if return_data:
        return fig, data_out
    return fig


def main():
    
    # Usage examples:
    # log_headers = [log(hdr) for hdr in headers]
    SPEC_HEADERS = [0.01,0.2,0.5,1.0,2.0,3.0,4.0,5.0,6.0]
    vector = [0.075, 0.148, 0.22, 0.25, 0.223, 0.126, 0.067, 0.047, 0.038]

    # 1. Just get interpolated values without plotting
    results = smooth_interpolate(SPEC_HEADERS, vector, interpolators=['Akima'])
    akima_x, akima_y = results['Akima']
    # print(f"Akima x-values: {akima_x}")
    # print(f"Akima y-values: {akima_y}")

    # 2. Get values at specific periods
    target_periods = [0.5, 1.0, 1.5, 2.0, 2.5]
    results = smooth_interpolate(SPEC_HEADERS, vector, new_x_values=target_periods, interpolators=['Akima'])
    akima_values = results['Akima'][1]  # y-values at target_periods
    print(f"\nAkima x-values: {np.array(target_periods)}")
    print(f"\nAkima y-values: {akima_values}")

    # 3. Plot with interpolation
    # ['Cubic (k=3)', 'Quadratic (k=2)', 'Linear (k=1)', 'Akima', 'PCHIP']
    results = smooth_interpolate(SPEC_HEADERS, vector, plot=True, interpolators=['Linear (k=1)', 'Akima', 'PCHIP'])

    # 4. Log plot with linear interpolation
    # ['Cubic (k=3)', 'Quadratic (k=2)', 'Linear (k=1)', 'Akima', 'PCHIP']
    results = smooth_interpolate(SPEC_HEADERS, vector, plot=True, interpolators=['Linear (k=1)', 'Akima', 'PCHIP'], log_x=False, log_y=False, plot_func=plt.loglog)


    # Usage examples:

    # 1. Just get interpolated values
    x_smooth, y_smooth = akima_interpolate(SPEC_HEADERS, vector)

    # 2. Get values at specific periods
    target_periods = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    x_target, y_target = akima_interpolate(SPEC_HEADERS, vector, new_x_values=target_periods)
    print(f"Akima values at target periods: {y_target}")

    # 3. Plot results
    x_plot, y_plot = akima_interpolate(SPEC_HEADERS, vector, plot=True, 
                                    xlim=(0.5, 3), ylim=(0.1, 0.3),
                                    title='Site-Specific Response Spectrum')
    
    # 4. Plot results on linear 
    x_plot, y_plot = akima_interpolate(SPEC_HEADERS, vector, plot=True, 
                                    xlim=(0, None), ylim=(0, None), 
                                    plot_func=plt.plot,
                                    title='Site-Specific Response Spectrum (log-interpolation, linear scale)')
    
    # 5. Plot results on log-log scale
    # the_dict = spec_dict.get(('ESF', 250), {})
    # vector = the_dict.get(zone, None)
    x_plot, y_plot = akima_interpolate(SPEC_HEADERS, vector, plot=True, 
                                    xlim=(0, None), ylim=(0, None), 
                                    plot_func=plt.plot,
                                    title='Site-Specific Response Spectrum (log-interpolation, linear scale)')

    # 6. Test with log-log plotting & annotation
    x_smooth, y_smooth = akima_interpolate(
        SPEC_HEADERS, vector,
        annotate_points=[0.5, 1.0, 2.0],  # Annotate these periods
        plot=True,
        log_x=True, 
        log_y=True,  # This will use plt.loglog()
        annotation_style={
            'text_offset': (15, 15),  # Points offset works consistently
            'color': 'red'
        }
    )

    #7. Test tripartite plot
    fig = tripartite_plot(SPEC_HEADERS, [vector], plot_pts=[(0.5, 0.2), (1.0, 0.25)]) #, plot=True)
    plt.show()

    print()
    print('*'*40)
    plot_storey_VM2(
        storey_elevs=(2,4,6,8), 
        storey_shears=[(50,40,30,20)], 
        storey_moments=[(100,92,70,42)], 
        base_level=0.5, debug = True,
        # **kwargs
        )
    plt.show()

if __name__ == "__main__":
    main()
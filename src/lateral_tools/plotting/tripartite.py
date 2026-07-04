"""
tripartite.py — Tripartite (log-log) response spectrum plots.

Part of the lateral-tools package for seismic and lateral load engineering.

Functions
---------
tripartite_plot         : Main tripartite plot with acceleration/displacement grid.
tripartite_plot_old     : Legacy tripartite implementation.
reorganise_plot_data    : Transpose nested (T, Sv) point data for plotting.
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
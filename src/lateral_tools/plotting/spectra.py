"""
spectra.py — Spectra plotting functions.

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
        from lateral_tools.codes.GB50011 import GB_get_alpha_max, GB_design_spectrum

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
        from lateral_tools.codes.EC8 import EC8_design_spectrum

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
        from lateral_tools.codes.ASCE7_16 import design_response_spectrum

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
        from lateral_tools.codes.ASCE7_16 import design_response_spectrum

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


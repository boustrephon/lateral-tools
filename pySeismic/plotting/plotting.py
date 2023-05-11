"""
These are plotting functions for seismic engineering. The data needs to 
be generated elsewhere, for example from eqsigs.

"""

from typing import Union, List, Tuple, Dict, Any
from math import log, pi, tau, sin, cos, tan, acos, asin, atan, atan2, log, log10
from itertools import accumulate
import matplotlib.pyplot as plt
import matplotlib.lines as lines


def reorganise_plot_data(list_of_lists):
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
    transposed_data = zip(*[new_list for new_list in [zip(*each_list) for each_list in list_of_lists]])
    return transposed_data


def tripartite_plot(
    T_ranges, 
    S_v_ranges, 
    series_names = [],
    plot_pts=[], 
    debug=False, 
    figsize=(10,8), 
    suptitle='', 
    minor_lines = [1, 2, 3, 4, 5, 6, 7, 8, 9],
    **kwargs
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
        T_ranges: collection of time values (each as a collection of floats)
            if there is only one item, it is assumed that the range is the same for all plots
        S_v_ranges: collection of spectral velocities (each as a collection of floats)
        series_names (list[str]): one name for each S_v_range
        plot_pts (list[tuple,list]): (T, S_v) points for plotting 
        kwargs: parameters for passing through to `plot`
    """
    # One of the challenges is the need to switch between plain and log10 values
    # 
    
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
    if debug: print(f'Figure size is {figsize}')
    if debug: print(f'Super-title is {suptitle}')
    
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
    a_str = [f'{a:5.1g}\,{a_unit}' for a in a_list]
    a_values = [a * a_fac for a in a_list]
    log_a_values = [log10(a) for a in a_values]
    # [print(a,b,s) for a,b,s in zip(a_list, a_values, a_str)]
    

    #== Displacement grid lines =========
    # generate the list of displacement line values provided they are within the limits
    d_list = [(m * 10**float(n)) for n in range(int(log10(d_min)), int(log10(d_max))+1) 
            for m in minor_lines 
            if (m * 10**float(n) * d_fac > d_min) and (m * 10**float(n) / d_fac < d_max)]
    d_str = [f'{d:5.1g}\,{d_unit}' for d in d_list]
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
    
    if debug: print(f'GradD = {slope_D:4.3f} radians ' +
                    f'({gradient_D:4.2f} degrees), ' +
                    f'GradA = {slope_A:4.3f} radians ' +
                    f'({gradient_A:4.2f} degrees)')
    
    
    # Add acceleration, displacement and origin lines to ('extend') existing lines in the plot
    ax0.lines.extend(a_lines + d_lines + ax_lines) #, l3, l4])

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


def response_spectrum_plots(T_ranges, PS_a_ranges, series_names, 
                            plot_pts=[], debug=False, **kwargs):
    """Generates a set of pseudo response spectrum plots based on 
    spectral acceleration (units of g)

    Args:
        T_ranges (list[float]: time values - it is assumed that the 
            range is the same for all plots 
        PS_a_ranges: pseudo spectral acceleration (assumed to be 
            in units of g)
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
        print(f'Number of series_names = len(series_names)\n')
    
    # Standardise period and check that numbers align
    if isinstance(T_ranges[0],(int, float)):
        if debug: print('(1) T_ranges is a single collection')
        T_ranges = len(PS_a_ranges) * [T_ranges]
    elif len(T_ranges) == 1 and isinstance(T_ranges[0],(tuple, list)):
        if debug: print('(2) T_ranges is a collection containing a single collection')
        T_ranges = len(PS_a_ranges) * T_ranges
    elif len(T_ranges) == len(PS_a_ranges):
        if debug: print('(3) T_ranges is a collection of collections')
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
        err_msg = (f'Lengths of curves do not match: ' +
                   f'\n\tlen(T_ranges) = {T_ranges_len} ' + 
                   f'\n\tlen(PS_a_ranges) = {PS_a_ranges_len}\n')
        raise ValueError(err_msg)

    # T_max: truncate collections as necessary
    if T_max is not None:
        if debug: print(f'\nTruncation based on T_max = {T_max}:')
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
    if debug: print(f'\n -- preparing plot --')

    # Preparing plot
    figs, (ax1, ax2, ax3, ax4) = plt.subplots(4,figsize=figsize)
    if len(suptitle) > 0:
        figs.suptitle(suptitle)
    figs.tight_layout()

    # DEBUG
    if debug: print(f'\n -- starting to plot --')

    # Plotting
    [ax1.plot(T_range, PS_a_range) 
     for T_range, PS_a_range in zip(T_ranges, PS_a_ranges)]
    ax1.grid(True)
    ax1.legend(series_names)
    ax1.set_ylabel(f'Acceleration (g)')
    ax1.set_xlim((0,T_max))
    ax1.set_ylim((0,None))
        
    [ax2.plot(T_range, S_v_range) 
     for T_range, S_v_range in zip(T_ranges, S_v_ranges)]
    ax2.grid(True)
    ax2.legend(series_names)
    ax2.set_ylabel(f'Velocity (m/s)')
    ax2.set_xlim((0,T_max))
    ax2.set_ylim((0,None))
    
    [ax3.plot(T_range, PS_d_range) 
     for T_range, PS_d_range in zip(T_ranges, PS_d_ranges)]
    ax3.grid(True)
    ax3.legend(series_names)
    ax3.set_xlabel(f'Period (s)')
    ax3.set_ylabel(f'Displacement (m)')
    ax3.set_xlim((0,T_max))
    ax3.set_ylim((0,None))
    
    [ax4.plot(PS_d_range, PS_a_range) 
     for PS_d_range, PS_a_range in zip(PS_d_ranges, PS_a_ranges)]
    ax4.grid(True)
    ax4.legend(series_names)
    ax4.set_ylabel(f'Acceleration (g)')
    ax4.set_xlabel(f'Displacement (m)')
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


def EQ_plots(EQ_dict, 
             event_key = 'EVENT_NAME',  # or EVENT_ID
             station_key = 'STATION_CODE', 
             eqdate_key = 'EVENT_DATE_YYYYMMDD',
             eqtime_key = 'EVENT_TIME_HHMMSS',
             eqdir_key = 'STREAM',
             datatype_key = 'DATA_TYPE', 
             data_key = 'Data',  
             events_list = [],
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
                    suptitle = suptitle
                ))
    
    return figs_list


def EQ_plot(TH_data, RS_data, 
            eq_name='EQ', 
            eq_dirs = ['X', 'Y', 'Z'], 
            eq_date = '', 
            suptitle = '', 
            rs_xlim = (0, None), 
            figsize = (18,18)):
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
    
    S_a_ranges = []
    series_names = []

    # th_dict = TH_dict[eq]
    # rs_dict = RS_dict[eq]
    # print(th_dict.keys())
    # print(rs_dict.keys())

    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=figsize) # this is the overall Figure and the Axes object
    # fig.set_label(eq)
    fig.suptitle(suptitle)
    
    # for ax in [ax11, ax21, ax12, ax22, ax13, ax23]:
    #     ax.grid(True)

    
    # axes = [[ax11, ax21], [ax12, ax22], [ax13, ax23]]
    
    lts = ['-' for _ in eq_dirs]
    # print('\nlts', lts)
    
    for i, (th_data, rs_data, eqd, ax, lt) in enumerate(zip(TH_data, RS_data, eq_dirs, axes, lts)):
        color = f'C{i}'
        # print(f'colour is {color}')
        # print('Calling plot function')
        # print('lengths', len(th_data), len(rs_data))
        # print('lengths', [len(th) for th in th_data][:5], [len(rs) for rs in rs_data][:5])
        # print('eq_dir:', eqd)
        # print('axs', len(ax))
            
        # Plotting the data...
        tE, aE = zip(*th_data)
        # print('lengths', len(tE), len(aE))
        # print('data', tE[:5], aE[:5])
        
        ax[0].hlines(0, 0, tE[-1], colors='k', lw=2)
        ax[0].plot(tE, aE, lt, color=color)
        ax[0].set_title(f'Acceleration - {eqd}')
        ax[0].set_xlabel('Time (s)')
        ax[0].set_ylabel('Acceleration')
        ax[0].set_xlim((0, tE[-1]))
        tE, aE = zip(*rs_data)
        S_a_ranges.append([a / 980.665 for a in aE]) # assumes cm/s^2
        ax[1].plot(tE, aE, lt, color=color)
        ax[1].set_title(f'Response Spectrum - {eqd}')
        series_names.append(f'{eq_date} ({eqd})')
        ax[1].set_xlabel('Period (s)')
        ax[1].set_ylabel('Acceleration')
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


def RS_plot(RS_sets:List[List[float]], 
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
    xlim = kwargs.get('xlim', (None, None))
    suptitle = kwargs.get('suptitle', None)
    
    if len(titles) != len(RS_sets):
        titles = [f'Spectrum {i:00d}' for i, _ in enumerate(RS_sets)]

    # Functions, scaling factor for acceleration
    # note: tau = 2.0 * pi
    rspec_ag_to_a = lambda T, a: (T, a * g)
    rspec_ag_to_v = lambda T, a: (T, a * g * T / tau)
    rspec_v_to_a = lambda T, v: (T, v * tau / T)
    rspec_a_to_v = lambda T, a: (T, a * T / tau)
    rspec_v_to_d = lambda T, v: (T, v * T / tau)
    rspec_d_to_v = lambda T, d: (T, d * tau / T)

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
#             if debug: print(f's_d is {s_d}')
    
    
#     ax0.legend(legend)

#     return fig0

# ========================================================
# =================== Storey Plots =======================
# ========================================================

from pySeismic.utilities import storey_force_distribution, \
    storey_shear_moment, shear_coords
from matplotlib.ticker import EngFormatter

def plot_storey_VM(
        base_shears: Union[float, List, Tuple], 
        storey_elevs: Union[List[float], Tuple[float]], 
        storey_wts: Union[List[float], Tuple[float]], 
        k:float = 1, 
        base_level: float = 0, 
        **kwargs: Dict):
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
    plot_title = kwargs.pop('plot_title', f'Building Shears & Moments')
    figsize = kwargs.pop('figsize', (6,8))
    force_units = kwargs.pop('force_units', 'kN')
    length_units = kwargs.pop('length_units', 'm')
    debug = kwargs.pop('debug', False)

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
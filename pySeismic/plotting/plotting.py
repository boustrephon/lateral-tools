"""
These are plotting functions for seismic engineering. The data needs to 
be generated elsewhere, for example from eqsigs.

"""

from math import log, pi, sin, cos, tan, acos, asin, atan, atan2, log, log10
from itertools import accumulate
import matplotlib.pyplot as plt
import matplotlib.lines as lines



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


def EQ_plots(TH_dict, RS_dict):
    """Returns a list of plots based on data in a pair of matching 
    dictionaries of time history and response spectrum

    These dictionaries should be generated by the ASC_import2 function. The data 
    for 'date', 'station', 'direction' etc assume a certain key format:
        key = "20230206011732_4628" -> Date:2023-02-06 01:17 Station:4628

    It is further assumed that the directions of the motion are 
    
    Args:
        TH_dict: a dictionary containing earthquake ground acceleration data under the following keys
            - ground motion as list of (t, a) under keys such as 'E', 'N', 'U'.
            - xlabel and ylabel
        RS_dict: 
            - response spectra as list of (T, Sa) under keys such as 'E', 'N', 'U'.
            - 
    Returns:
        a list of figure objects
    """

    # eq_dates = [f'{eq[:4]}-{eq[4:6]}-{eq[6:8]} {eq[8:10]}:{eq[10:12]}' for eq in TH_dict.keys()]
    # eq_stations = [f'{eq[-4:]}' for eq in TH_dict.keys()]

    figs_list = []
    for eq in TH_dict.keys():
        eq_date = f'{eq[:4]}-{eq[4:6]}-{eq[6:8]} {eq[8:10]}:{eq[10:12]}'
        eq_station = f'{eq[-4:]}'
        suptitle = f'Station: {eq_station},    Date: {eq_date}'
        
        figs_list.append(EQ_plot(
            eq, TH_dict, RS_dict, ('E', 'N', 'U'), 
            rs_xlim = (0,None), 
            eq_date = eq_date,
            suptitle = suptitle
        ))
    
    return figs_list



def EQ_plot(eq, TH_dict, RS_dict, eq_dirs, eq_date = '', 
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

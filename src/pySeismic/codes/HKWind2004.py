"""2004 Wind Code

Note that a building is considered to be one with resonant dynamic response if:
  a) the height exceeds five times the least dimension
  b) the height of the building is greater than 100m
  Unless the fundamental natural frequency of the building can be shown to be
  greater than 1.0Hz
"""

from math import log  #, pi, sin, cos, tan, acos, asin, atan, atan2, log10
from itertools import accumulate
import doctest
from sys import path as sys_path
from os.path import abspath, join as path_join, dirname
import warnings
# warnings.filterwarnings("ignore", category=DeprecationWarning)
# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

from pySeismic.utilities import list_interp
# from pySeismic.utilities import print_table, interp, \
#     storey_force_coeffs, storey_force_distribution, storey_forces, \
#     storey_shear_moment, shear_coords


# ===============================
# ====== General Functions ======
# ===============================

def get_storey_shear(storey_forces, elevations_rising=True):
    """elevations rising means forces are ordered with rising storey levels
    
    """
    if elevations_rising:
        return list(accumulate(storey_forces[::-1]))[::-1]
    else:
        return list(accumulate(storey_forces))


def get_storey_moment(storey_shears, storey_heights, elevations_rising=True):
    
    """elevations rising means forces are ordered with rising storey levels
    
    """
    if elevations_rising:
        return list(accumulate(V * h for V, h in zip(storey_shears[::-1], storey_heights[::-1])))[-2::-1] + [0]
    else:
        return [0] + (accumulate(V * h for V, h in zip(storey_shears, storey_heights)))[:-1]


def gen_ff_heights(h, h_list):
    """
    
    """
    # [5.5] + [3.5 for i in range(int((100-5.5)//3.5))]
    return h_list[:-1] + [h_list[-1] for i in range(int((h - sum(h_list[:-1])) // h_list[-1]))]


def gen_elevs(list_of_storey_heights):
    """
    
    """
    return list(accumulate(list_of_storey_heights))


def gen_trib_heights(list_of_storey_heights):
    """
    
    """
    if list_of_storey_heights[0] == 0:
        return [0.5 * (a + b) for b, a in zip(list_of_storey_heights,
                                        list_of_storey_heights[1:] + [0.0])]
    else:
        return [0.5 * (a + b) for b, a in zip([0.0] + list_of_storey_heights,
                                        list_of_storey_heights + [0.0])]


def gen_trib_areas(list_of_storey_heights, list_of_storey_widths):
    """
    
    """
    if isinstance(list_of_storey_widths, list):
        h_len = len(list_of_storey_heights)
        w_len = len(list_of_storey_widths)
        if h_len != w_len:
            err_msg = f'number of storey heights ({h_len}) does not match number of storey widths({w_len})'
            raise ValueError(err_msg)
        # TODO - fix 'list_of_zip'
        list_of_storey_areas = [h * w for h, w in zip(list_of_storey_heights, list_of_storey_widths)]
        if list_of_storey_heights[0] == 0:
            return [0.5 * (a + b) for b, a in zip(list_of_storey_areas,
                                                list_of_storey_areas[1:] + [0.0])]
        else:
            return [0.5 * (a + b) for b, a in zip([0.0] + list_of_storey_areas,
                                            list_of_storey_areas + [0.0])]
    #elif isinstance(list_of_storey_widths,(float,int)):
    else:
        w = list_of_storey_widths
        if list_of_storey_heights[0] == 0:
            return [0.5 * w * (a + b) for b, a in zip(list_of_storey_heights,
                                                list_of_storey_heights[1:] + [0.0])]
        else:
            return [0.5 * w * (a + b) for b, a in zip([0.0] + list_of_storey_heights,
                                          list_of_storey_heights + [0.0])]


# =================================
# ====== Wind Code Functions ======
# =================================


def get_1989_q(z, zone='General'):
    """Returns wind pressures (kPa) from HK 1983 Code
    
    The two zones are 'General' & 'Built-Up'

    z: height above artificial base (m)

    Returns:
        Design Wind Pressure (kPa)
    """
    if zone == 'General':
        heights = (0,10,30,50,100,150,200,250,600)
        pressures = (0.0,1.2,2.2,2.5,3.0,3.5,3.8,4.1,4.3)
    else:
        heights = (0,30,50,100,150,200,250,300, 600)
        pressures = (0.0,1.2,1.9,2.4,3.0,3.4,3.7,4.0,4.3)
    if z >= heights[-1]: 
        return pressures[-1]
    else: 
        return [p for p, h in zip(pressures, heights) if h > z][0]


def get_2004_q(z):
    """Wind pressures (kPa) from HK 2004 Code
    """
    heights = (0,5,10,20,30,50,75,100,150,200,250,300,400,500,600)
    pressures = (1.82,1.82,2.01,2.23,2.37,2.57,2.73,
                 2.86,3.05,3.2,3.31,3.41,3.58,3.72,3.72)
    if z < 0:
        return pressures[0]
    elif z >= heights[-1]: 
        return pressures[-1]
    else: 
        return list_interp(z, heights, pressures, extend=True)



def get_v_z_mean(z, z_g=500, alpha=0.11, v_g_mean=56.6):
    """hourly-mean wind velocity at height z (eq 1.5)"""
    z = min(max(z, 5), 500)
    return 1.05 * v_g_mean * (z / z_g)**alpha


def get_v_z(z, z_g=500, alpha=0.11, v_g_mean=56.6, g_v=3.7, I_g=0.087):
    """design 3 second gust wind speed v_z at height z (eq 1.2, 1.3, 1.5 combined)
    
    See section 1.21 of explanatory materials

    z - elevation in metres
    z_g - reference height in metres (500m)
    alpha - the power coefficient defining the wind gradient (0.11)
    v_g_mean - hourly-mean wind velocity at reference height in m/s (56.6m/s)
    g_v - the peak factor, normally taken between 3.4 and 3.7
    I_g - turbulence intensity at gradient height = 0.087
    """
    z_ratio = min(max(z, 5), 500) / z_g

    if z > 0:
        return 1.05 * v_g_mean * z_ratio**alpha * \
        (1 + g_v * I_g * z_ratio**(-alpha))
    else:
        return 1.05 * v_g_mean * g_v * I_g


def get_q_z(z, rho=1.2):
    """design wind pressure
    z - elevation in metres
    rho - air density in kg / m^3 (1.2)
    """
    return 0.0005 * rho * get_v_z(z)**2


def get_F_z(C_f, q_z, A_z):
    """Along-wind force for an enclosed building at level z (static) in kN
    C_f - force coefficient (product of C_h and C_s, Appendix D for enclosed buildings)
    q_z - design pressure (kPa)
    A_z - area perpendicular to wind velocity (m^2)
    """
    return C_f * q_z * A_z


def get_F_p(C_p, q_z, A_m):
    """Along-wind pressure for individual elements at level z (static) in kN
    C_p - pressure coefficient (for individual elements, Appendix E)
    q_z - design pressure (kPa)
    A_m - element projected area perpendicular to wind velocity (m^2)
    """
    return C_p * q_z * A_m


def get_G(h, b, n_a=None, zeta = 0.02, g_v=3.7):
    """Dynamic Magnification Factor, G (Appendix F)
    h - building height in metres
    b - building breadth (dimension perpendicular to the wind) in metres
    n_a - building frequency in Hz
    """
    # h = min(max(h, 100), 200) # should h be limited?
    V_h_mean = get_v_z_mean(h)
    n_a = (46 / h) if (n_a is None) else n_a
    I_h = 0.1055 * (h / 90)**(-0.11)
    g_f = (2 * log(3600 * n_a))**0.5
    L_h = 1000 * (h / 10)**0.25
    B = 1 / (1 + (36 * h**2 + 64 * b**2)**0.5 / L_h)
    N = n_a * L_h / V_h_mean
    E = 0.47 * N / (2 + N**2)**(5/6)
    S = 1 / (1 + 3.5 * n_a * h / V_h_mean) / (1 + 4.0 * n_a * b / V_h_mean)
    return 1 + 2 * I_h * (g_v**2 * B + g_f**2 * S * E / zeta)**0.5


def get_F_z_dyn(C_f, q_z_mean, A_z, G):
    """Along-wind force for an enclosed building at level z (dynamic effects)
    G -  Dynamic Magnification Factor
    Note that area force reduction should not be applied
    """
    return G * C_f * q_z_mean * A_z


def get_C_h(h_b):
    """Shape factors (Table D2)
    h_b - building height divided by 
        building breadth (perpendicular to wind) (dimensionless)
    """
    h_bs = [1.0, 2.0, 4.0, 6.0, 10.0, 20.0]
    C_hs = [0.95, 1.0, 1.05, 1.1, 1.2, 1.4]
    return list_interp(h_b, h_bs, C_hs, extend=True)


def get_C_s(b_d):
    """Shape factors (Table D2)
    b_d - building breadth (perpendicular to wind) divided by 
        building depth in direction of wind (dimensionless)
    """
    b_ds = [1.0, 2.0, 3.0]
    C_ss = [1.0, 1.1, 1.3]
    return list_interp(b_d, b_ds, C_ss, extend=True)


def get_R_A(A):
    """Reduction factors (Table D3)
    A - frontal projected area (in metres squared)
    """
    Areas = [500.0, 800.0, 1000.0, 3000.0, 5000.0, 8000.0, 10000.0, 15000.0]
    R_As = [1.00, 0.97, 0.96, 0.92, 0.89, 0.86, 0.84, 0.80]
    return list_interp(max(500,min(A,15000)), Areas, R_As, extend=True)


def full_wind_calcs(storey_heights, 
                    storey_widths_x, storey_widths_y, 
                    print_data=False, raise_error=False,
                    old_format=False, debug=False):
    """Returns full set of calculations for HK Wind Code 2004"""
    # This could be speeded up significantly by using numpy
    #print(f'Bottom storey height is {elevs[0]} m, top storey elevation is {elevs[-1]} m, sum of storeys is {sum(storey_heights)} m')
    dd = {}
    dd['storey_heights'] = storey_heights
    dd['storey_elevs'] = list(accumulate(dd['storey_heights']))
    dd['force_elevs'] = [z - h / 2 for z, h in zip(dd['storey_elevs'], dd['storey_heights'])]

    # if storey widths (x) are not lists, make them lists
    if not isinstance(storey_widths_x, list):
        bx_ave = float(storey_widths_x)
        dd['storey_widths_x'] = [bx_ave for _ in storey_heights]
    else:
        if len(storey_widths_x) != len(storey_heights):
            raise ValueError('Number of storey heights does not match number of storey widths (x)')
        bx_ave = sum(storey_widths_x) / len(storey_widths_x)
        dd['storey_widths_x'] = storey_widths_x
        
    # if storey widths (y) are not lists, make them lists
    if not isinstance(storey_widths_y, list):
        by_ave = float(storey_widths_y)
        dd['storey_widths_y'] = [by_ave for _ in storey_heights]
    else:
        if len(storey_widths_y) != len(storey_heights):
            raise ValueError('Number of storey heights does not match number of storey widths (y)')
        by_ave = sum(storey_widths_y) / len(storey_widths_y)
        dd['storey_widths_y'] = storey_widths_y
        
    # Calculate storey areas (tributary areas), velocities and pressures

    dd['trib_areas_x'] = gen_trib_areas(dd['storey_heights'], dd['storey_widths_x'])
    dd['trib_areas_y'] = gen_trib_areas(dd['storey_heights'], dd['storey_widths_y'])
    dd['storey_mean_vels'] = [get_v_z(z) for z in dd['storey_elevs']] 
    dd['storey_pressures'] = [get_q_z(v) for v in dd['storey_mean_vels']]

    # building height (based on top storey)
    h = max(dd['storey_elevs'])

    # Get height shape factors based on average width
    C_hx = get_C_h(h / by_ave)
    C_hy = get_C_h(h / bx_ave)

    # Get horizontal plane shape factors by storey
    dd['C_sx'] = [get_C_s(by / bx) for bx, by in zip(dd['storey_widths_x'], dd['storey_widths_y'])]
    dd['C_sy'] = [get_C_s(bx / by) for bx, by in zip(dd['storey_widths_x'], dd['storey_widths_y'])]
    if debug:
        print(f'bx_ave: {bx_ave:4.4g}; by_ave: {by_ave:4.4g}')
        print(f'C_hx: {C_hx:4.4g}; C_hy: {C_hy:4.4g}')
        print(f'C_sx: {dd["C_sx"]}; C_sy: {dd["C_sy"]}')

    # Calculate horizontal forces
    dd['storey_forces_x'] = [get_F_z(C_hx + C, p, a_y) for p, a_y, C in zip(dd['storey_pressures'], dd['trib_areas_y'], dd['C_sx'])]
    dd['storey_forces_y'] = [get_F_z(C_hy + C, p, a_x) for p, a_x, C in zip(dd['storey_pressures'], dd['trib_areas_x'], dd['C_sy'])]
    if debug:
        # 'storey_pressures,   trib_areas,     C_s,     storey_forces'
        _ = [print(f'{col:16s}| len: {len(dd[col])}') for col in ('storey_pressures', 'trib_areas_y', 'C_sx', 'storey_forces_x' )]
        _ = [print(f'{col:16s}| len: {len(dd[col])}') for col in ('storey_pressures', 'trib_areas_x', 'C_sy', 'storey_forces_y' )]

    dd['storey_shears_x'] = get_storey_shear(dd['storey_forces_x'])
    dd['storey_shears_y'] = get_storey_shear(dd['storey_forces_y'])
    dd['storey_moments_yy'] = get_storey_moment(dd['storey_shears_x'], dd['storey_heights'])
    dd['storey_moments_xx'] = get_storey_moment(dd['storey_shears_y'], dd['storey_heights'])

    ## END OF CALCULATIONS #######

    # Check that maximum shear matches sum of forces and 
    # maximum moment is the sum of the products of the storey forces and storey elevations
    sum_of_storey_forces_x = sum(dd['storey_forces_x'])
    sum_of_storey_moments_yy = sum(F * h for F, h in zip(dd['storey_forces_x'], dd['storey_elevs']))
    sum_of_storey_forces_y = sum(dd['storey_forces_y'])
    sum_of_storey_moments_xx = sum(F * h for F, h in zip(dd['storey_forces_y'], dd['storey_elevs']))

    err_list = []
    max_storey_shear_x = dd['storey_shears_x'][0]
    max_storey_shear_y = dd['storey_shears_y'][0]
    max_storey_moment_xx = dd['storey_moments_xx'][0] + max_storey_shear_y * dd['storey_elevs'][0]
    max_storey_moment_yy = dd['storey_moments_yy'][0] + max_storey_shear_x * dd['storey_elevs'][0]

    if round(sum_of_storey_forces_x) != round(max_storey_shear_x):
        err_msg = f'Error in sum of x-dir forces ({sum_of_storey_forces_x} vs {max_storey_shear_x})'
        err_list.append(err_msg)
    if round(sum_of_storey_moments_yy) != round(max_storey_moment_yy):
        err_msg = f'Error in sum of yy-dir moments ({sum_of_storey_moments_yy} vs {max_storey_moment_yy})'
        err_list.append(err_msg)
    if round(sum_of_storey_forces_y) != round(max_storey_shear_y):
        err_msg = f'Error in sum of y-dir forces ({sum_of_storey_forces_y} vs {max_storey_shear_y})'
        err_list.append(err_msg)
    if round(sum_of_storey_moments_xx) != round(max_storey_moment_xx):
        err_msg = f'Error in sum of xx-dir moments ({sum_of_storey_moments_xx} vs {max_storey_moment_xx})'
        err_list.append(err_msg)
    if len(err_list) > 0:
        err_msg = '\n'.join(err_list)
        if raise_error:
            raise ValueError(err_msg)
        else:
            warnings.warn(err_msg)
            # print(err_msg)
    
    if debug:

        trib_x = dd['trib_areas_x'][0]
        trib_y = dd['trib_areas_y'][0]
        columns = ['storey_elevs', 'storey_heights', 
                    'storey_widths_x', 'trib_areas_x', 'C_sx', 'storey_forces_x', 'storey_shears_x', 'storey_moments_yy', 
                    'storey_widths_y', 'trib_areas_y', 'C_sy', 'storey_forces_y', 'storey_shears_y', 'storey_moments_xx']
        new_row = [0, 0, 0, trib_x, 0, 0, max_storey_shear_x, max_storey_moment_yy, 0, trib_y, 0, 0, max_storey_shear_y, max_storey_moment_xx]
        # columns = ['storey_elevs', 'storey_heights', 
        #             'storey_widths_x', 'trib_areas_x', 'storey_forces_x', 'storey_shears_x', 'storey_moments_yy', 
        #             'storey_widths_y', 'trib_areas_y', 'storey_forces_y', 'storey_shears_y', 'storey_moments_xx']
        # new_row = [0, 0, 0, trib_x, 0, max_storey_shear_x, max_storey_moment_yy, 0, trib_y, 0, max_storey_shear_y, max_storey_moment_xx]
        # columns = ['storey_elevs', 'storey_heights', 
        #            'storey_forces_x', 'storey_shears_x', 'storey_moments_yy', 
        #             'storey_forces_y', 'storey_shears_y', 'storey_moments_xx']
        # new_row = [0, 0, 0, max_storey_shear_x, max_storey_moment_yy, 0, max_storey_shear_y, max_storey_moment_xx]
        try:
            import pandas as pd
            from IPython.display import display
            df = pd.DataFrame({k:(v[1:] if k.startswith('trib_areas_') else v) for k, v in dd.items() if k in columns})
            # New row data (as a dictionary or list)
            new_row_data = {k:v for k, v in zip(columns, new_row)}

            # Convert the new row to a DataFrame
            # The index can be explicitly set to [0] for clarity before concatenation
            new_row_df = pd.DataFrame([new_row_data])

            # Concatenate the new row DataFrame with the original DataFrame
            # The new_row_df comes first to be at the top
            df_updated = pd.concat([new_row_df, df], ignore_index=True)
            display(df_updated)

        except (ImportError, ValueError) as e:
            print(f'keys: {dd.keys()}')
            print(f'storey_heights: {dd["storey_heights"]}')
            print(f'storey_elevs: {dd["storey_elevs"]}')
            print(f'storey_widths_x: {dd["storey_widths_x"]}')
            print(f'trib_areas_x: {dd["trib_areas_x"]}')
            print(f'storey_forces_x: {dd["storey_forces_x"]}')
            print(f'storey_shears_x: {dd['storey_shears_x']}')
            print(f'storey_moments_yy: {dd["storey_moments_yy"]}')
            print(f'Base Myy: {max_storey_moment_yy}')
            print(f'storey_widths_y: {dd["storey_widths_y"]}')
            print(f'trib_areas_y: {dd["trib_areas_y"]}')
            print(f'storey_forces_y: {dd["storey_forces_y"]}')
            print(f'storey_shears_y: {dd["storey_shears_y"]}')
            print(f'storey_moments_xx: {dd["storey_moments_xx"]}')
            print(f'Base Mxx: {max_storey_moment_xx}')
            print()
            print(f'Error: {e}')
            if debug:
                _ = [print(f'{col:18s}| len: {len(dd[col])}') for col in columns]


    ## PRINT RESULTS if print_data is True #######
    if print_data:

        print(f'C_hx : {C_hx:4.3f}   C_hy : {C_hy:4.3f}\n')

        [print(f'{a:4.2f} m, {ax:4.2f} m, {ay:4.2f} m, {v:4.2f} m/s, {p:4.2f} kPa') \
        for a, ax, ay, v, p in zip(dd['storey_heights'], dd['trib_areas_x'], 
                                    dd['trib_areas_y'], dd['storey_mean_vels'], dd['storey_pressures'])]
        print()

        [print(f'{a:4.2f} m, {ax:4.2f} m, {ay:4.2f} m, {v:4.2f} m/s, {p:4.2f} kPa, {x:7,.0f} kN, {y:7,.0f} kN') \
        for a, ax, ay, v, p, x, y in zip(dd['storey_heights'], dd['trib_areas_x'], dd['trib_areas_y'], 
                                        dd['storey_mean_vels'], dd['storey_pressures'], dd['storey_forces_x'], dd['storey_forces_y'])]
        print()

        print('Vx,base = {:7,.0f} kN,    Myy,base = {:7,.0f} kNm'.format(sum_of_storey_forces_x, sum_of_storey_moments_yy))
        print('Vy,base = {:7,.0f} kN,    Mxx,base = {:7,.0f} kNm'.format(sum_of_storey_forces_y, sum_of_storey_moments_xx))
        print()

        [print(f'{h:4.2f} m, {z:4.2f} m, {Fx:7,.0f} kN, {Vx:7,.0f} kN, {Myy:7,.0f} kNm, {Vy:7,.0f} kN, {Mxx:7,.0f} kNm') \
        for h, z, Fx, Vx, Myy, Vy, Mxx in zip(dd['storey_heights'], dd['storey_elevs'], dd['storey_forces_x'], 
                                            dd['storey_shears_x'], dd['storey_moments_yy'], dd['storey_shears_y'], dd['storey_moments_xx'])]
        print()

    max_storey_moment_xx = dd['storey_moments_xx'][0] + dd['storey_shears_y'][0] * dd['storey_elevs'][0]
    max_storey_moment_yy = dd['storey_moments_yy'][0] + dd['storey_shears_x'][0] * dd['storey_elevs'][0]

    if not old_format:
        dd['storey_moments_yy'] = [max_storey_moment_yy] + dd['storey_moments_yy'][:-1]
        dd['storey_moments_xx'] = [max_storey_moment_xx] + dd['storey_moments_xx'][:-1]
    
    ### PRINT RESULTS if print_data is True #######
    if print_data:
    
        print(f'C_hx : {C_hx:4.3f}   C_hy : {C_hy:4.3f}\n')

        [print(f'{a:4.2f} m, {ax:4.2f} m, {ay:4.2f} m, {v:4.2f} m/s, {p:4.2f} kPa') \
            for a, ax, ay, v, p in zip(dd['storey_heights'], dd['trib_areas_x'], 
                                        dd['trib_areas_y'], dd['storey_mean_vels'], dd['storey_pressures'])]
        print()

        [print(f'{a:4.2f} m, {ax:4.2f} m, {ay:4.2f} m, {v:4.2f} m/s, {p:4.2f} kPa, {x:7,.0f} kN, {y:7,.0f} kN') \
            for a, ax, ay, v, p, x, y in zip(dd['storey_heights'], dd['trib_areas_x'], dd['trib_areas_y'], 
                                            dd['storey_mean_vels'], dd['storey_pressures'], dd['storey_forces_x'], dd['storey_forces_y'])]
        print()

        print('Vx,base = {:7,.0f} kN,    Myy,base = {:7,.0f} kNm'.format(sum_of_storey_forces_x, sum_of_storey_moments_yy))
        print('Vy,base = {:7,.0f} kN,    Mxx,base = {:7,.0f} kNm'.format(sum_of_storey_forces_y, sum_of_storey_moments_xx))
        print()
        # print('?'*20)
        # for key in ['storey_heights', 'storey_elevs', 'storey_forces_x', 'storey_shears_x', 'storey_moments_yy', 'storey_shears_y', 'storey_moments_xx']:
        #     val = dd[key][0] # Check the first item of each
        #     print(f"{len(dd[key])} | {key}: {type(val)} - {val}")
        # print('?'*20)
        [print(f'{h:4.2f} m, {z:4.2f} m, {Fx:7,.0f} kN, {Vx:7,.0f} kN, {Myy:7,.0f} kNm, {Vy:7,.0f} kN, {Mxx:7,.0f} kNm') \
            for h, z, Fx, Vx, Myy, Vy, Mxx in zip(dd['storey_heights'], dd['storey_elevs'], dd['storey_forces_x'], 
                                                dd['storey_shears_x'], dd['storey_moments_yy'], dd['storey_shears_y'], dd['storey_moments_xx'])]
        print()

    return dd


heights = [0, 5, 10, 20, 30, 
           50, 75, 100, 150, 
           200, 250, 300, 400, 
           500, 600]


def print_v_z_mean_table(heights = heights):
    print('Table F3: Design hourly-mean wind velocity')
    v_z_mean = [get_v_z_mean(z) for z in heights]
    [print(f'{z:4.0f} m | {v:4.2f} m/s') for z, v in zip(heights, v_z_mean)]
    pass


def print_table_G(zeta,n_table=''):
    """

    """
    if n_table != '':
        table_txt = f'Table F{n_table} :'
    print(f'{table_txt} Dynamic Magnification Factor G for ζ = {100*zeta:4.1f}%')
    [print([round(get_G(h, b, zeta = zeta),3) 
    for b in [20, 30, 40]])
    for h in [200, 180, 160, 140, 120, 100]]


def print_table_F1():
    """

    """
    zeta = 0.015
    print_table_G(zeta,1)


def print_table_F2():
    """

    """
    zeta = 0.02
    print_table_G(zeta,2)



heights = [0, 5, 10, 20, 30, 
            50, 75, 100, 150, 
            200, 250, 300, 400, 
            500, 600]


def print_table_F3(heights=heights):
    """

    """
    print('Table F3: Design hourly-mean wind velocity')
    v_z_mean = [get_v_z_mean(z) for z in heights]
    [print(f'{z:4.0f} m | {v:4.2f} m/s') for z, v in zip(heights, v_z_mean)]


def print_table_3sec_wind(heights=heights):
    """

    """
    print('Design 3 second gust wind speed at height z')
    v_z = [get_v_z(z) for z in heights]
    [print(f'{z:4.0f} m | {v:4.2f} m/s') for z, v in zip(heights, v_z)]


def print_table_1(heights=heights):
    """

    """
    v_z = [get_v_z(z) for z in heights]
    print('Table 1: Design Wind Pressure (kPa)')
    [print(f'{z:4.0f} m | {v:4.2f} m/s ' + \
        f'| {get_q_z(v):4.2f} kPa') 
        for z, v in zip(heights, v_z)]


def main():
    pass


if __name__ == "__main__":
    doctest.testmod()
    print('== Doctest Complete ==')
    main()

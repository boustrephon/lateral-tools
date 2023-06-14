"""2019 Wind Code

Note that a building is considered to be one with resonant dynamic response if:
  a) the height exceeds five times the least dimension
  b) the height of the building is greater than 100m
  Unless the fundamental natural frequency of the building can be shown to be
  greater than 1.0Hz
"""

from math import pi, sin, cos, tan, acos, asin, atan, atan2, \
    log, log10, exp
from itertools import accumulate

from pySeismic.utilities import printt, interp, list_interp, \
    storey_force_coeffs, storey_force_distribution, storey_forces, \
    storey_shear_moment, shear_coords


# ===============================
# ====== General Functions ======
# ===============================




# =================================
# ====== Wind Code Functions ======
# =================================

def get_v_z_mean(z, z_g=500, alpha=0.11, v_g_mean=56.6):
    """hourly-mean wind velocity at height z (eq 1.5??)"""
    z = min(max(z, 5), 500)
    return 1.05 * v_g_mean * (z / z_g)**alpha


def get_v_z(z, z_g=500, alpha=0.11, v_g_mean=56.6, g_v=3.7, I_g=0.087):
    """design 3 second gust wind speed v_z at height z (eq 1.2, 1.3, 1.5 combined)
    z - elevation in metres
    z_g - reference height in metres
    alpha - the power coefficient defining the wind gradient
    v_g_mean - hourly-mean wind velocity at reference height in m/s
    g_v - the peak factor, normally taken between 3.4 and 3.7
    I_g - 
    """
    z = min(max(z, 5), 500)
    if z > 0:
        return 1.05 * v_g_mean * (z / z_g)**alpha * (1 + g_v * I_g * (z / z_g)**(-alpha))
    else:
        return 1.05 * v_g_mean * g_v * I_g


def get_q_z(v_z, rho=1.2):
    """design wind pressure
    z - elevation in metres
    rho - air density in kg / m^3
    """
    return 0.0005 * rho * v_z**2


def get_G(h, b, n_a=None, zeta = 0.02, g_v=3.7):
    """Dynamic Magnification Factor, G (Appendix F)
    h - building height in metres
    b - building breadth (dimension perpendicular to the wind) in metres
    n_a - building frequency in Hz
    """
    # h = min(max(h, 100), 200) # should h be limited?
    V_h_mean = get_v_z_mean(h)
    n_a = n_a if n_a else (46 / h)
    I_h = 0.1055 * (h / 90)**(-0.11)
    g_f = (2 * log(3600 * n_a))**0.5
    L_h = 1000 * (h / 10)**0.25
    B = 1 / (1 + (36 * h**2 + 64 * b**2)**0.5 / L_h)
    N = n_a * L_h / V_h_mean
    E = 0.47 * N / (2 + N**2)**(5/6)
    S = 1 / (1 + 3.5 * n_a * h / V_h_mean) / (1 + 4.0 * n_a * b / V_h_mean)
    return 1 + 2 * I_h * (g_v**2 * B + g_f**2 * S * E / zeta)**0.5


def get_table_3_1(Ze):
    """Ze is the effective height for open exposure,
    conservatively taken as the height above ground level"""
    Zes = (2.5,5,10,20,30,50,75,100,150,200,250,300,400,500)
    Q_ozs = (1.59,1.77,1.98,2.21,2.36,2.56,2.73,2.86,3.05,3.20,3.31,3.41,3.57,3.70)
    if Ze <= 500: 
        raise ValueError('Seek specialist advice for Ze > 500m')
    return list_interp(Ze, Zes, Q_ozs) 


def get_Q_oz(Ze):
    """Ze is the effective height for open exposure,
    conservatively taken as the height above ground level"""
    if Ze > 500: 
        raise ValueError('Seek specialist advice for Ze > 500m')
    return 3.7 * (Ze / 500)**0.15 # Eq 3-2


def get_I_oz(Ze):
    """Ze is the effective height for open exposure,
    conservatively taken as the height above ground level"""
    if Ze > 500: 
        raise ValueError('Seek specialist advice for Ze > 500m')
    return 0.087 * (Ze / 500)**-0.11 # Eq 3-3


def get_C_f(He_D, B_D):
    """Returns 
    He_D 
    B_D
    """
    He_D = min(12, He_D)
    return 1.1 + 0.55 * He_D / (
        exp(log(0.6 * B_D * (1 - 0.011 * He_D)**(1.7 - 0.0013 * (He_D**2))))
    )


table_B3_1_dict = {
    3: {'A': 2.3, 'B': 1.4, 'C':1.2, 'D': None},
    5: {'A': 2.9, 'B': 1.8, 'C':1.4, 'D': 1.2},
    10: {'A': 3.4, 'B': 2.1, 'C':1.7, 'D': 1.2},
    }


def freestanding_wall_pressure(phi, l_h, lc_h):
    """
    TO DO - review - what is l_h??
    """
    if phi < 0.8:
        raise ValueError('')
    if lc_h >= 1.0:
        Cp_dict = {'A': 2.1, 'B': 1.8, 'C':1.4, 'D': 1.2}
        return Cp_dict
    return None


def get_table_4_2(phi):
    """"""
    phis = (0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 0.9, 1.0)
    C_ps = (2.0, 1.8, 1.7, 1.6, 1.5, 1.5, 1.6, 2.0)
    return list_interp(phi, phis, C_ps)


def get_Ss(L05p, zone=None):
    """Appendix C1, eq C1-1"""
    Ss = exp(0.17 - 0.07 * L05p**0.32)
    if zone == 'Edge':
        return max(Ss, 1.3 - log(L05p) / 9.0)
    elif zone == 'Corner':
        return max(Ss, 1.5 - log(L05p) / 5.4)
    else:
        return Ss

def cantilever_deflection(w, L, E, I):
    return 0.125 * w * L**4 / E / I

def cantilever_forces(p, B, L):
    V = p * B * L
    return V, 0.5 * V * L


# =============================
# ====== Plot Functions =======
# =============================


import matplotlib.pyplot as plt

def plot_Ss(Lmax = 3, L05p = 6.1, n=30, debug=False):
    """
    Appendix C1, eq C1-1
    """
    
    L_range = [10**(i/n*Lmax) for i in range(n+1)]
    if debug:
        print(f'Q_oz(6.1) is {get_Q_oz(6.1):.3f}')
        print(f'S_s(6.1) is {get_Ss(6.1):.3f}')
        print(get_Ss(7.5), exp(0.17 - 0.07 * 7.5**0.32), 
            get_Ss(12), exp(0.17 - 0.07 * 12**0.32))
    S_range = [get_Ss(L,'Corner') for L in L_range]
    plt.semilogx(L_range, S_range)
    S_range = [get_Ss(L, 'Edge') for L in L_range]
    plt.semilogx(L_range, S_range)
    S_range = [get_Ss(L) for L in L_range]
    plt.semilogx(L_range, S_range)
    plt.legend(['Corner', 'Edge', 'Other'])
    plt.grid(True)
    plt.scatter(L05p, get_Ss(L05p), color = '#88c999')


def print_table_F1():
    """

    """
    print('Table F1 : Dynamic Magnification Factor G for ζ = 2.0%')
    [[round(get_G(h, b, zeta = 0.015),3) 
    for b in [20, 30, 40]] 
    for h in [200, 180, 160, 140, 120, 100]]

heights = [0, 5, 10, 20, 30, 
            50, 75, 100, 150, 
            200, 250, 300, 400, 
            500, 600]


def print_table_F2():
    """

    """
    print('Table F2 : Dynamic Magnification Factor G for ζ = 1.5%')
    [[round(get_G(h, b, zeta = 0.02),3) 
      for b in [20, 30, 40]] 
      for h in [200, 180, 160, 140, 120, 100]]


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


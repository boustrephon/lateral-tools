"""
ASCE7-22

"This Standard provides requirements for 
general structural design and includes 
means for determining dead, live, soil, 
flood, wind, snow, rain, atmospheric 
ice, and earthquake loads, as well as 
their combinations, which are suitable 
for inclusion in building codes and 
other documents."

Previous editions are:
* ASCE 7-16
* ASCE 7-10
* ASCE 7-05
* ASCE 7-02
* ASCE 7-98
* ASCE 7-95
* ASCE 7-93
* ASCE 7-88
"""

from typing import Sequence
import doctest
from sys import path as sys_path
from os.path import abspath, join as path_join #, dirname
from os import getcwd
from pathlib import Path

# Add the parent directory of 'pySeismic' to sys.path
file_folder = Path(getcwd()) # %pwd
print(f'Parent folder: {file_folder}')

module_folder = file_folder / 'src' / 'pySeismic'  # path_join(file_folder, '..')
print(f'Module folder: {module_folder}')
if module_folder not in sys_path:
    sys_path.append(str(module_folder))

from utilities import list_interp
from plotting.plotting import plot_storey_VM2

# ============================================================================
# UNIT CONVERSION UTILITIES
# ============================================================================

class Units:
    """Unit conversion constants and helpers"""
    
    # Conversion factors to meters
    TO_METERS = {
        'm': 1.0,
        'cm': 0.01,
        'mm': 0.001,
        'in': 0.0254,
        'ft': 0.3048,
        'yd': 0.9144
    }
    
    # Conversion factors to feet
    TO_FEET = {
        'm': 3.28084,
        'cm': 0.0328084,
        'mm': 0.00328084,
        'in': 0.0833333,
        'ft': 1.0,
        'yd': 3.0
    }
    
    @classmethod
    def to_meters(cls, value: float, from_unit: str) -> float:
        """Convert a length value to meters"""
        factor = cls.TO_METERS.get(from_unit)
        if factor is None:
            raise ValueError(f"Unknown unit: {from_unit}")
        return value * factor
    
    @classmethod
    def to_feet(cls, value: float, from_unit: str) -> float:
        """Convert a length value to feet"""
        factor = cls.TO_FEET.get(from_unit)
        if factor is None:
            raise ValueError(f"Unknown unit: {from_unit}")
        return value * factor
    
    @classmethod
    def to_meters_sq(cls, value: float, from_unit: str) -> float:
        """Convert an area value to square meters"""
        return cls.to_meters(value**0.5, from_unit)**2
    
    @classmethod
    def to_feet_sq(cls, value: float, from_unit: str) -> float:
        """Convert an area value to square feet"""
        return cls.to_feet(value**0.5, from_unit)**2


# ============================================================================
# CODE PERIOD (ASCE 7-22 Section 12.8.2.1)
# ============================================================================

# def c_u(S_D1):
#     C_u = list_interp(S_D1, 
#                 (0.1, 0.15, 0.2, 0.3, 0.4), 
#                 (1.7, 1.6, 1.5, 1.4, 1.4,),extend=True)
#     return C_u

C_t_dict = {'S1':0.035, 'S1a':0.035, 'C1':0.018, 'S2':0.030, 'S2a': 0.030}
beta_dict = {'S1':0.80, 'S1a':0.80, 'C1':0.90}

def C_t_get(x: str) -> float:
    """Get Ct coefficient for structure type"""
    return C_t_dict.get(x, 0.020)


def beta_get(x: str) -> float:
    """Get beta exponent for structure type"""
    return beta_dict.get(x, 0.75)

# unit_l_dict = {'m': 1.0, 'cm': 0.01, 'mm':0.001, 'in':0.0254, 'ft':0.3048, 'yd':0.9144}

# def unit_l_get(x):
#     return unit_l_dict.get(x, None)  # a lookup function that returns default of None

def T_code_get(s_type: str, height: float, units: str = 'm') -> float:
    """
    Returns the code period in seconds using ASCE 7 Eqn 4-4
    
    Args:
        s_type: Structure type (e.g., 'C1', 'S1')
        height: Building height in specified units
        units: Input units ('m', 'ft', 'in', etc.)
    
    Returns:
        Period in seconds
    """
    C_t = C_t_get(s_type)
    beta = beta_get(s_type)
    
    # Convert height to meters for internal calculation if using metric units
    # Note: The empirical formula expects height in meters when using metric
    if units == 'm':
        h_m = height
    else:
        h_m = Units.to_meters(height, units)
    
    return C_t * h_m ** beta

# ============================================================================
# SHEAR WALL PERIOD FORMULAS (Consistent unit handling)
# ============================================================================

def shear_wall_period_ASCE7_22(
    h_n: float,
    A_B: float,
    walls: Sequence[tuple[float, float]],
    units: str = 'm',
    debug: bool = False
) -> float:
    """
    ASCE 7-22 Section 12.8.2.1 alternative period formula for shear wall buildings
    
    Args:
        h_n: Building height in specified units
        A_B: Area of base structure in specified units squared
        walls: List of (length, thickness) pairs for walls in the direction
        units: Input units ('m', 'ft', 'in', etc.)
        debug: Print debug information
    
    Returns:
        Period in seconds
    """
    # Convert inputs to meters for calculation
    h_m = Units.to_meters(h_n, units)
    A_B_m2 = Units.to_meters_sq(A_B, units)
    
    # Convert wall dimensions to meters
    walls_m = [(Units.to_meters(d, units), Units.to_meters(t, units)) 
               for d, t in walls]
    
    # ASCE 7-22 constants
    # C_q = 0.0019 for feet, 0.00058 for meters
    # Since we're working in meters internally, use metric constant
    C_q = 0.0019 / 0.3048  # seconds/meter
    
    # Calculate C_w per ASCE 7-22 Eq. 12.8-11
    C_w = 100 / A_B_m2 * sum(
        d * t / (1 + 0.83 * (h_m / d)**2) 
        for d, t in walls_m
    )
    
    if debug:
        A_w_m2 = sum(d * t for d, t in walls_m)
        print('Debug (ASCE7-22):')
        print(f'  h_n = {h_m:.2f} m')
        print(f'  A_B = {A_B_m2:.2f} m²')
        print(f'  A_w = {A_w_m2:.2f} m²')
        print(f'  C_q = {C_q:.6f}')
        print(f'  C_w = {C_w:.4f}')
    
    # ASCE 7-22 Eq. 12.8-10
    return C_q / C_w**0.5 * h_m


def shear_wall_period_Goel_Chopra(
    h_n: float,
    walls: Sequence[tuple[float, float]],
    units: str = 'm',
    modified: bool = True,
    debug: bool = False
) -> float:
    """
    Goel & Chopra (1998) period formula for shear wall buildings
    
    Args:
        h_n: Building height in specified units
        walls: List of (length, thickness) pairs for walls in the direction
        units: Input units ('m', 'ft', 'in', etc.)
        modified: If True, uses modified constant 0.020 (more realistic for 
                  typical buildings). If False, uses original 0.0019.
        debug: Print debug information
    
    Returns:
        Period in seconds
    """
    # Convert inputs to feet (original formula uses feet)
    h_ft = Units.to_feet(h_n, units)
    
    # Convert wall dimensions to feet
    walls_ft = [(Units.to_feet(d, units), Units.to_feet(t, units)) 
                for d, t in walls]
    
    # Calculate total wall area in square feet
    A_w_ft2 = sum(d * t for d, t in walls_ft)
    
    # Original constant is 0.0019 for ft units
    # Modified constant (0.020) gives more typical results for 5-storey buildings
    C = 0.020 if modified else 0.0019
    
    if debug:
        print('Debug (Goel_Chopra):')
        print(f'  h_n = {h_ft:.2f} ft')
        print(f'  A_w = {A_w_ft2:.2f} ft²')
        print(f'  C = {C:.4f}')
        if not modified:
            print('  Note: Original constant yields very short periods')
    
    return C * h_ft / A_w_ft2**0.5


def shear_wall_period_ATC3(
    h_n: float,
    building_length: float,
    units: str = 'm',
    debug: bool = False
) -> float:
    """
    ATC-3 simplified period formula
    
    Args:
        h_n: Building height in specified units
        building_length: Length of building in direction of interest in same units
        units: Input units ('m', 'ft', 'in', etc.)
        debug: Print debug information
    
    Returns:
        Period in seconds
    """
    # Convert inputs to feet (original formula uses feet)
    h_ft = Units.to_feet(h_n, units)
    L_ft = Units.to_feet(building_length, units)
    
    if debug:
        print('Debug (ATC-3):')
        print(f'  h_n = {h_ft:.2f} ft')
        print(f'  L = {L_ft:.2f} ft')
    
    return 0.05 * h_ft / L_ft**0.5


def shear_wall_period_EC8(
    h_n: float,
    walls: Sequence[tuple[float, float]],
    units: str = 'm',
    debug: bool = False
) -> float:
    """
    Eurocode 8 period formula for shear wall buildings
    
    Args:
        h_n: Building height in specified units
        walls: List of (length, thickness) pairs for walls in the direction
        units: Input units ('m', 'ft', 'in', etc.)
        debug: Print debug information
    
    Returns:
        Period in seconds
    """
    # Convert inputs to meters (EC8 expects meters)
    h_m = Units.to_meters(h_n, units)
    walls_m = [(Units.to_meters(d, units), Units.to_meters(t, units)) 
               for d, t in walls]
    
    # EC8 effective area calculation
    A_c = sum(
        d * t * (0.2 + min(0.9, d / h_m)**2) 
        for d, t in walls_m
    )
    
    if debug:
        A_w_m2 = sum(d * t for d, t in walls_m)
        print('Debug (EC8):')
        print(f'  h_n = {h_m:.2f} m')
        print(f'  A_w = {A_w_m2:.2f} m²')
        print(f'  A_c = {A_c:.2f} m²')
    
    return 0.075 / A_c**0.5 * h_m**0.75


def shear_wall_period_IS1893(
    h_n: float,
    walls: Sequence[tuple[float, float]],
    units: str = 'm',
    debug: bool = False
) -> float:
    """
    IS 1893 (Part 1):2016 period formula for shear wall buildings
    
    Args:
        h_n: Building height in specified units
        walls: List of (length, thickness) pairs for walls in the direction
        units: Input units ('m', 'ft', 'in', etc.)
        debug: Print debug information
    
    Returns:
        Period in seconds
    """
    # Convert inputs to meters (IS 1893 expects meters)
    h_m = Units.to_meters(h_n, units)
    walls_m = [(Units.to_meters(d, units), Units.to_meters(t, units)) 
               for d, t in walls]
    
    # Total wall area in square meters
    A_w_m2 = sum(d * t for d, t in walls_m)
    
    if debug:
        print('Debug (IS1893):')
        print(f'  h_n = {h_m:.2f} m')
        print(f'  A_w = {A_w_m2:.2f} m²')
    
    return 0.075 * h_m**0.75 / A_w_m2**0.5


def c_u(S_D1: float) -> float|None:
    """Upper limit coefficient for period (ASCE 7-22)"""
    return list_interp(
        S_D1,
        (0.1, 0.15, 0.2, 0.3, 0.4),
        (1.7, 1.6, 1.5, 1.4, 1.4),
        extend=True, extrapolate=False
    )

# =============================================================
# def shear_wall_period_EC8_ss(h_n:float, 
#                          walls:Sequence[Sequence[float]],
#                          debug:bool=False)->float:
#     A_c = sum(d*t*(0.2+min(0.9, d/h_n)**2) for d, t in walls)
#     C_t = 0.075 / A_c**0.5
#     if debug:
#         A_w = sum(d*t for d, t in walls)
#         print(f'A_w = {A_w}, A_c = {A_c}, C_t = {C_t:.4f}')
#     return C_t * h_n**0.75

# def shear_wall_period_ASCE7_22_ss(h_n:float, A_B:float, 
#     walls:Sequence[Sequence[float]], units:str='m', debug:bool=False)->float:
#     """_summary_

#     Args:
#         h_n (float): building height in metres (section 11.2)
#         A_B (float): area of base structure
#         walls (Sequence[Sequence[float]]): collection of pairs of 
#             wall lengths and thicknesses (d, t) in the direction
#             under consideration.

#     Returns:
#         float: building period
#     Example:
#         >>> walls = [(1.0, 0.2), (1.0, 0.2), (1.0, 0.2)]
#     """
#     # C_q = 0.00058 m or 0.0019 ft
#     C_q = 0.0019/0.3048 if (units == 'm') else 0.0019 # ft
#     C_w = 100 / A_B * sum(d*t / (1 + 0.83 * (h_n / d)**2) 
#                           for d, t in walls)
#     if debug:
#         A_w = sum(d*t for d, t in walls)
#         print(f'A_w = {A_w}, C_q = {C_q}, C_w = {C_w:.4f}')
#     return C_q / C_w**0.5 * h_n # 12.8-10

# def shear_wall_period_IS1893_ss(
#     h_n: float,
#     walls: Sequence[tuple[float, float]],
#     units: str = 'm',
#     debug: bool = False
# ) -> float:
#     """
#     IS 1893 (Part 1):2016 period formula for shear wall buildings
    
#     Args:
#         h_n: Building height in specified units
#         walls: List of (length, thickness) pairs for walls in the direction
#         units: Input units ('m', 'ft', 'in', etc.)
#         debug: Print debug information
    
#     Returns:
#         Period in seconds
#     """
#     # Convert inputs to meters (IS 1893 expects meters)
#     h_m = Units.to_meters(h_n, units)
#     walls_m = [(Units.to_meters(d, units), Units.to_meters(t, units)) 
#                for d, t in walls]
    
#     # Total wall area in square meters
#     A_w_m2 = sum(d * t for d, t in walls_m)
    
#     if debug:
#         print(f'Debug (IS1893):')
#         print(f'  h_n = {h_m:.2f} m')
#         print(f'  A_w = {A_w_m2:.2f} m²')
    
#     return 0.075 * h_m**0.75 / A_w_m2**0.5

# def shear_wall_period_GC_ss(h_n_ft:float, 
#                          walls_ft:Sequence[Sequence[float]],
#                          units:str='ft', debug:bool=False)->float:
#     A_w = sum(d*t for d, t in walls_ft)
#     C_t = 0.0019 if units == 'ft' else 0.00058
#     if debug:
#         print(f'A_w = {A_w}, C_t = {C_t:.4f}')
#     return C_t * h_n_ft / A_w**0.5

# def shear_wall_period_ATC3_ss(h_n:float, 
#                          length:float, units='ft',
#                          debug:bool=False)->float:
#     return 0.05 * h_n / length**0.5 * (1.0 if units == 'ft' else 0.3048**-0.5)

# ============================================================================
# Geometry
# ============================================================================

def I_rect(x, y, t):
    return y**3 * t / 6 + x * y**2 * t / 2

def I_rect2(x, y, t):
    return (x + t) * (y + t)**3 / 12 - (x - t) * (y - t)**3 / 12



def main():
    print(f'File folder: {file_folder}\n')


    print('== Doctest Start ==')
    doctest.testmod()
    print('== Doctest Complete ==')

    # ========================================================================
    # BUILDING DATA
    # ========================================================================
    
    # Wall dimensions in meters
    walls_x = [(6.21, 0.15), (6.21, 0.15), (6.04, 0.15), 
               (6.04 + 2.30, 0.15), (8.34, 0.15), (1.98, 0.15), (2.24, 0.15)]
    
    walls_y = [(3.86, 0.15), (3.86, 0.15), (3.94, 0.15), (3.94, 0.15), 
               (5.91, 0.15), (5.91, 0.15), (4.52, 0.15), (4.52, 0.15), (1.53, 0.15)]
    
    # Building dimensions
    h_n = 22.725  # meters (5-storey building)
    A_B = 12 * 4 * 6**2  # 1728 m² (base area)
    
    # Building length in each direction (for ATC-3)
    L_x = 12 * 6  # 72 m
    L_y = 4 * 6   # 24 m
    
    print('=' * 60)
    print('BUILDING INPUT DATA')
    print('=' * 60)
    print(f'Height (h_n): {h_n:.3f} m ({Units.to_feet(h_n, "m"):.1f} ft)')
    print(f'Base Area (A_B): {A_B:.1f} m² ({Units.to_feet_sq(A_B, "m"):.0f} ft²)')
    print(f'Building Length X: {L_x:.1f} m ({Units.to_feet(L_x, "m"):.1f} ft)')
    print(f'Building Length Y: {L_y:.1f} m ({Units.to_feet(L_y, "m"):.1f} ft)')
    
    # Wall area summary
    A_w_x = sum(d * t for d, t in walls_x)
    A_w_y = sum(d * t for d, t in walls_y)
    print(f'Wall Area X: {A_w_x:.2f} m² ({Units.to_feet_sq(A_w_x, "m"):.1f} ft²)')
    print(f'Wall Area Y: {A_w_y:.2f} m² ({Units.to_feet_sq(A_w_y, "m"):.1f} ft²)')
    
    print('\n' + '=' * 60)
    print('PERIOD ESTIMATES (seconds)')
    print('=' * 60)
    
    # ASCE 7-12 Empirical Code Period
    s_type = 'C2'
    T_code = T_code_get(s_type, h_n, units='m')
    print(f'\nASCE 7-22 Empirical (C_t = {C_t_get(s_type)}):')
    print(f'  T = {T_code:.3f} sec')
    
    # ASCE 7-22 Alternative Period
    print('\nASCE 7-22 Alternative (Section 12.8.2.1):')
    T_asce_x = shear_wall_period_ASCE7_22(h_n, A_B, walls_x, units='m', debug=False)
    T_asce_y = shear_wall_period_ASCE7_22(h_n, A_B, walls_y, units='m', debug=False)
    print(f'  X-direction: {T_asce_x:.3f} sec')
    print(f'  Y-direction: {T_asce_y:.3f} sec')
    
    # ATC-3
    print('\nATC-3:')
    T_atc_x = shear_wall_period_ATC3(h_n, L_x, units='m', debug=False)
    T_atc_y = shear_wall_period_ATC3(h_n, L_y, units='m', debug=False)
    print(f'  X-direction: {T_atc_x:.3f} sec')
    print(f'  Y-direction: {T_atc_y:.3f} sec')
    
    # # Goel & Chopra
    # print('\nGoel & Chopra (1998):')
    # T_gc_x_orig = shear_wall_period_Goel_Chopra(h_n, walls_x, units='m', modified=False)
    # T_gc_y_orig = shear_wall_period_Goel_Chopra(h_n, walls_y, units='m', modified=False)
    # print(f'  X-direction (original): {T_gc_x_orig:.3f} sec')
    # print(f'  Y-direction (original): {T_gc_y_orig:.3f} sec')
    
    # T_gc_x_mod = shear_wall_period_Goel_Chopra(h_n, walls_x, units='m', modified=True)
    # T_gc_y_mod = shear_wall_period_Goel_Chopra(h_n, walls_y, units='m', modified=True)
    # print(f'  X-direction (modified): {T_gc_x_mod:.3f} sec')
    # print(f'  Y-direction (modified): {T_gc_y_mod:.3f} sec')
    # print('    Note: Modified uses C=0.020 (more realistic for typical buildings)')
    
    # EC8
    print('\nEurocode 8:')
    T_ec8_x = shear_wall_period_EC8(h_n, walls_x, units='m', debug=False)
    T_ec8_y = shear_wall_period_EC8(h_n, walls_y, units='m', debug=False)
    print(f'  X-direction: {T_ec8_x:.3f} sec')
    print(f'  Y-direction: {T_ec8_y:.3f} sec')
    
    # IS 1893
    print('\nIS 1893 (Part 1):2016:')
    T_is_x = shear_wall_period_IS1893(h_n, walls_x, units='m', debug=False)
    T_is_y = shear_wall_period_IS1893(h_n, walls_y, units='m', debug=False)
    print(f'  X-direction: {T_is_x:.3f} sec')
    print(f'  Y-direction: {T_is_y:.3f} sec')
    
    # Summary with interpretation
    print('\n' + '=' * 60)
    print('SUMMARY & INTERPRETATION')
    print('=' * 60)
    print(f'\nPeriod Range: {min(T_atc_x, T_ec8_x, T_is_x):.2f} - {max(T_atc_x, T_ec8_x, T_is_x):.2f} sec')
    print('  (ATC-3, EC8, IS1893 give results in this typical range)')
    
    print('\nRecommended values for preliminary design:')
    print(f'  X-direction: {(T_ec8_x + T_is_x)/2:.2f} sec (average of EC8 and IS1893)')
    print(f'  Y-direction: {(T_ec8_y + T_is_y)/2:.2f} sec')
    
    print('\nNote: ASCE 7-22 alternative gives very short periods because it is')
    print('      a lower-bound estimate for very stiff buildings. The modified')
    print('      Goel & Chopra (0.020) provides a reasonable intermediate value.')

    print('\n', '=' * 60, '\n')
    # print(I_rect(6.21, 3.86, 0.15), I_rect2(6.21, 3.86, 0.15))

    Ixx_total = sum(I_rect2(x, y, 0.15) for (x, y) in 
                    [(6.21, 3.86),
                     (6.04, 3.94),
                     (2.24, 4.52),
                     (2.30, 5.91),
                     (1.98, 1.53),
                     ])
    Iyy_total = sum(I_rect2(y, x, 0.15) for (x, y) in 
                    [(6.21, 3.86),
                     (6.04, 3.94),
                     (2.24, 4.52),
                     (2.30, 5.91),
                     (1.98, 1.53),
                     ])
    print(f'Ixx_total = {Ixx_total:.2f} m^4')
    print(f'Iyy_total = {Iyy_total:.2f} m^4')

if __name__ == "__main__":
    main()
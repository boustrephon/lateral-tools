"""
lateral_tools — Seismic and lateral load engineering utilities.

This package provides tools for structural engineers working with:
    - Earthquake response spectra (GB50011, EC8, ASCE 7, ASCE 41)
    - Wind loading (HK Wind Code 2004/2019, EC1, MS1553)
    - Seismic hazard and risk analysis
    - Storey-based structural analysis
    - Ground motion processing and importing
    - Plotting (spectra, storey diagrams, tripartite, time-histories)

Sub-packages
------------
codes       : Code-specific seismic and wind loading functions
plotting    : Plotting functions for spectra, storey diagrams, EQ data
importing   : Ground motion data import (.asc format)

Modules
-------
math_utils   : Interpolation, type conversion, unit scaling
risk_utils   : Hazard recurrence, PGA estimation, remaining-life calculations
analysis     : Storey force distribution, shear/moment diagrams
print_utils  : Table formatting, storey data display, rebar detailing
"""
from .analysis import (
    shear_coords,
    shear_interp,
    split_shear_list,
    storey_elevations,
    storey_force_coeffs,
    storey_force_distribution,
    storey_forces,
    storey_shear_moment,
)
from .math_utils import (
    as_float,
    as_numeric,
    interp,
    list_interp,
    scale_to_SI_units,
    units_SI_scaling_dict,
)
from .print_utils import (
    import_storey_data,
    print_bar_areas_given_spacing,
    print_dir,
    print_storey_data,
    print_storey_force_data,
    print_table,
    printt,
)
from .risk_utils import (
    PGA2S475,
    GB_remaining_life_pga,
    S_calc,
    exceedance2return,
    remaining_life_pga,
    return2exceedance,
)

__all__ = [
    "as_float",
    "as_numeric",
    "exceedance2return",
    "GB_remaining_life_pga",
    "import_storey_data",
    "interp",
    "list_interp",
    "PGA2S475",
    "print_bar_areas_given_spacing",
    "print_dir",
    "print_storey_data",
    "print_storey_force_data",
    "print_table",
    "printt",
    "remaining_life_pga",
    "return2exceedance",
    "S_calc",
    "scale_to_SI_units",
    "shear_coords",
    "shear_interp",
    "split_shear_list",
    "storey_elevations",
    "storey_force_coeffs",
    "storey_force_distribution",
    "storey_forces",
    "storey_shear_moment",
    "units_SI_scaling_dict",
]

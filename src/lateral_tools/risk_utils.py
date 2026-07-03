"""
risk_utils.py — Seismic hazard and risk analysis functions.

Part of the lateral-tools package for seismic and lateral load engineering.

This module provides:
- Hazard recurrence: exceedance2return, return2exceedance
- Spectral acceleration: PGA2S475, S_calc
- Remaining-life PGA estimates: remaining_life_pga, GB_remaining_life_pga
"""
from __future__ import annotations
from math import exp, log
from typing import Any
def exceedance2return(exceedance_probability:float, period:float=50)->float:
    """Returns the return period for ground motion with the given probability of
    being exceeded over the period considered. For example, the return period 
    for a ground motion with a 10% probability of being exceeded in 50 years
    is around 475 years. 
  
    Args:
        exceedance_probability (float): _description_
        period (float, optional): _description_. Defaults to 50.

    Returns:
        float: _description_
    
    Examples:
        >>> isclose(exceedance2return(exceedance_probability=0.1, period=50), 474.56107905149526, rel_tol=1e-9)
        True
        >>> isclose(exceedance2return(exceedance_probability=0.2, period=50), 224.07100588622754, rel_tol=1e-9)
        True
  

    """
    return - period / log(1 - exceedance_probability)




def return2exceedance(return_period, period=50):
    """Returns the probability of a ground motion being exceeded over a period  
    when given a ground motion return period. For example, the probability of a  
    ground motion with return period of 475 years being exceeded over a time  
    period of 50 years is around 10%. 
    >>> isclose(return2exceedance(return_period=475, period=50), 0.09991237374774076)
    True
    >>> isclose(return2exceedance(return_period=225, period=50), 0.19926259708319194)
    True
    """
    return 1 - exp(-period/return_period)




def PGA2S475(a_g_475):
    """"Returns a dictionary containing S_s & S_1 values 
    from a reference 475-year return PGA.
    
    Reference: 15WCEE, 'Deriving S_s and S_1 parameters 
        from PGA maps'
    """
    def S_s_calc(PGA): return 0.3386 * PGA + 2.1696
    def S_1_calc(PGA): return 0.5776 * PGA + 0.5967
    return {'S_s':S_s_calc(a_g_475), 'S_1':S_1_calc(a_g_475)}




def S_calc(S_s_475, S_1_475, T):
    """"Returns a dictionary containing S_s, S_1, 
    k_s & k_1 values for a given return period based 
    on a reference PGA and reference return period.
    
    Reference: 15WCEE, 'Deriving S_s and S_1 parameters 
        from PGA maps'
    """
    def S_func(S_475, k, T): 
        return S_475 * (T / 475)**k
    
    S_s_475_S_keys = (0.20, 0.55, 0.90, 1.20, 1.50)
    S_s_475_k_vals = (0.46, 0.42, 0.39, 0.35, 0.30)
    S_1_475_S_keys = (0.01, 0.02, 0.05, 0.10, 0.20)
    S_1_475_k_vals = (0.40, 0.40, 0.40, 0.40, 0.40)
    
    k_s = list_interp(S_s_475, S_s_475_S_keys, S_s_475_k_vals)
    S_s = S_func(S_s_475, k_s, T)
    
    k_1 = list_interp(S_1_475, S_1_475_S_keys, S_1_475_k_vals)
    S_1 = S_func(S_1_475, k_1, T)
    
    return {'S_s':S_s, 'S_1':S_1, 'k_s':k_s, 'k_1':k_1}




def remaining_life_pga(pga_50:float, years:float, 
                   k:float=0.46, basis_years:float=50)->float:
    """Returns the remaining life of a building given a 50-year PGA and a  
    number of years.

    Optionally can specify a different design life basis with the associated PGA.

    Args:
        pga_50 (float): _description_
        years (float): _description_
        k (float, optional): _description_. Defaults to 0.46.
        basis_years (float, optional): _description_. Defaults to 50.

    Returns:
        float: _description_
    """

    return pga_50 * (years/basis_years)**k



def GB_remaining_life_pga(intensity:str, soil_class:str, 
                          years_list:Sequence[float], k:float=0.46, 
                          as_df:bool=False, basis_years:float=50)->dict[str,float]|Any:
    """Returns the remaining life of a building given a 50-year PGA and a  
    number of years.

    Optionally can specify a different design life basis with the associated PGA.    

    """
    pga_50 = 0.15 # g
    pga_50 = 0.09
    pga_freq_50 = pga_50 / 0.15 * 0.054 / 0.45
    pga_des_50 = pga_50 / 0.15 * 0.15 / 0.45
    pga_rare_50 = pga_50 / 0.15 * 0.324 / 0.45

    if len(years_list) < 1:
        years_list = [10, 20, 30, 40, 50]
    
    output_dict = {}
    ratios = [(years / basis_years)**k for years in years_list]
    output_dict['ratios'] = ratios
    output_dict['pga_freq'] = [pga_freq_50 * ratio for ratio in ratios]
    output_dict['pga_des']  = [pga_des_50 * ratio for ratio in ratios]
    output_dict['pga_rare'] = [pga_rare_50 * ratio for ratio in ratios]
    output_dict['intensities'] = [get_intensity(pga_des_50 * ratio) for ratio in ratios]

    if as_df:
        import pandas as pd
        df = pd.DataFrame({
            "Remaining Life (years)": years_list,
            "Ratio": [f"{100 * x:.1f} %" for x in ratios],
            "Frequent Pk (g)": [f"{x:.4f}" for x in output_dict['pga_freq']],
            "Design Pk (g)":   [f"{x:.4f}" for x in output_dict['pga_des']],
            "Rare Pk (g)":     [f"{x:.4f}" for x in output_dict['pga_rare']],
            "Intensity":        [f"{x:s}" for x in output_dict['intensities']],
        })
        return df

    return output_dict

# =======================================
# ===== Storey-based Functions ==========
# =======================================

# Functions:
# storey_force_coeffs, storey_force_distribution, storey_forces, storey_shear_moment, shear_coords
# storey_force_coeffs, storey_force_distribution, storey_shear_moment, shear_coords




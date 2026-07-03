"""
_summary_
"""


def deep_beam_shear(k_1: float, f_cd: float, 
                    b_1: float, gamma_b: float):
    """Design Shear Capacity Calculation (Eq 5.17)

    Design shear capacity of deep beam carried by concrete of one-way action, kg

    Args:
        k_1 (float): coefficient for one-way deep-beam shear
        f_cd (float): design compressive strength of concrete, kg/cm2
        b_1 (float): width of the cross section, cm
        gamma_b (float): partial safety factor for member

    Returns:
        float: design shear capacity
    """
    V_cdb = k_1 * (f_cd)**0.5 * b_1 * gamma_b
    return V_cdb

def deep_beam_shear_coefficient():
    pass


"""Malaysian Wind Code MS 1553 - 2002

The following are wrong(!)
Terrain type 1: Mz,cat = 1.0 + 8.0 (z/10)**0.13
Terrain type 2: Mz,cat = 1.0 + 7.0 (z/10)**0.2
Terrain type 3: Mz,cat = 1.0 + 6.0 (z/10)**0.07
Terrain type 4: Mz,cat = 1.0 + 5.0 (z/10)**0.3
"""


def get_V_sit(V_s, M_d, M_zcat, M_s = 1, M_h = 1):
    """Site wind speed
    
    Args:
        V_s: 33.5 (Zone 1) 32.5 (Zone 2) for 50 year return
        M_d: direction factor (default 1.0)
        M_zcat: terrain / height multiplier (section 4)
        M_s: shielding multiplier (section 4)
        
    Returns:
        V_sit: site wind speed (multiply by I for V_des)
    """
    V_sit = V_s * M_d * M_zcat * M_s * M_h
    return V_sit


def get_V_des(V_sit:float, I:float)-> float:
    """Design force for given site wind speed and project importance

    Args:
        V_sit (float): _description_
        I (float): _description_

    Returns:
        float: _description_
    """
    V_des = I * V_sit
    return V_des


def get_pressure(V_des, C_fig, C_dyn=1, rho_air = 1.225):
    """Design wind pressure
    
    Args:
        V_des: for 50 year return
        C_fig: aerodynamic shape factor (section 5)
        C_dyn: dynamic response factor (section 6)
        rho_air: 1.225 kg/m3
        
    Returns:
        p: design wind pressure (Pa)
    """
    p = 0.5 * rho_air * V_des ** 2 * C_fig * C_dyn
    return p 
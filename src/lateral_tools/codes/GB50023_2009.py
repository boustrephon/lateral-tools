""" 

"""
import numpy as np

# Remaining Service Life
rsl_dict = {'A': 30, 'B': 40, 'C': 50,}

# GB 50223-2008 A, B, C, and D are the abbreviated terms for the special
# fortification category, key fortification category, standard fortification category, and moderate
# fortification category

seismic_fortification_category_dict = {
'A': {'Name': 'Special Fortification Category (Category A / 甲类)',
'Definition': '''Buildings with special facilities, major projects involving national
public security, and those where severe secondary disasters (e.g., radioactive
leaks, massive explosions) may occur.''',
'Seismic Standard': '''Highest level. Seismic details (detailing and construction requirements) are enhanced by one degree higher than the local intensity (or more for Intensity 9). Seismic actions (design forces) are determined by a site-specific seismic safety evaluation higher than the local intensity.'''},
'B': {'Name': 'Key Fortification Category (Category B / 乙类)', 
'Definition': '''Lifeline buildings whose function must not be interrupted (or must be quickly restored), and buildings where failure could lead to mass casualties.''',
'Seismic Standard': '''Enhanced level. Seismic details are enhanced by one degree (or more for Intensity 9). Seismic actions are based on the local intensity.
This is the most common category for critical power infrastructure.'''},
'C': {'Name': 'Standard Fortification Category (Category C / 丙类)',
'Definition': '''The vast majority of buildings not falling into the other categories.''',
'Seismic Standard': '''Standard level. Both seismic details and actions are based on the local seismic fortification intensity. The objective is “no collapse under the rare earthquake (MCE).”'''},
'D':  {'Name': 'Moderate Fortification Category (Category D / 丁类)', 
'Definition': '''Buildings with few occupants where damage is unlikely to cause secondary disasters.''',
'Seismic Standard': '''Reduced level. Seismic details may be reduced by one degree (except for Intensity 6). '
'Seismic actions are based on the local intensity.'''},
}

seismic_classification_dict = {
    'A': 'Special',
    'B': 'Key',
    'C': 'Standard',
    'D': 'Moderate',
    'E': 'No Classification'
    }  


site_class_dict = {
    'Group 1': {'I0':0.20, 'I1': 0.25, 'II': 0.35, 'III': 0.45, 'IV':0.65},
    'Group 2': {'I0':0.25, 'I1': 0.30, 'II': 0.40, 'III': 0.55, 'IV':0.75},
    'Group 3': {'I0':0.30, 'I1': 0.35, 'II': 0.45, 'III': 0.65, 'IV':0.90},
}

# Table 5.1.4-2
T_G_dict_GB50011_2001 = {
    1:{'I': 0.25, 'II': 0.35, 'III': 0.45, 'IV': 0.65},
    2:{'I': 0.30, 'II': 0.40, 'III': 0.55, 'IV': 0.75},
    3:{'I': 0.35, 'II': 0.45, 'III': 0.65, 'IV': 0.90}
    }

# Table 5.1.4-2
T_G_dict_GB50011_2010 = {
    1:{'I0': 0.20, 'I1': 0.25, 'II': 0.35, 'III': 0.45, 'IV': 0.65},
    2:{'I0': 0.25, 'I1': 0.30, 'II': 0.40, 'III': 0.55, 'IV': 0.75},
    3:{'I0': 0.30, 'I1': 0.35, 'II': 0.45, 'III': 0.65, 'IV': 0.90}
    }

# Table 3.0.5
T_G_dict_GB50023_2009 = {
    1:{'I': 0.20, 'II': 0.30, 'III': 0.40, 'IV': 0.65},
    2:{'I': 0.20, 'II': 0.30, 'III': 0.40, 'IV': 0.65},
    3:{'I': 0.25, 'II': 0.40, 'III': 0.55, 'IV': 0.85}
}


def GB_spectrum(alpha_max, T, 
                site_class='I0', group=1, 
                zeta=0.05, T_G_dict=T_G_dict_GB50023_2009):
    """
    Returns spectral acceleration Sa(g) for given PGA(g) and period array T(s).
    Eurocode 8 Type 1 spectrum, ground type B.
    
    Args:
        alpha_max (_type_): _description_
        T (_type_): _description_
        site_class (str, optional): _description_. Defaults to 'I0'.
        group (str, optional): _description_. Defaults to 'Group 1'.
        zeta (float, optional): _description_. Defaults to 0.05.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    gamma = 0.9 + (0.05 - zeta) / (0.3 + 6.0 * zeta)
    eta_1 = max(0.0, 0.02 + (0.05 - zeta) / (4.0 + 32.0 * zeta))
    eta_2 = max(0.55, 1.0 + (0.05 - zeta) / (0.08 + 1.6 * zeta))
    
    T_g = T_G_dict[group][site_class]
    # pga = alpha_max * 0.45
    Sa = np.zeros_like(T)
    for i, t in enumerate(T):
        if t <= 0.1:
            Sa[i] = (0.45 + (eta_2 - 0.45) * 10 * t) * alpha_max
        elif t <= T_g:
            Sa[i] = eta_2 * alpha_max
        elif t <= 5 * T_g:
            Sa[i] = (T_g / t)**gamma * eta_2 * alpha_max
        elif t <= 6.0:
            Sa[i] = (eta_2 * 0.2**gamma - eta_1 * (t - 5 * T_g)) * alpha_max
        else:
            raise ValueError(f'Period ({t:.4f}) is greater than 6.0')
    return Sa


"""
These are some old functions that are mainly based on the2001 version of the code, 
since there is a readily available traonslation.
"""



# China Code
# Table 5.1.4-1
seismic_zone_dict = {
    '6': {'intensity':'6 (0.05g)', 'PGA': 0.05, 'alphamax_L1': 0.04, 'alphamax_L3': None}, 
    '7': {'intensity':'7 (0.10g)', 'PGA': 0.10, 'alphamax_L1': 0.08, 'alphamax_L3': 0.5},
    '7+': {'intensity':'7 (0.15g)', 'PGA': 0.05, 'alphamax_L1': 0.12, 'alphamax_L3': 0.72},
    '8': {'intensity':'8 (0.20g)', 'PGA': 0.05, 'alphamax_L1': 0.16, 'alphamax_L3': 0.9},
    '8+': {'intensity':'8 (0.30g)', 'PGA': 0.05, 'alphamax_L1': 0.24, 'alphamax_L3': 1.2},
    '9': {'intensity':'9 (0.40g)', 'PGA': 0.05, 'alphamax_L1': 0.32, 'alphamax_L3': 1.4},
    }

fortification_dict = {
    '6':0.05, '7':0.10, '7+':0.15, '8':0.20, '8+':0.30, '9':0.40, 
}
#def 


T_G_dict = {
    1:{'I':0.25,'II':0.35,'III':0.45,'IV':0.65},
    2:{'I':0.30,'II':0.40,'III':0.55,'IV':0.75},
    3:{'I':0.35,'II':0.45,'III':0.65,'IV':0.90}
    }

def T_G_freq_calc(MD_Group, soil_class):
  """GB50011-2001 Table 5.1.4-2
  >>> T_G_freq_calc(2, 'IV')
  0.75
  """
  return T_G_dict.get(MD_Group,{}).get(soil_class, None)


def T_G_rare(MD_Group, soil_class, T):
  return T_G_freq_calc(MD_Group, soil_class) + (0.05 if T > 2 else 0.05)


alpha_max_dict = {
    1: {'6':18, '7':35, '7+':55, '8':70, '8+':110, '9':140, },
    #2: {'6':0.04, '7':None, '7+':None, '8':None, '8+':None, '9':None, },
    3: {'6':None, '7':220, '7+':310, '8':400, '8+':510, '9':620, },
}


alpha_dict = {
    1: {'6':0.04, '7':0.08, '7+':0.12, '8':0.16, '8+':0.24, '9':0.32, },
    #2: {'6':0.04, '7':0.08, '7+':0.12, '8':0.16, '8+':0.24, '9':0.32, },
    3: {'6':None, '7':0.50, '7+':0.72, '8':0.90, '8+':1.2, '9':1.4, },
}


def alpha_calc(zeta = 0.05):
    # 5.1.5-1
    gamma = 0.9 + (0.05 - zeta) / (0.5 + 5 * zeta)
    # 5.1.5-2
    eta_1 = max(0, 0.02 + 0.125 * (0.05 - zeta))
    # 5.1.5-3
    eta_2 = max(0.55, 1 + (0.05 - zeta) / (0.06 + 1.7 * zeta))
    
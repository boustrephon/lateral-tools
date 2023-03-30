"""
Structural Types according to ASCE41-06

ASCE/SEI 41-06, Seismic Rehabilitation of Existing Buildings

*** INCOMPLETE ***

The full list of codes referenced may be found in the commentary

"""

# ASCE41 ch 4
# M_s_dict = {'C1': {'CP':2.0, 'LS': 1.5, 'IO': 1.0},
#             'C2': {'CP':4.5, 'LS': 3.0, 'IO': 1.5},
#             'URM': {'CP':1.75, 'LS': 1.25, 'IO': 1.0}}



struct_types_dict = {
'W1': 'Light Wood Frame',
'W1A': 'Multistory Multi-Unit Residential Wood Frame',
'W2': 'Commercial and Industrial Wood Frame',
'S1': 'Steel MRF',
'S2': 'Steel Braced Frame',
'S2A': 'Steel Braced Frame Flexible Diaphragm',
'S3': 'Steel Light Frame',
'S4': 'Steel Frame RC SW',
'S5': 'Steel Frame URM INF',
'S5A': 'Steel Frame URM INF Flexible Diaphragm',
'C1': 'RC MRF',
'C2': 'RC Shear Wall',
'C2A': 'RC Shear Wall Flexible Diaphragm',
'C3': 'RC Frame URM Infill',
'C3A': 'RC Frame URM Infill Flexible Diaphragm',
'PC1': 'Precast/Tilt-Up RC SW',
'PC1A': 'Precast/Tilt-Up RC SW Flexible Diaphragm',
'PC2': 'Precast RC Frame with Shear Walls',
'PC2A': 'Precast RC Frame without Shear Walls',
'RM1': 'Reinforced Masonry Bearing Walls',
'RM2': 'Reinforced Masonry Bearing Walls Flexible Diaphragm',
'URM': 'Unreinforced Masonry Bearing Walls',
'URM': 'Unreinforced Masonry Bearing Walls Flexible Diaphragm'
}


# Simplified Rehabilitation
# Table 10-1  ASCE41-06 Limitations on Use of the Simplified Rehabilitation Method (number of storeys)
struct_types_full_dict = {
'W1': {'Description':'Light Wood Frame', 'L':3, 'M':3, 'H':2, },
'W1A': {'Description':'Multistory Multi-Unit Residential Wood Frame', 'L':3, 'M':3, 'H':2, },
'W2': {'Description':'Commercial and Industrial Wood Frame', 'L':3, 'M':3, 'H':2, },
'S1': {'Description':'Steel MRF', 'L':6, 'M':4, 'H':3, },
'S1': {'Description':'Steel MRF', 'L':4, 'M':4, 'H':3, },
'S2': {'Description':'Steel Braced Frame', 'L':6, 'M':4, 'H':3, },
'S2A': {'Description':'Steel Braced Frame Flexible Diaphragm', 'L':3, 'M':3, 'H':3, },
'S3': {'Description':'Steel Light Frame', 'L':2, 'M':2, 'H':2, },
'S4': {'Description':'Steel Frame RC SW', 'L':6, 'M':4, 'H':3, },
'S5': {'Description':'Steel Frame URM INF', 'L':3, 'M':3, 'H':0, },
'S5A': {'Description':'Steel Frame URM INF Flexible Diaphragm', 'L':3, 'M':3, 'H':0},
'C1': {'Description':'RC MRF', 'L':3, 'M':0, 'H':0},
'C2': {'Description':'RC Shear Wall', 'L':6, 'M':4, 'H':3},
'C2A': {'Description':'RC Shear Wall Flexible Diaphragm', 'L':3, 'M':3, 'H':3},
'C3': {'Description':'RC Frame URM Infill', 'L':3, 'M':0, 'H':0},
'C3A': {'Description':'RC Frame URM Infill Flexible Diaphragm', 'L':3, 'M':0, 'H':0},
'PC1': {'Description':'Precast/Tilt-Up RC SW', 'L':3, 'M':2, 'H':2},
'PC1A': {'Description':'Precast/Tilt-Up RC SW Flexible Diaphragm', 'L':3, 'M':2, 'H':2},
'PC2': {'Description':'Precast RC Frame with Shear Walls', 'L':3, 'M':2, 'H':0},
'PC2A': {'Description':'Precast RC Frame without Shear Walls', 'L':0, 'M':0, 'H':0},
'RM1': {'Description':'Reinforced Masonry Bearing Walls', 'L':3, 'M':3, 'H':3},
'RM2': {'Description':'Reinforced Masonry Bearing Walls Flexible Diaphragm', 'L':6, 'M':4, 'H':3},
'URM': {'Description':'Unreinforced Masonry Bearing Walls', 'L':3, 'M':3, 'H':2},
'URMA': {'Description':'Unreinforced Masonry Bearing Walls Flexible Diaphragm', 'L':3, 'M':3, 'H':2}
}

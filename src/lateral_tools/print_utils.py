"""
print_utils.py — Display, printing, and data import functions.

Part of the lateral-tools package for seismic and lateral load engineering.

This module provides:
- Package inspection: print_dir
- Table formatting: printt, print_table
- Storey data display: print_storey_data, print_storey_force_data
- Data import: import_storey_data
- Rebar detailing: print_bar_areas_given_spacing
"""
from __future__ import annotations
from csv import reader
from tkinter import filedialog, Tk
from os.path import exists
from typing import Any, Sequence
def print_dir(pkg, n=4):
    """Print all non-private members of a package, n items per line.

    Useful for quickly seeing what members are in a package, without having to
    manually open the package and look.

    Parameters
    ----------
    pkg : package
        The package to print members of.
    n : int, optional
        Number of items per line. Default is 4.
    """
    
    [print(x) if i % n == (n - 1) else print(x, end=' | ') 
     for i, x in enumerate(dir(pkg)) if not x.startswith('__')]


# def as_float(x):
#     try:
#         return float(x)        
#     except ValueError:
#         return x




def printt(table, h1 = 20, h2 = 12):
    """
    Printing tables where first column is text and the remainder are floats
    h1 is the 
    """
    for line in table:
        print('{:20s}'.format(line[0]), end = ' ')
        [print('{:8.4f}'.format(item), end = ' ') for item in line[1:]]
    return None





def import_storey_data(file_path: str|None, delimiter:str =','):
    """
    Importing 4-column data from CSV file
    First column is assumed to be keys (storey names) and
    other three are required to be numeric. Any row containing
    non-numeric data in the last three columns will be ignored
    """
    if file_path is None:
        root = Tk()
        root.withdraw()  # Hide the main window
        file_path = filedialog.askopenfilename(
            initialdir = "/", title = "Select file", 
            filetypes = (("csv files","*.csv"),("all files","*.*"))
        )
        root.destroy()  # Destroy the main window

    print(f'Building data file path is "{file_path}"')
    print('File EXISTS' if exists(file_path) else 'File NOT found')

    with open(file_path, newline='') as csvfile:
        csv_reader = reader(csvfile, dialect = 'excel', delimiter=delimiter)
        data = [r for r in csv_reader]
    
    # _, b, c, d = data[0]
    # storey_headers = data[0] if not all([isinstance(as_float(x), float) 
    #                                         for x in (b, c, d)]) else []
    
    storey_data = [[a, float(b), float(c), float(d)] 
                for a,b,c,d in data if all([isinstance(as_float(x), float) 
                                            for x in (b, c, d)])]
    
    return storey_data




def print_table(table, h1 = 20, h2 = 12, h3 = 4):
    """
    Printing tables where first column is text and the remainder are floats
    h1 is the width of the text column, 
    h2 is the width of all other columns
    h3 is the number of decimal places
    """
    for line in table:
        print(f'{line[0]:{h1}s}', end = ' ')
        [print(f'{item:{h2-h3}.{h3}f}', end = ' ') for item in line[1:]]
    return None




def print_storey_data(
        storey_elevs_m, 
        storey_weights_kN, 
        storey_names = [], 
        base_level_m = 0,
        ascending = True):
    """Prints out storey data
    All data should be for ascending height. The base level should 
    not be included in the data, except as the `base_level_m` value.
    
    """
    if ascending:
        storey_elevs_m.reverse()
        storey_weights_kN.reverse()
        storey_names.reverse()
    
    storey_heights_m = [h2 - h1 for h2, h1 in zip(storey_elevs_m, list(storey_elevs_m[1:]) + [base_level_m])]
    storey_masses_kg = [1000 * w / 9.80665 for w in storey_weights_kN] # convert kN to kg

    if len(storey_weights_kN) == len(storey_elevs_m) == len(storey_names):
        pass
    else:
        storey_names = [f'L{i:04d}' for i, _ in enumerate(storey_elevs_m)]
        storey_names.reverse()
        
    basic_f = 'Height is {:.2f}m, total mass is {:8,.0f} tonne, total weight is {:8,.0f} kN\n'
    print(basic_f.format(max(storey_elevs_m), 0.001 * sum(storey_masses_kg), sum(storey_weights_kN)))

    tit_f1 = '{:20s} {:9s} {:10s} {:10s} {:10s}'
    tit_f2 = '{:20s} {:9s} {:10s} {:>10,.0f} {:8,.0f}'
    print(tit_f1.format('Storey', 'Ht. (m)', 'Elev (m)', 'Mass (kg)', 'Weight (kN)'))
    num_f = '{:20s} {:6.3f}   {:8.3f}   {:>10,.0f}   {:8,.0f}'
    [print(num_f.format(n, h, e, m, w)) for n, h, e, m, w in zip(storey_names, storey_heights_m, 
                                                    storey_elevs_m, storey_masses_kg, storey_weights_kN)]
    # base_level = f'{base_level_m:8.3f}'
    # print(base_level)
    print(tit_f2.format('BASE', '------', f'{base_level_m:7.3f}  ', sum(storey_masses_kg), sum(storey_weights_kN)))
    print()

    return None




def print_storey_force_data(
        V_kN, 
        storey_elevs_m, 
        storey_weights_kN, 
        storey_names = [], 
        k = 1,
        base_level_m = 0,
        ascending = True
        ):
    """Calculates and prints out storey force data
    All data should be for ascending height. The base level should 
    not be included in the data, except as the `base_level_m` value.
    
    """
    if ascending:
        list(storey_elevs_m).reverse()
        list(storey_weights_kN).reverse()
        list(storey_names).reverse()
    
    storey_heights_m = [h2 - h1 for h2, h1 in zip(storey_elevs_m, list(storey_elevs_m[1:]) + [base_level_m])]
    # storey_masses_kg = [1000 * w / 9.80665 for w in storey_weights_kN] # convert kN to kg

    if len(storey_weights_kN) == len(storey_elevs_m) == len(storey_names):
        pass
    else:
        storey_names = [f'L{i:04d}' for i, _ in enumerate(storey_elevs_m)]
        storey_names.reverse()

    k = k
    storey_force_coeff_list = storey_force_coeffs(storey_elevs_m, storey_weights_kN, k)
    storey_loads_kN = storey_force_distribution(V_kN, storey_elevs_m, storey_weights_kN, k)

    storey_shears_kN, storey_moments_kNm = storey_shear_moment(storey_loads_kN, storey_heights_m, False)

    tit_f1 = '{:20s} {:9s} {:10s} {:10s} {:10s} {:10s} {:10s}'
    print(tit_f1.format('Storey', 'Elevs (m)', 'Weights (kN)', ' Coeffs ', 'Forces (kN)', 'Shears (kN)', 'Moments (kNm)'))
    num_f = '{:20s} {:8.2f}   {:8,.0f}   {:8.4f}   {:8,.0f}   {:>10,.0f}   {:>10,.0f}'
    [print(num_f.format(n, h, wt, c, f, s, m)) for n, h, wt, c, f, s, m in zip(storey_names, storey_elevs_m, storey_weights_kN, 
                            storey_force_coeff_list, storey_loads_kN, storey_shears_kN, [0] + storey_moments_kNm[:-1])]
    tit_f2 = '{:20s} {:8.2f} {:10s} {:12s} {:10s} {:12s} {:>10,.0f}'
    print(tit_f2.format('BASE', base_level_m, '    ------', '    ------', ' ------   ', '  ------  ', storey_moments_kNm[-1]))
    tit_f3 = '{:20s} {:10s} {:8,.0f} {:10.4f} {:10,.0f} {:14s} {:12s}'
    print(tit_f3.format('SUM ', '  ------', sum(storey_weights_kN), sum(storey_force_coeff_list), sum(storey_loads_kN), '     ------', '     ------'))
    
    return None




def print_bar_areas_given_spacing(bar_spacing_mm, 
                       bar_list = [16, 20, 25, 32, 36, 40]):
    
    [print(f'T{D}@200 = {1000 * D**2 / 4 / bar_spacing_mm}mm2/m') 
        for D in bar_list]
    




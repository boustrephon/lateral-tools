"""
Functions for importing ground motion data

The only functions currently provided are for a format provided from the recent 2023 Turkey EQ

Other resources for future reference:
https://ds.iris.edu/ds/nodes/dmc/data/formats/seed/
https://www.geonet.org.nz/data/supplementary/strong_motion_file_formats
http://www.cosmos-eq.org/menu/2_Publications/cosmos_format_1_20.pdf
https://www.usgs.gov/programs/earthquake-hazards/science/seismic-design-data-sets-text-file-format
https://github.com/pham-hung/PEER-to-ASCII
https://peer.berkeley.edu/peer-strong-ground-motion-databases


"""

from typing import Union, Any, Dict, List, Tuple # , PathLike
import doctest
from sys import path as sys_path
from os.path import abspath, join as path_join, dirname, exists, basename, splitext

# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))


from pySeismic.utilities import scale_to_SI_units, as_float #, as_numeric
import csv


def read_csv(csv_file, header=1, as_dict=True):
    with open(csv_file, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        data = [row if i < header else [float(r) for r in row] for i, row in enumerate(csv_reader)]
    axes_txt = ['x', 'y', 'z']
    if as_dict:
        x_data, y_data = zip(*data)
        x_key = 'x' if header == 0 else (x_data[0] if header == 1 else tuple(x_data[:header]))
        y_key = 'y' if header == 0 else (y_data[0] if header == 1 else tuple(y_data[:header]))
        return {x_key: x_data[header:], y_key: y_data[header:]}
    else:
        return data


def ASC_record_pair_SI_scaling(
        record: Union[
            List[Union[float, int]], 
            Tuple[Union[float, int]]], 
        units: Union[str, float, int]
        ) -> Tuple[float]:
    """Scales a time, value pair to SI units.

    This scales a time, value pair to SI units based
    on the value of `units`. If `units` is numeric, then
    second item is scaled by that number, if `units` is
    text then a scaling factor to SI is looked up.

    Args:
        record: a data pair where the second needs scaling
        units: either 
    
    Result:
        tuple containing the data pair where second has been 
        scaled to SI
    """
    x, y = record
    if isinstance(units, str):
        y_scaled = scale_to_SI_units(y, units)
    elif isinstance(units, (float, int)):
        y_scaled = y * units
    return x, y_scaled


def read_ASC(
        file_path: str, # Union[str, PathLike]
        seed_dict: Dict[str, Any]={}
        ) -> Dict[str, Any]:
    """Imports ASC seismic data into a dictionary
    
    Imports an ASC files and returns a dictionary 
    containing all the meta data (on the basis that the 
    data is compatible with DYNA 1.2 header format). 
    It places the time-domain data or frequency-domain 
    data under the "Data" keyword, and also adds 
    'xlabel' and 'ylabel' data for plotting.

    The dictionary can easily be read into a Pandas dataframe for
    further processing.

    Args:
        file_path: path to the ASC file
        seed_dict: additional metadata for the dictionary

    Returns:
        asc_dict: 
    """
    if exists(file_path):
        t = 0
        th_list = []
        asc_dict = seed_dict
        asc_dict['Filename'] = basename(file_path)
        asc_dict['Filepath'] = file_path
        with open(file_path, 'r') as asc_file:
            sampling_interval = 0.0
            # header_format = ''
            data_type = ''
            for line in asc_file:
                if len(line.strip()) == 0:   #  empty line
                    pass
                elif ':' in line:
                    key, value = line.split(':', 1)
                    asc_dict[key.strip()] = value.strip()
                    if key.strip() == 'SAMPLING_INTERVAL_S':
                        sampling_length = len(value.strip())
                        sampling_interval = as_float(value.strip(), 0.0)
                    elif key.strip() == 'HEADER_FORMAT':
                        # header_format = value.strip() # DYNA 1.2
                        pass
                    elif key.strip() == 'DATA_TYPE':
                        data_type = value.strip() # e.g. Acceleration Response Spectrum
                    elif key.strip() == 'UNITS':
                        units = value.strip() # e.g. cm/s^2
                        scale = scale_to_SI_units(units, 1)
                else:
                    data = [as_float(n, 'None') for n in line.strip().split()]
                    data = [(None if a=='None' else a) for a in data]
                    if sampling_interval <= 0:
                        pass
                    elif len(data) == 1:
                        th_list.append((t, scale * as_float(line.strip())))
                        t = round(t + sampling_interval, sampling_length)
                    elif len(data) == 2:
                        th_list.append((float(data[0]), scale * float(data[1])))
        asc_dict['Data'] = th_list
        if len(data_type) > 0:
            asc_dict['xlabel'] = 'Period (s)' if 'Spectrum' in data_type else 'Time (s)'
            asc_dict['ylabel'] = f'{data_type.split()[0]} ({units})'
        return asc_dict
    
    else:
        return None


def ASC_import(asc_files, folder=r'./', 
               export=False, debug=False):
    """
    Imports a list of ASC files from a single folder and returns 
    a dictionary containing one entry for each file. In other words, 
    it does not identify pairs or triplets of files. 

    The dictionary can easily be read into a Pandas dataframe for
    further processing.

     
    """
    EQ_dict = {}
    for asc in asc_files:
        if basename(asc) != asc:
            raise ValueError('List of ASC files should be basenames only')
        name, _ = splitext(asc)
        asc_list = name.strip().split('_')
        event = asc_list[0].strip()
        station = asc_list[1].strip()
        # ev_st = f'{event}_{station}'
        eq_dir = asc_list[-1].strip()
        
        infile = path_join(folder, asc)
        outfile = path_join(folder, name + '.csv')

        if not exists(infile):
            raise IOError('File is not present')
        
        # print(infile) ; print(outfile) ; print(exists(infile))
        
        if name in EQ_dict.keys():
            print(f'File {asc} is already in the dictionary')
        elif not exists(infile):
            print(f'File {asc} does not exist')
        else:
            asc_dict = read_ASC(infile)
            units = asc_dict.get('UNITS', '')
            data_type= asc_dict.get('DATA_TYPE', '')            
            th = asc_dict.get('Data', [])
            asc_dict['Event'] = event
            asc_dict['Station'] = station
            asc_dict['Eq_dir'] = eq_dir
            EQ_dict[name] = asc_dict.copy()
            if debug:
                print(f'{data_type:40s}: File {asc}')
            
            # Export data to CSV files
            if export:
                with open(outfile, 'w') as csv_file:
                    csv_file.write(f'{name}\n')
                    if '_Acc_' in name:
                        csv_file.write(f'Time (s), Acceleration ({units})\n')
                    elif '_Vel_' in name:
                        csv_file.write(f'Time (s), Velocity ({units})\n')
                    elif '_Disp_' in name:
                        csv_file.write(f'Time (s), Displacement ({units})\n')
                    elif '_ResSpecAcc_' in name:
                        csv_file.write(f'Period (s), Acceleration ({units})\n')
                    elif '_ResSpecVel_' in name:
                        csv_file.write(f'Period (s), Velocity ({units})\n')
                    elif '_ResSpecDisp_' in name:
                        csv_file.write(f'Period (s), Displacement ({units})\n')
                    [csv_file.write(f'{x},{y}\n') for x, y in th]    
    
    return EQ_dict


def ASC_import_ss(asc_files, folder=r'./', export=False):
    """
    Imports a list of ASC files from a single folder and returns 
    a dictionary containing one entry for each file. In other words, 
    it does not identify pairs or triplets of files. 

    The dictionary can easily be read into a Pandas dataframe for
    further processing.

    This data is best suited for plotting with the 
    """
    EQ_dict = {}
    for asc in asc_files:
        name, _ = splitext(asc)
        asc_list = name.strip().split('_')
        event = asc_list[0].strip()
        station = asc_list[1].strip()
        # ev_st = f'{event}_{station}'
        eq_dir = asc_list[-1].strip()

        asc_dict = {'Input File': asc, 'EQ_DIR': eq_dir}
        
        infile = path_join(folder, asc)
        outfile = path_join(folder, name + '.csv')
        # print(infile) ; print(outfile) ; print(exists(infile))
        
        th = []
        t = 0
        record_type = ''
        
        #
        with open(infile, 'r') as asc_file:
            for line in asc_file:
                if len(line.strip()) == 0:
                    pass
                elif ':' in line:
                    key, value = line.split(':', 1)
                    asc_dict[key] = value.strip()
                elif ' ' in line.strip():
                    record_type = 'RS'
                    T, S = line.strip().split()
                    th.append((T, S))
                else:
                    samp_int = asc_dict.get('SAMPLING_INTERVAL_S', '').strip()
                    dt = float(samp_int)
                    if len(samp_int) > 0:
                        record_type = 'Acc'
                        th.append((t, line.strip()))
                        t = round(t + dt, len(samp_int))
        data = [(float(x), float(y)) for x, y in th]
        units = asc_dict.get('UNITS', '')
        asc_dict['RecType'] = record_type
        data_type = asc_dict.get('DATA_TYPE', '').strip()
            
        if record_type == 'Acc':            
            asc_dict['xlabel'] = 'Time (s)'
            asc_dict['ylabel'] = f'Acceleration ({units})'
        elif record_type == 'RS':
            data_type = asc_dict.get('DATA_TYPE', '').strip()
            if data_type == 'Acceleration Response Spectrum':
                asc_dict['xlabel'] = 'Period (s)'
                asc_dict['ylabel'] = f'Acceleration ({units})'
            elif data_type == 'Velocity Response Spectrum':
                asc_dict['xlabel'] = 'Period (s)'
                asc_dict['ylabel'] = f'Velocity ({units})'
            elif data_type == 'Displacement Response Spectrum':
                asc_dict['xlabel'] = 'Period (s)'
                asc_dict['ylabel'] = f'Displacement ({units})'
            else:
                asc_dict['xlabel'] = 'Period (s)'
                asc_dict['ylabel'] = f'Unknown ({units})'
        asc_dict['data'] = data
        EQ_dict[(event, station, record_type, eq_dir)] = asc_dict.copy()

        # Export data to CSV files
        if export:
            with open(outfile, 'w') as csv_file:
                
                csv_file.write(f'{name}\n')
                if '_Acc_' in name:
                    csv_file.write(f'Time (s), Acceleration ({units})\n')
                elif '_Vel_' in name:
                    csv_file.write(f'Time (s), Velocity ({units})\n')
                elif '_ResSpecAcc_' in name:
                    csv_file.write(f'Period (s), Acceleration ({units})\n')
                elif '_ResSpecVel_' in name:
                    csv_file.write(f'Period (s), Velocity ({units})\n')
                elif '_ResSpecDisp_' in name:
                    csv_file.write(f'Period (s), Displacement ({units})\n')
                [csv_file.write(f'{x},{y}\n') for x, y in th]
    return EQ_dict



def ASC_import2_ss(asc_files, folder=r'./'):
    """
    Returns a dictionary of two dictionaries - one for ground motions (key 'TH') and 
    one for response spectra (key 'RS'). The keys of each sub-dictionary are the truncated
    filenames. 

    This function ignores all metadata, apart from extracting 'SAMPLING_INTERVAL_S' 
    and 'UNITS' to get the delta-t for the acceleration records ('TH') and to 
    report the recorded measurement units respectively. 

    Note that it is hardwired for the filename formats in the resources folder:
        20230206011732_4628_ap_ResSpecAcc_E.asc -> key = "20230206011732_4628"

    This data is best suited for direct plotting using the EQ_Plots function.
    """
    TH_dict = {}
    RS_dict = {}

    for asc in asc_files:
        # print(splitext(asc))
        asc_list = asc.strip().replace('.asc','').split('_')
        event = asc_list[0].strip()
        station = asc_list[1].strip()
        ev_st = f'{event}_{station}'
        direction = station = asc_list[-1].strip()

        infile = path_join(folder, asc)
        # outfile = path_join(folder, asc.replace('.asc', '.csv'))
        # print(infile)
        # print(outfile)
        # print(exists(infile))
        
        if any(txt in asc for txt in ('_Disp_', '_Vel_', '_Acc_')):
            # print('Acc', asc)
            t = 0
            th = []
            asc_dict = {}
            th_dict = TH_dict.get(ev_st,{})
            
            with open(infile, 'r') as asc_file:
                for line in asc_file:
                    if len(line.strip()) == 0:
                        pass
                    elif ':' in line:
                        key, value = line.split(':', 1)
                        asc_dict[key] = value.strip()
                    elif ' ' in line.strip():
                        pass
                    else:
                        samp_int = asc_dict.get('SAMPLING_INTERVAL_S', '').strip()
                        dt = float(samp_int)
                        if len(samp_int) > 0:
                            th.append((t, line.strip()))
                            t = round(t + dt, len(samp_int))
            
            units = asc_dict.get('UNITS', '')
            data_type = asc_dict.get('DATA_TYPE', '').strip()
            if any(txt in data_type for txt in ('Displacement', 'Velocity', 'Acceleration')):
                dt = [txt in data_type for txt in ('Displacement', 'Velocity', 'Acceleration')][0]
            else:
                dt = ''
            th_dict['xlabel'] = 'Time (s)'
            th_dict['ylabel'] = f'{dt} ({units})'
            data = [(float(x), float(y)) for x, y in th]
            th_dict[direction] = data
            TH_dict[ev_st] = th_dict.copy()
            asc_dict = {}  # asc_dict is discarded

        elif any(txt in asc for txt in ('_ResSpecDisp_', '_ResSpecVel_', '_ResSpecAcc_')):
        
            # print('ResSpecAcc', asc)
            th = []
            asc_dict = {}
            rs_dict = RS_dict.get(ev_st,{})
            th = []
            
            with open(infile, 'r') as asc_file:
                for line in asc_file:
                    if len(line.strip()) == 0:
                        pass
                    elif ':' in line:
                        key, value = line.split(':', 1)
                        asc_dict[key] = value.strip()
                    elif ' ' in line.strip():
                        T, S = line.strip().split()
                        th.append((T, S))
                    
            units = asc_dict.get('UNITS', '')
            data_type = asc_dict.get('DATA_TYPE', '').strip()
            if any(txt in data_type for txt in ('Displacement', 'Velocity', 'Acceleration')):
                dt = [txt in data_type for txt in ('Displacement', 'Velocity', 'Acceleration')][0]
            else:
                dt = ''
            rs_dict['xlabel'] = 'Period (s)'
            
            rs_dict['ylabel'] = f'Spectral {dt} ({units})'
            data = [(float(x), float(y)) for x, y in th]
            rs_dict[direction] = data
            RS_dict[ev_st] = rs_dict.copy()
    
    return {'TH':TH_dict, 'RS':RS_dict}


def main():
    pass


if __name__ == "__main__":
    doctest.testmod()
    print('== Doctest Complete ==')
    main()

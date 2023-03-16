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

from os.path import join, splitext
# from os.path import exists, isdir, isfile, basename

def ASC_import(asc_files, folder=r'./', export=False):
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
        
        infile = join(folder, asc)
        outfile = join(folder, name + '.csv')
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

        if record_type == 'Acc':            
            asc_dict['xlabel'] = 'Time (s)'
            asc_dict['ylabel'] = f'Acceleration ({units})'
        elif record_type == 'RS':
            asc_dict['xlabel'] = 'Period (s)'
            asc_dict['ylabel'] = f'Acceleration ({units})'
        asc_dict['data'] = data
        EQ_dict[(event, station, record_type, eq_dir)] = asc_dict.copy()

        # Export data to CSV files
        if export:
            with open(outfile, 'w') as csv_file:
                
                csv_file.write(f'{name}\n')
                if '_Acc_' in name:
                    csv_file.write(f'Time (s), Acceleration ({units})\n')
                elif '_ResSpecAcc_' in name:
                    csv_file.write(f'Period (s), Acceleration ({units})\n')
                [csv_file.write(f'{x},{y}\n') for x, y in th]
    return EQ_dict


def ASC_import2(asc_files, folder=r'./'):
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

        infile = join(folder, asc)
        # outfile = join(folder, asc.replace('.asc', '.csv'))
        # print(infile)
        # print(outfile)
        # print(exists(infile))
        
        if '_Acc_' in asc:
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
            th_dict['xlabel'] = 'Time (s)'
            th_dict['ylabel'] = f'Acceleration ({units})'
            data = [(float(x), float(y)) for x, y in th]
            th_dict[direction] = data
            TH_dict[ev_st] = th_dict.copy()
            asc_dict = {}  # asc_dict is discarded

        elif '_ResSpecAcc_' in asc:
            print('ResSpecAcc', asc)
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
            rs_dict['xlabel'] = 'Period (s)'
            rs_dict['ylabel'] = f'Spectral Acceleration ({units})'
            data = [(float(x), float(y)) for x, y in th]
            rs_dict[direction] = data
            RS_dict[ev_st] = rs_dict.copy()
    
    return {'TH':TH_dict, 'RS':RS_dict}


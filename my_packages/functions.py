
'''
Basic Functions for MSanalyst
'''

import os
import ast
import json
import pandas as pd
import numpy as np
import spectrum_utils.spectrum as sus
from collections import namedtuple

def arrary2list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

TopK=namedtuple('topk',['index','number'])
def ex_startswith(file, start_txt):
    '''
    Extract lines starting with a specific keyword as **float** or **str**.
    :param filename: Path including suffix of the **text** file you intend to slice
    :param start_txt: Starting keyword
    :return: A list containing content after keywords
    '''
    with open(file, 'r') as f:
        content = [line[len(start_txt):].rstrip() for line in f if line.startswith(start_txt)]
    return content

def create_result_folders(args):
    '''
    Create result folders based on the quant file
    '''
    df = pd.read_csv(args.quant_file)
    parent_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result'# output/_quant_result/**
    os.makedirs(parent_folder, exist_ok=True)
    for _, row in df.iterrows():
        folder_name = f"{parent_folder}/{int(row['row ID'])}"
        os.makedirs(folder_name, exist_ok=True)
    print('Result folders have been created!')

def create_subresults(args):
    '''
    Split the results  after the MS1 match
    Create a separate CSV for each row ID, writing the corresponding information to facilitate detailed inspection
    '''
    parent_folder =  f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result' # filename ''output/_quant_result/**''
    npms1_result_path =os.path.join(parent_folder, f'IS_MS1match_{os.path.basename(args.quant_file)}')
    edbms1_result_path = os.path.join(parent_folder, f'E_MS1match_{os.path.basename(args.quant_file)}')

    quant_df = df_preprocess(args.quant_file)
    npms1_match_df = df_preprocess(npms1_result_path)
    edbms1_match_df = df_preprocess(edbms1_result_path)


    for i in range(len(quant_df)):
        id = quant_df['row ID'][i]
        folder_name = os.path.join(parent_folder, str(id))

        npcsv_file = os.path.join(folder_name, f'IS_MS1match_{str(id)}.csv') # isdb results
        if not os.path.exists(npcsv_file):
            pd.DataFrame(columns=npms1_match_df.columns).to_csv(npcsv_file, index=False)
        selected_rows =npms1_match_df.loc[npms1_match_df['row ID'] == id]
        with open(npcsv_file, 'a', newline='') as f1:
            selected_rows.to_csv(f1, index=False, header=False)

        edbcsv_file = os.path.join(folder_name, f'E_MS1match_{str(id)}.csv') # edb result
        if not os.path.exists(edbcsv_file):
            pd.DataFrame(columns=edbms1_match_df.columns).to_csv(edbcsv_file, index=False)
        selected_rows = edbms1_match_df.loc[edbms1_match_df['row ID'] == id]
        with open(edbcsv_file, 'a', newline='') as f2:
            selected_rows.to_csv(f2, index=False, header=False)

def ex_spectra(file, start_txt, end_txt, skip_words=None):
    '''
    Horizontal and vertical coordinates of tandem mass
    :param file:
    :param start_txt:
    :param end_txt:
    :param skip_words:
    :return: A list contain lists of sliced content, like[[],[],...,[]],and converting to an array
    '''
    if skip_words == None:
        skip_words = []
    spectra = []
    with open(file, 'r') as f:
        lines = f.readlines()
        start_idx = 0
        for i in range(len(lines)):
            if start_txt in lines[i]:
                if any(word in lines[i + 1] for word in skip_words):
                    start_idx = i+2
                else:
                    start_idx = i+1
            elif end_txt in lines[i]:
                spectrum = ''.join(lines[start_idx:i])
                spectra_list = spectrum.split('\n')[:-1]
                temp=[]
                for s in spectra_list:
                    m_z, intensity = s.split()
                    temp.append([float(m_z), float(intensity)])
                temp = np.array(temp,dtype=np.float64)
                spectra.append(temp)
    return spectra

def mgf_process(mgf_file):
    '''
    Process MGF file to extract relevant information.
    :param mgf_file: '.mgf'
    :return: id<str> pepmass<str>, ms2<np array>
    '''
    id_txt = 'FEATURE_ID='
    id = ex_startswith(mgf_file, id_txt)

    pepmass_txt = 'PEPMASS='
    pepmass = ex_startswith(mgf_file, pepmass_txt)

    charge_txt = 'CHARGE='
    charge = ex_startswith(mgf_file, charge_txt)
    charge = [s.replace('+', '') for s in charge]

    start_txt = 'MSLEVEL=2'
    end_txt = 'END'
    ms2 = ex_spectra(mgf_file, start_txt, end_txt, skip_words=['MERGED'])

    exp_info = pd.DataFrame({
        'id': id
        ,'pepmass': pepmass
        ,'charge' : charge
        ,'ms2': ms2
    })
    exp_info = exp_info[exp_info['ms2'].apply(len) > 1]  # delete empty list
    exp_info = exp_info.reset_index(drop=True)  # reindex
    return exp_info

def spectra_process(qms1,qms2):
    '''
    Directly process input query MS1 and MS2 spectra
    :param qms1: e.g. ''
    :param qms2: e.g. '381.2958'
    :return: e.g. '381.2284 1.0E2 381.2344 1.1E2 381.2822 1.1E2 381.2842 1.3E2 381.2862 5.2E2'
    '''
    id='1'
    pepmass= qms1
    charge='1'
    try:
        spectra = []
        temp = []
        lines = qms2.strip().split('\n')
        for line in lines:
            m_z, intensity = line.split()
            temp.append([float(m_z), float(intensity)])
        temp = np.array(temp, dtype=np.float64)
        spectra.append(temp)
    except:
        spectra = []
        temp = []
        elements = qms2.split()
        for i in range(0,len(elements),2):
            temp.append([float(elements[i]),float(elements[i+1])])
        temp = np.array(temp,dtype=np.float64)
        spectra.append(temp)

    exp_info = pd.DataFrame({
        'id': id
        , 'pepmass': pepmass
        , 'charge': charge
        , 'ms2': spectra
    })
    exp_info = exp_info[exp_info['ms2'].apply(len) > 1]  # delete empty list
    exp_info = exp_info.reset_index(drop=True)
    return exp_info

def get_mgf_info(mgf_info,mgf_id):
    '''
    Retrieve information from MGF file based on ID.
    :param mgf_info:
    :param id:
    :return:pepmass<float>, spec<np.adarray>, spectrum<spectrum_utils object>
    '''
    if not mgf_info.empty:
        pepmass = float(mgf_info[mgf_info['id'] == mgf_id]['pepmass'].iloc[0])
        charge = int(mgf_info[mgf_info['id'] == mgf_id]['charge'].iloc[0])
        spec = mgf_info[mgf_info['id'] == mgf_id]['ms2'].iloc[0]
        mz = np.array(spec[:, 0])
        spectrum = sus.MsmsSpectrum(identifier=mgf_id
                                    , precursor_mz=pepmass
                                    , precursor_charge=charge
                                    , mz=mz
                                    , intensity=spec[:, 1])
        return {'pepmass': pepmass, 'spec': spec, 'spectrum': spectrum, 'charge': charge, 'id': mgf_id}
    else:
        raise ValueError(f"No data found for mgf_id: {mgf_id}")

def get_edb_info(gnps_info, gnps_id):
    '''

    :param isdb_info:
    :param id:
    :return:
    '''
    keys_to_retrieve = ['smiles', 'pepmass', 'ms2','charge']
    values = [gnps_info[gnps_id][key] for key in keys_to_retrieve]
    smiles, pepmass, spec, charge = values
    # string convertion
    pepmass = float(pepmass)
    charge = int(charge)
    spec = np.asarray(ast.literal_eval(spec))
    mz = np.array(spec[:, 0])
    spectrum = sus.MsmsSpectrum(identifier=f'{gnps_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=charge
                                 , mz=mz
                                 , intensity=spec[:, 1])

    return {'smiles': smiles, 'pepmass': pepmass
        , 'spec': spec, 'spectrum': spectrum,'charge': charge}

def get_isdb_info(isdb_info, is_id):
    '''

    :param isdb_info:
    :param id:
    :return:
    '''
    keys_to_retrieve = ['smiles', 'pepmass', 'energy0_ms2', 'energy1_ms2', 'energy2_ms2']
    values = [isdb_info[is_id][key] for key in keys_to_retrieve]
    smiles, pepmass, e0spec, e1spec, e2spec = values
    # string convertion
    pepmass = float(pepmass)
    e0spec = np.asarray(ast.literal_eval(e0spec))
    e1spec = np.asarray(ast.literal_eval(e1spec))
    e2spec = np.asarray(ast.literal_eval(e2spec))

    mz0 = np.array(e0spec[:, 0])
    spectrum0 = sus.MsmsSpectrum(identifier=f'e0_{is_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=1
                                 , mz=mz0
                                 , intensity=e0spec[:, 1])
    mz1 = np.array(e1spec[:, 0])
    spectrum1 = sus.MsmsSpectrum(identifier = f'e1_{is_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=1
                                 , mz=mz1
                                 , intensity=e1spec[:, 1])
    mz2 = np.array(e2spec[:, 0])
    spectrum2 = sus.MsmsSpectrum(identifier=f'e2_{is_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=1
                                 , mz=mz2
                                 , intensity=e2spec[:, 1])

    return {'smiles': smiles, 'pepmass': pepmass
        , 'e0spec': e0spec, 'e1spec': e1spec, 'e2spec': e2spec
        , 'e0spectrum': spectrum0, 'e1spectrum': spectrum1, 'e2spectrum': spectrum2}

def df_preprocess(filename):
    '''
    Preprocess DataFrame by removing empty columns and resetting index.
    '''
    if filename.endswith('.csv'):
        df = pd.read_csv(filename, low_memory=False)
    elif filename.endswith('.tsv'):
        df = pd.read_csv(filename, sep='\t', low_memory=False)
    elif filename.endswith('.xlsx') or filename.endswith('.xls'):
        df = pd.read_excel(filename)
    else:
        raise ValueError("Unsupported file format. Please use .csv, .tsv, or .xlsx files.")

    if  df.index[-1] != len(df)-1:
        df.index.name = ''
        df.reset_index(inplace=True)
    return df

def calculate_ppm(query_mass_value: float, reference_mass_value: float) -> float:
    '''
    Calculate parts per million (ppm) for mass values.
    '''
    if not isinstance(query_mass_value, (int, float)) or not isinstance(reference_mass_value, (int, float)):
        raise TypeError('Input parameters must be numbers.')
    if reference_mass_value != 0:
        return abs((query_mass_value - reference_mass_value) / reference_mass_value * 1e6)
    return float('inf')

def db_parsing():
    '''
    Parse default databases of MSanalyst.
    '''
    isdb_file = './msdb/isdb_info.json'
    edb_file = './msdb/edb_info.json'
    with open(isdb_file, 'r') as f:
        isdb_info = json.load(f)
    with open(edb_file, 'r') as f1:
        gnps_info = json.load(f1)
    return isdb_info, gnps_info



if __name__ == '__main__':
    print('')
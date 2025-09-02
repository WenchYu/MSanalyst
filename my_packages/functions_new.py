
'''
Basic Functions for MSanalyst
'''

import os,re,ast,json,ujson,glob,spectral_entropy
import spectral_entropy as se
import pandas as pd
import numpy as np
import spectrum_utils.spectrum as sus
import matchms.filtering as ms_filters
import pyopenms as oms
from matchms import Spectrum
from matchms.importing import load_from_json, load_from_mgf,load_from_msp, load_from_mzml, load_from_mzxml, load_from_usi
from matchms.similarity import ModifiedCosine,NeutralLossesCosine,CosineGreedy
from spectral_entropy import similarity
from ms_entropy.file_io import spec_file
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from collections import namedtuple
from my_packages import peaktools,config,ms2tools_new
from tqdm import trange,tqdm

TopK = namedtuple('topk',['index','number']) # e.g. topk(index=2, number=8)

class MgfProcess:
    def __init__(self, MGF_FILE = None):
        self.MGF_FILE = MGF_FILE
        if 'GNPS' not in MGF_FILE:
            self.MGF_INFO = self.mgf_extract()
        else:
            self.MGF_INFO = self.gnps_mgf_extract()
        # if MGF_FILE is not None:
        #     self.MGF_FILE = MGF_FILE
        #     # self.MGF_FILE = self.mgf_process(MGF_FILE)
        # else:
        #     self.MGF_FILE = pd.DataFrame()

    def line_extract(self, START_TEXT):
        '''
        Extract all the lines starting with a specific keyword and save as a idlist
        :param MGF_FILE: Path including suffix of the **text** file you intend to slice
        :param START_TEXT: Starting keyword
        :return: A idlist containing content <str> after keywords
        '''
        with open(self.MGF_FILE, 'r') as f:
            CONTENT = [line[len(START_TEXT):].rstrip() for line in f if line.startswith(START_TEXT)]
        return CONTENT

    def spectra_extract(self, START_TEXT, END_TEXT, skip_words=None):
        '''
        Extracting Horizontal and vertical coordinates of tandem mass
        :param file:
        :param start_txt:
        :param end_txt:
        :param skip_words:
        :return: A idlist contain lists of sliced content, like[[],[],...,[]],and converting to an array
        '''
        if skip_words is None:
            skip_words = []
        spectra = []
        with open(self.MGF_FILE, 'r') as f:
            lines = f.readlines()
            start_idx = 0
            for i in range(len(lines)):
                if START_TEXT in lines[i]:
                    if any(word in lines[i + 1] for word in skip_words):
                        start_idx = i + 2
                    else:
                        start_idx = i + 1
                elif END_TEXT in lines[i]:
                    spectrum = ''.join(lines[start_idx:i])
                    spectra_list = spectrum.split('\n')[:-1]
                    temp = []
                    for s in spectra_list:
                        m_z, intensity = s.split()
                        temp.append([float(m_z), float(intensity)])
                    # temp = np.array(temp, dtype=np.float64) # this led to the introduction of ... and incomplete spectra
                    spectra.append(temp)
        return spectra

    def mgf_extract(self):
        '''
        Process MGF file to extract relevant information.
        :param mgf_file: '.mgf'
        :return: id<str> pepmass<str>, ms2<np array>
        '''
        ID_TEXT = 'FEATURE_ID='
        PEPMASS_TEXT = 'PEPMASS='
        CHARGE_TEXT = 'CHARGE='

        id = self.line_extract(ID_TEXT)
        pepmass = self.line_extract(PEPMASS_TEXT)
        charge = self.line_extract(CHARGE_TEXT)
        charge = [s.replace('+', '') for s in charge]

        start_txt = 'MSLEVEL=2'
        end_txt = 'END'
        ms2 = self.spectra_extract(start_txt, end_txt, skip_words=['MERGED'])

        exp_info = pd.DataFrame({
            'id': id,
            'pepmass': pepmass,
            'charge': charge,
            'ms2': ms2
        })
        exp_info = exp_info[exp_info['ms2'].apply(len) > 1]  # delete empty idlist
        exp_info = exp_info.reset_index(drop=True)  # reindex
        return exp_info

    def feature_extract(self, FEATURE_ID):
        '''
        Retrieve information from MGF file based on ID.
        :param mgf_id:
        :return: pepmass<float>, spec<np.array>, spectrum<vsl object>
        '''
        mgf_info = self.MGF_INFO
        try:
            pepmass = float(mgf_info[mgf_info['id'] == FEATURE_ID]['pepmass'].iloc[0])
            charge = int(mgf_info[mgf_info['id'] == FEATURE_ID]['charge'].iloc[0])
            spec = mgf_info[mgf_info['id'] == FEATURE_ID]['ms2'].iloc[0]
            spec = se.clean_spectrum(spec, max_mz=pepmass + 0.01)

            mz = np.array(spec[:, 0])
            intensity = np.array(spec[:, 1])
            spectrum = sus.MsmsSpectrum(
                identifier=FEATURE_ID,
                precursor_mz=pepmass,
                precursor_charge=charge,
                mz=mz,
                intensity=intensity
            )
            return {
                'pepmass': pepmass,
                'spec': spec,
                'spectrum': spectrum,
                'charge': charge,
                'id': FEATURE_ID
            }
        except:
            raise ValueError(f"No data found for mgf_id: {FEATURE_ID}")

    def gnps_mgf_extract(self):
        '''
        Process MGF file to extract relevant information.
        :param mgf_file: '.mgf'
        :return: a dict used to generate json file
        id<str> pepmass<str>, ms2<np array>
        '''

        SPECTRUMID = self.line_extract('SPECTRUMID=')
        PEPMASS = self.line_extract('PEPMASS=')
        CHARGE = self.line_extract('CHARGE=')
        MSLEVEL = self.line_extract('MSLEVEL=')
        IONMODE = self.line_extract('IONMODE=')
        NAME = self.line_extract('NAME=')
        SMILE = self.line_extract('SMILES=')
        INSTRUMENT = self.line_extract('SOURCE_INSTRUMENT=')
        # charge = [s.replace('+', '') for s in charge]

        START_TEXT = 'SCAN'
        END_TEXT = 'END'
        MS2_SPECTRUM = self.spectra_extract(START_TEXT, END_TEXT, skip_words=['MERGED'])

        dict = {}
        print('Start to convert to dict')
        for i in trange(len(SPECTRUMID)):
            dict[f'{SPECTRUMID[i]}'] = {
                'PEPMASS': f'{PEPMASS[i]}',
                'CHARGE': f'{CHARGE[i]}',
                'MSLEVEL': f'{MSLEVEL[i]}',
                'IONMODE': f'{IONMODE[i]}',
                'NAME': f'{NAME[i]}',
                'SMILE': f'{SMILE[i]}',
                'INSTRUMENT': f'{INSTRUMENT[i]}',
                'MS2_SPECTRUM': f'{MS2_SPECTRUM[i]}'
                }

        return dict

    def query_mass_process(self, QMS1, QMS2):
        '''
        Specifically designed for MSanalyst
        Directly process input query MS1 and MS2 spectra
        :param qms1: e.g. '381.2958'
        :param qms2: e.g. '381.2284 1.0E2 381.2344 1.1E2 381.2822 1.1E2 381.2842 1.3E2 381.2862 5.2E2'
        :return: e.g. '381.2284 1.0E2 381.2344 1.1E2 381.2822 1.1E2 381.2842 1.3E2 381.2862 5.2E2'
        '''
        id = '1'
        pepmass = QMS1
        charge = '1'
        try:
            spectra = []
            temp = []
            lines = QMS2.strip().split('\n')
            for line in lines:
                m_z, intensity = line.split()
                temp.append([float(m_z), float(intensity)])
            temp = np.array(temp, dtype=np.float64)
            spectra.append(temp)
        except:
            spectra = []
            temp = []
            elements = QMS2.split()
            for i in range(0, len(elements), 2):
                temp.append([float(elements[i]), float(elements[i + 1])])
            temp = np.array(temp, dtype=np.float64)
            spectra.append(temp)

        exp_info = pd.DataFrame({
            'id': [id],
            'pepmass': [pepmass],
            'charge': [charge],
            'ms2': [spectra[0]]  # spectra is a idlist containing one numpy array
        })
        return exp_info

def load_spectra_from_file(INPUT_FILE):
    '''
    It will return a list containing spectral objects
    Properties can be accessed by (.mz, .intensities, .metadata)
    '''
    try:
        FILE_TYPE = spec_file.guess_file_type_from_file_name(INPUT_FILE) # msp, mgf, mzml, mzml, hdf5, raw, lbm2
    except:
        FILE_TYPE = INPUT_FILE


    if FILE_TYPE == 'json':
        SPECTRA = list(load_from_json(INPUT_FILE))
    if FILE_TYPE == 'mgf':
        SPECTRA = list(load_from_mgf(INPUT_FILE))
    if FILE_TYPE == 'msp':
        SPECTRA = list(load_from_msp(INPUT_FILE))
    if FILE_TYPE == 'mzml':
        SPECTRA = list(load_from_mzml(INPUT_FILE))
    if FILE_TYPE == 'mzxml':
        SPECTRA = list(load_from_mzxml(INPUT_FILE))
    if 'mzspec:' in INPUT_FILE:
         SPECTRA = list(load_from_usi(INPUT_FILE))

    return SPECTRA

def spec_str2array(SPEC_STR,PEPMASS):
    '''
    Convert MS2_spectrum_str to MS2_spectrum_array
    :param SPEC_STR:
    :param PEPMASS:
    :return:
    '''
    SPEC_STR = SPEC_STR.replace(',', ' ').replace('[', ' ').replace(']', ' ')  # Formatting the spec_str[mz ins mz ins mz ins ...]
    SPEC_STR = SPEC_STR.split()  # Split by space
    SPEC_FLOAT = [float(x) for x in SPEC_STR]  # floating
    SPEC_ARRAY = np.array(SPEC_FLOAT).reshape(-1, 2)  # mz ins\n mz ins\n ... <np.array>
    SPEC_ARRAY = se.clean_spectrum(SPEC_ARRAY, max_mz=PEPMASS + 0.01)
    return SPEC_ARRAY

def gnps_info_format(GNPS_INFO,CCMSID):
    '''
    Json only support text format, conversion to np.array is needed.
    :param GNPS_INFO:
    :param CCMSID:
    :return:
    '''
    EX_INFO = GNPS_INFO[CCMSID]
    spec_str = EX_INFO['MS2_SPECTRUM']
    PEPMASS = float(EX_INFO['PEPMASS'])

    spec_str = spec_str.replace(',', '').replace('[', '').replace(']', '')  # Formatting the spec_str[mz ins mz ins mz ins ...]
    spec_str = spec_str.split()  # Split by space
    spec_float = [float(x) for x in spec_str]  # floating
    spec_array = np.array(spec_float).reshape(-1, 2)  # mz ins\n mz ins\n ... <np.array>
    spec_array = se.clean_spectrum(spec_array, max_mz=PEPMASS + 0.01)

    SPECTRUM = Spectrum(mz=np.array(spec_array[:, 0], dtype=float),
                          intensities=np.array(spec_array[:, 1], dtype=float),
                          metadata={"precursor_mz": PEPMASS})
    return spec_array,SPECTRUM

def json_load(JSONFILE):
    with open(JSONFILE,'r') as f:
        return ujson.load(f)

def json_dump(JSONFILE,TODUMP):
    with open(JSONFILE,'r') as f:
        return ujson.dump(TODUMP,f)

def refms_compare(GNPS_INFO, ALGORITHM, CCMSID1, CCMSID2):

    SPEC1, SPECTRUM1 = gnps_info_format(GNPS_INFO, CCMSID1)
    SPEC2, SPECTRUM2 = gnps_info_format(GNPS_INFO, CCMSID2)

    if ALGORITHM == 'cosine':
        return peaktools.cosine(SPECTRUM1, SPECTRUM2, 0.05)
    elif ALGORITHM == 'modified_cosine':
        return peaktools.modified_cosine(SPECTRUM1, SPECTRUM2, 0.05)
    elif ALGORITHM == 'neutral_loss':
        return peaktools.neutral_loss(SPECTRUM1, SPECTRUM2, 0.05)
    else:
        return se.similarity(SPEC1, SPEC2, method=ALGORITHM, ms2_da=0.05)

def arrary2list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

def create_result_folders(args):
    '''
    Create result folders based on the input file, subfolders based on the feature ID
    args.output
    args.spectra_file
    :param args:
    :return:
    '''
    parent_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.spectra_file))[0]}_result'# output/_quant_result/**
    os.makedirs(parent_folder, exist_ok=True)
    SPECTRA = load_spectra_from_file(args.spectra_file)
    print(len(SPECTRA))
    for SPECTRUM in SPECTRA:  # STD data
        FID = SPECTRUM.metadata['feature_id']
        folder_name = f"{parent_folder}/{FID}"
        os.makedirs(folder_name, exist_ok=True)
    print('Result folders have been created!')

def create_subresults(args):
    '''
    Split the results  after the MS1 match
    Create a separate CSV for each row ID, writing the corresponding information to facilitate detailed inspection
    '''

    basename = os.path.splitext(os.path.basename(args.spectra_file))[0]
    parent_folder = f'{args.output}/{basename}_result'  # filename ''output/_result/**''

    npms1_result_path =os.path.join(parent_folder, f'IS_MS1match_{basename}.csv')
    edbms1_result_path = os.path.join(parent_folder, f'E_MS1match_{basename}.csv')
    npms1_match_df = df_preprocess(npms1_result_path)
    edbms1_match_df = df_preprocess(edbms1_result_path)

    result_dir = os.path.join(args.output, f'{basename}_result')
    edb_result_path = os.path.join(result_dir, f'E_MS1match_{basename}.csv')
    quant_df = df_preprocess(edb_result_path)



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

def list_files(directory,keyword):
    '''
    idlist files with keyword
    :param directory: dirt
    :return: A idlist containing all .graphml file paths
    '''
    graphml_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(keyword):
                graphml_files.append(os.path.join(root, file))
    return graphml_files

def ex_algorithm_name(file,prefix):
    '''
    extract algorithm name from '.graphml' files
    e.g. dot_product from KutzOsmac_dot_product_0.7_3.graphml
    :return:
    '''
    match = re.search(r'^(\w+)_|(\w+)_\d+\.\d+', file)
    if match:
        pattern = match.group(1) if match.group(1) else match.group(2)
        pattern = pattern.replace(f'{prefix}_', '')
        return pattern

def spectra_to_fsspectra(SPECTRA):
    '''

    :param SPECTRA: output of load_spectra_from_file
    :return:
    '''
    FS_SPECTRA = []
    for SPECTRUM in SPECTRA:  # STD data
        FID = SPECTRUM.metadata['feature_id']
        PM = SPECTRUM.metadata['precursor_mz']
        PEAKs = np.column_stack((SPECTRUM.mz, SPECTRUM.intensities))
        FS_SPECTRA.append(
            {
                "id": FID,
                "precursor_mz": PM,
                "peaks": PEAKs
            })

    return FS_SPECTRA

def get_compare_dict(FS_SPECTRA,LIBRARY):
    '''

    param:FS SPECTRA
    param:FS library
    all_ions_mz: all the m/z in flash lib
    all_ions_mz_idx_start: all the start idx of spectrum
    all_ions_spec_idx: all the start idx of spectra
    library_intensity: all the intensities in flash lib

    '''

    # all_ions_mz,all_ions_mz_idx_start,all_ions_spec_idx = e_ions_mz,e_ions_mz_idx_start,e_ions_spec_idx
    # library_intensity = e_ions_intensity

    # Param
    ms2_tolerance_in_da = 0.02
    mz_index_step = 0.0001
    index_number_in_one_da = int(1 / mz_index_step)

    ems2_search = FlashEntropySearchCore()
    (all_ions_mz_idx_start, all_ions_mz, all_ions_intensity, all_ions_spec_idx,
     all_nl_mass_idx_start, all_nl_mass, all_nl_intensity, all_nl_spec_idx, all_ions_idx_for_nl
     )= ems2_search.build_index(LIBRARY) # Load experimental library

    # Get Compare_dict
    temp_compare_dict = {}
    for idx, spectrum in enumerate(FS_SPECTRA):
        compare_ids = []
        for mz_query, intensity_query in spectrum['peaks']:
            # Determine the mz index range
            product_mz_idx_min = ems2_search._find_location_from_array_with_index(
                mz_query - ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "left", index_number_in_one_da)
            product_mz_idx_max = ems2_search._find_location_from_array_with_index(
                mz_query + ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "right", index_number_in_one_da)

            ref_ids = all_ions_spec_idx[product_mz_idx_min:product_mz_idx_max]  # <class 'numpy.ndarray'>
            compare_ids.extend(ref_ids)
        compare_ids = np.array(list(set(compare_ids)))
        temp_compare_dict[idx] = compare_ids

    # Remove replicate pairs in compare_dict(e.g. remove (8,0), becaues of former exist of (0,8))
    keys_to_remove = []
    compare_dict = {}  # {query_idx: [ref_idx1,ref_idx2,ref_idx2, ...]} {<int>:<list>}
    for key, values in temp_compare_dict.items():
        keys_to_remove.append(key)
        filtered_arr = values[~np.isin(values, keys_to_remove)]
        compare_dict[key] = filtered_arr

    return compare_dict

def spectra_process(SPECTRA_FILE):
    '''
    Process SPECTRA file to extract relevant information.
    :param
    :return: DICT{ID:{}}
    '''
    SPECTRA = load_spectra_from_file(SPECTRA_FILE) #.mz, .intensities, .metadata['feature_id'], .metadata['precursor_mz'],.metadata['charge']

    SPECTRA_INFO = {}
    for SPECTRUM in SPECTRA:
        ID = SPECTRUM.metadata['feature_id']
        PEAK = np.column_stack((SPECTRUM.mz, SPECTRUM.intensities))
        SPECTRA_INFO[ID]= {'pepmass':SPECTRUM.metadata['precursor_mz'],
                           'charge':SPECTRUM.metadata['charge'],
                           'ms2':PEAK
                           }

    return SPECTRA_INFO

def get_spectra_from_info(SPEC_INFO,ID):
    '''

    :param SPECTRA: output of load_from_spectra
    :param id:
    :return:
    '''
    TARGET_SPECTRUM = SPEC_INFO[ID]
    exp_pm = float(TARGET_SPECTRUM['pepmass']) # pepmass of query feature
    SPEC =TARGET_SPECTRUM['ms2'] # ms2 of query feature
    SPEC = spectral_entropy.clean_spectrum(SPEC, max_mz=exp_pm+0.01) # MS2 spectrum clean by normalizing and removing signals with intensity less than 1% of the base peak
    SPECTRUM = Spectrum(mz=np.array(SPEC[:, 0],dtype = float),
                  intensities=np.array(SPEC[:, 1],dtype = float),
                  metadata={"precursor_mz": exp_pm})
    return SPEC, SPECTRUM

def get_spectra_from_einfo(EINFO,ID):
    '''

    :param einfo:
    :param id:
    :return:
    '''
    edb_pm = float(EINFO[ID]['pepmass'])
    ESMILE = EINFO[ID]['smiles']
    temp = EINFO[ID]['ms2'] # <list>[[],[],...]
    if type(temp) == str:
        ESPEC = np.asarray(ast.literal_eval(temp)) # "[[],[],[],...]"
    else:
        ESPEC = temp
    ESPEC = spectral_entropy.clean_spectrum(ESPEC, max_mz=edb_pm+0.01)
    ESPECTRUM = Spectrum(mz=np.array(ESPEC[:, 0],dtype = float),
                  intensities=np.array(ESPEC[:, 1],dtype = float),
                  metadata={"precursor_mz": edb_pm,"smile":ESMILE})

    return ESPEC,ESPECTRUM

def get_spectra_from_isinfo(ISINFO,ID):
    '''

    :param ISINFO:
    :param ID:
    :return:
    '''
    is_pm = ISINFO[ID]['pepmass']
    is_smile = ISINFO[ID]['smiles']
    # 0, 1, 2 represent the in-silico MS2 predicted by CFM-ID at three collision energy (10, 20 and 40 eV)
    e0_ms2 = np.asarray(ast.literal_eval(ISINFO[ID]['energy0_ms2']))
    e0_spectrum = Spectrum(mz=np.array(e0_ms2[:, 0],dtype = float),
                  intensities=np.array(e0_ms2[:, 1],dtype = float),
                  metadata={"precursor_mz": is_pm+0.1,"smile":is_smile})

    e1_ms2 = np.asarray(ast.literal_eval(ISINFO[ID]['energy1_ms2']))
    e1_spectrum = Spectrum(mz=np.array(e1_ms2[:, 0],dtype = float),
                  intensities=np.array(e1_ms2[:, 1],dtype = float),
                  metadata={"precursor_mz": is_pm+0.1,"smile":is_smile})

    e2_ms2 = np.asarray(ast.literal_eval(ISINFO[ID]['energy2_ms2']))
    e2_spectrum = Spectrum(mz=np.array(e2_ms2[:, 0],dtype = float),
                  intensities=np.array(e2_ms2[:, 1],dtype = float),
                  metadata={"precursor_mz": is_pm+0.1,"smile":is_smile})

    return e0_ms2,e0_spectrum,e1_ms2,e1_spectrum,e2_ms2,e2_spectrum

def clac_spec_sim(SPEC1,SPECTRUM1,SPEC2,SPECTRUM2,ALGO):
    ''''''
    cosine = CosineGreedy(tolerance=0.05)
    neutral_loss = NeutralLossesCosine(tolerance=0.05)
    modified_cosine = ModifiedCosine(tolerance=0.05)

    if ALGO == 'modified_cosine':
        score = modified_cosine.pair(SPECTRUM1, SPECTRUM2)
        SPEC_SIM = score['score'].item()
        N_PEAK = score['matches'].item()

    elif ALGO == 'cosine':
        score = cosine.pair(SPECTRUM1, SPECTRUM2)
        SPEC_SIM = score['score'].item()
        N_PEAK = score['matches'].item()

    elif ALGO == 'neutral_loss':
        score = neutral_loss.pair(SPECTRUM1, SPECTRUM2)
        SPEC_SIM = score['score'].item()
        N_PEAK = score['matches'].item()

    elif ALGO == 'peak_percentage':
        score = cosine.pair(SPECTRUM1, SPECTRUM2)
        N_PEAK = score['matches'].item()
        SPEC_SIM = N_PEAK / min(len(SPEC1), len(SPEC2))

    else:
        score = cosine.pair(SPECTRUM1, SPECTRUM2)
        N_PEAK = score['matches'].item()
        SPEC_SIM = similarity(SPEC1, SPEC2, method=ALGO, ms2_da=0.05)

    return SPEC_SIM,N_PEAK

def preprocess_mzml(args):
    mzML_files = glob.glob(os.path.join(args.input_folder, '*.mzML'))
    feature_maps = []
    for file in mzML_files:
        # load mzML file into MSExperiment
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file, exp)  # load each mzML file to an OpenMS file format (MSExperiment)

        # mass trace detection
        mass_traces = ([])  # introduce an empty list where the mass traces will be loaded
        mtd = oms.MassTraceDetection()
        mtd_par = (mtd.getDefaults())  # get the default parameters in order to edit them
        mtd_par.setValue("mass_error_ppm", args.allowed_mass_tolerance)  # high-res instrument, orbitraps
        mtd_par.setValue("noise_threshold_int", args.noise_level)  # data-dependent (usually works for orbitraps)
        mtd.setParameters(mtd_par)  # set the new parameters
        mtd.run(exp, mass_traces, 0)  # run mass trace detection

        # elution peak detection
        mass_traces_deconvol = []
        epd = oms.ElutionPeakDetection()
        epd_par = epd.getDefaults()
        epd_par.setValue("width_filtering", "fixed")  # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
        epd.setParameters(epd_par)
        epd.detectPeaks(mass_traces, mass_traces_deconvol)

        # feature detection
        feature_map = oms.FeatureMap()  # output features
        chrom_out = []  # output chromatograms
        ffm = oms.FeatureFindingMetabo()
        ffm_par = ffm.getDefaults()
        ffm_par.setValue("remove_single_traces", "true")  # remove mass traces without satellite isotopic traces
        ffm.setParameters(ffm_par)
        ffm.run(mass_traces_deconvol, feature_map, chrom_out)
        feature_map.setUniqueIds()  # Assigns a new, valid unique id per feature
        feature_map.setPrimaryMSRunPath([file.encode()])  # Sets the file path to the primary MS run (usually the mzML file)
        feature_maps.append(feature_map)


    ''' 
    use as reference for alignment, the file with the largest number of features, (works well if you have a pooled QC for example)
    '''
    ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    trafos = {}
    # parameter optimization
    aligner_par = aligner.getDefaults()
    aligner_par.setValue("max_num_peaks_considered", -1)  # infinite
    aligner_par.setValue(
        "pairfinder:distance_MZ:max_difference", args.allowed_mass_tolerance
    )  # Never pair features with larger m/z distance
    aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
    aligner.setParameters(aligner_par)
    aligner.setReference(feature_maps[ref_index])

    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1:]:
        trafo = oms.TransformationDescription()  # save the transformed data points
        aligner.align(feature_map, trafo)
        trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
        transformer = oms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(feature_map, trafo, True)



    '''Align mzML files based on FeatureMap alignment (optional, only for GNPS).'''
    for file in mzML_files:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file, exp)
        exp.sortSpectra(True)
        exp.setMetaValue("mzML_path", file)
        if file not in trafos.keys():
            oms.MzMLFile().store(file[:-5] + "_aligned.mzML", exp)
            continue
        transformer = oms.MapAlignmentTransformer()
        trafo_description = trafos[file]
        transformer.transformRetentionTimes(exp, trafo_description, True)
        oms.MzMLFile().store(file[:-5] + "_aligned.mzML", exp)
    mzML_files = [file[:-5] + "_aligned.mzML" for file in mzML_files]

    '''Map MS2 to features'''
    feature_maps_mapped = []
    use_centroid_rt = False
    use_centroid_mz = True
    mapper = oms.IDMapper()
    for file in mzML_files:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file, exp)
        for i, feature_map in enumerate(feature_maps):
            if feature_map.getMetaValue("spectra_data")[
                0
            ].decode() == exp.getMetaValue("mzML_path"):
                peptide_ids = []
                protein_ids = []
                mapper.annotate(
                    feature_map,
                    peptide_ids,
                    protein_ids,
                    use_centroid_rt,
                    use_centroid_mz,
                    exp,
                )
                fm_new = oms.FeatureMap(feature_map)
                fm_new.clear(False)
                # set unique identifiers to protein and peptide identifications
                prot_ids = []
                if len(feature_map.getProteinIdentifications()) > 0:
                    prot_id = feature_map.getProteinIdentifications()[0]
                    prot_id.setIdentifier(f"Identifier_{i}")
                    prot_ids.append(prot_id)
                fm_new.setProteinIdentifications(prot_ids)
                for feature in feature_map:
                    pep_ids = []
                    for pep_id in feature.getPeptideIdentifications():
                        pep_id.setIdentifier(f"Identifier_{i}")
                        pep_ids.append(pep_id)
                    feature.setPeptideIdentifications(pep_ids)
                    fm_new.push_back(feature)
                feature_maps_mapped.append(fm_new)
    feature_maps = feature_maps_mapped

    '''Link features in a ConsensusMap.'''
    feature_grouper = oms.FeatureGroupingAlgorithmKD()

    consensus_map = oms.ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()

    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, oms.ColumnHeader())
        file_description.filename = os.path.basename(
            feature_map.getMetaValue("spectra_data")[0].decode()
        )
        file_description.size = feature_map.size()
        file_descriptions[i] = file_description

    feature_grouper.group(feature_maps, consensus_map)
    consensus_map.setColumnHeaders(file_descriptions)
    consensus_map.setUniqueIds()
    oms.ConsensusXMLFile().store(os.path.join(args.input_folder, "FeatureMatrix.consensusXML"), consensus_map)

    '''Create consensus file'''
    consensusXML_file = os.path.join(args.input_folder, "FeatureMatrix.consensusXML")
    consensus_map = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(consensusXML_file, consensus_map)
    filtered_map = oms.ConsensusMap(consensus_map)
    filtered_map.clear(False)
    for feature in consensus_map:
        if feature.getPeptideIdentifications():
            filtered_map.push_back(feature)

    consensusXML_file = os.path.join(args.input_folder, "filtered.consensusXML")
    oms.ConsensusXMLFile().store(consensusXML_file, filtered_map)


    '''Export for GNPS'''
    mgf_file = os.path.join(args.input_folder, "MS2data.mgf")
    quant_file = os.path.join(args.input_folder, "quant_ms2data.txt")
    oms.GNPSMGFFile().store(oms.String(consensusXML_file),[file.encode() for file in mzML_files],
        oms.String(mgf_file),)

    oms.GNPSQuantificationFile().store(consensus_map, quant_file)

if __name__ == '__main__':
    # input
    # test_spectra_file = '../msdb/data/input/240421-Kutz43850OSMAC-GSM.mzML'
    # output = '../msdb/data/kutz/'
    #
    # args = config.args
    # args.spectra_file = test_spectra_file
    # args.input_folder = '../msdb/data/input/'
    #
    # args.noise_level =1e4
    # args.allowed_mass_tolerance = 10.0
    # args.output = output
    # args.edbms2_file = '../msdb/edb_info.json' # format check{'CCMSID':...}
    # args.isms2_file = '../msdb/isdb_info.json'

    '''Preprocess of mzML file'''
    # preprocess_mzml(args)

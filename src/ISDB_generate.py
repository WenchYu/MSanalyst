
'''
Converting cfm-id output '.log' files to json files
Due to the large number of files, a demo is provided here('./data/xxx.log')

'''
import sys
sys.path.append('../')

import os,ast,time,json,re
import pandas as pd
import numpy as np
from tqdm import tqdm,trange
from my_packages import functions, cheminfo_tools

def ex_CFMIDspectra(file, start_txt, end_txt, skip_words=None):
    '''
    Horizontal and vertical coordinates of tandem mass
    :param file:
    :param start_txt:
    :param end_txt:
    :param skip_words:
    :return: A idlist contain lists of sliced content, like[[],[],...,[]],and converting to an array
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
            elif re.match(end_txt, lines[i]):
                spectrum = ''.join(lines[start_idx:i])
                spectra_list = spectrum.split('\n')[:-1]
                temp=[]
                for s in spectra_list:
                    text = s.split()
                    m_z, intensity = text[0],text[1]
                    temp.append([float(m_z), float(intensity)])
                temp = np.array(temp,dtype=np.float64)
                spectra.append(temp)
    return spectra

def cfmid_process(file):
    '''
    Preprocessing cfm-id output '.log' files
    :param file:
    :return:id,smiles,pepmass,energy0_ms2,energy1_ms2,energy2_ms2
    '''
    id_txt = '#ID='
    id = functions.ex_startswith(file, id_txt)[0]

    smiles_txt = '#SMILES='
    smiles = functions.ex_startswith(file, smiles_txt)[0]

    pepmass_txt = '#PMass='
    pepmass = [float(x) for x in functions.ex_startswith(file, pepmass_txt)][0]

    start_txt0 = 'energy0'
    end_txt0 = 'energy1'
    energy0_ms2 = ex_CFMIDspectra(file, start_txt0, end_txt0)

    start_txt1 = 'energy1'
    end_txt1 = 'energy2'
    energy1_ms2 = ex_CFMIDspectra(file, start_txt1, end_txt1)

    start_txt2 = 'energy2'
    end_txt2 = r'^\s*$'
    energy2_ms2 = ex_CFMIDspectra(file, start_txt2, end_txt2)

    return (id,smiles,pepmass,energy0_ms2,energy1_ms2,energy2_ms2)

def cfmid_json_preprocess(json_file):
    '''

    :param json_file:
    :return:
    '''
    with open(json_file, "r") as f:
        data = json.load(f)

    id = []
    smiles = []
    energy0_ms2 = []
    energy1_ms2 = []
    energy2_ms2 = []
    pepmass = []
    for i in trange(len(data)):
        id.append(data[i]['id'])
        smiles.append(data[i]['smiles'])
        pepmass.append(data[i]['pepmass'])
        energy0_ms2.append(ast.literal_eval(data[i]['energy0_ms2']))
        energy1_ms2.append(ast.literal_eval(data[i]['energy1_ms2']))
        energy2_ms2.append(ast.literal_eval(data[i]['energy2_ms2']))
    library_info = pd.DataFrame({
        'id': id,
        'smiles': smiles,
        'pepmass': pepmass,
        'energy0_ms2': energy0_ms2,
        'energy1_ms2': energy1_ms2,
        'energy2_ms2': energy2_ms2
    })
    print('The json file is parsed successfully！')
    return library_info

if __name__ =='__main__':
    t = time.time()

    '''1. CFM-ID .log file preprocessing'''
    cfmid_output_dir = '../msdb/GNPSLIBRARY/ISspec_demo/'
    cfmid_output_files = [cfmid_output_dir + f for f in os.listdir(cfmid_output_dir) if '.log' in f]
    data=[]
    for file in tqdm(cfmid_output_files,total=len(cfmid_output_files)):
        data.append(cfmid_process(file))

    '''isdb MS2 library'''
    ids, smiles = [],[]
    info_dict = {}
    for row in tqdm(data,total=len(data)):
        info_dict[row[0]] = {
            'smiles': row[1],
            'pepmass': row[2],
            'energy0_ms2': f'{row[3][0]}',
            'energy1_ms2': f'{row[4][0]}',
            'energy2_ms2': f'{row[5][0]}'
        }
        ids.append(row[0])
        smiles.append(row[1])
    with open('../msdb/GNPSLIBRARY/isdb_info.json', 'w') as f:
        json.dump(info_dict, f)

    '''isdb MS1 library'''
    MS1_library_df = pd.DataFrame({'id':ids,'smiles':smiles})
    MS1_library_df['formula'] = np.nan
    MS1_library_df['exactmass'] = np.nan
    for i in trange(len(MS1_library_df.index)):
        try:
            smile = MS1_library_df.smiles[i]
            MS1_library_df.loc[i, 'formula'] = cheminfo_tools.Smile2Formula(smile)
            MS1_library_df.loc[i, 'exactmass'] = cheminfo_tools.MyChemInfo.MolWt(MS1_library_df.formula[i])
            MS1_library_df.loc[i, 'm+h'] = MS1_library_df.loc[i, 'exactmass'] + 1.007276
            # MS1_library_df.loc[i, 'm+nh4'] = MS1_library_df.loc[i, 'exactmass'] + 18.033823
            # MS1_library_df.loc[i, 'm+na'] = MS1_library_df.loc[i, 'exactmass'] + 22.989218
        except:
            pass
    MS1_library_df.to_csv('../msdb/GNPSLIBRARY/isdbMS1.csv',index=None)

    print(f'Finish in {(time.time() - t) / 60:.2f}min')

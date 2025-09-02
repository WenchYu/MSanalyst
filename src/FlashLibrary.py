# -*- coding: utf-8 -*-
# @Time :2025/6/18 22:34
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
from FPSim2.io import create_db_file
from FPSim2 import FPSim2Engine
import os,sys,json,ms_entropy,matchms,json,time
sys.path.append('../')
import pandas as pd
import numpy as np
from my_packages import functions_new
import spectral_entropy as se
from ms_entropy.file_io import spec_file
from matchms.importing import load_from_mgf
from tqdm import tqdm, trange
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from sklearn.metrics import roc_curve,auc,confusion_matrix
from my_packages.peaktools import neutral_loss,modified_cosine,find_match_peaks_efficient,convert_to_peaks


if __name__ == '__main__':
    t = time.time()
    GNPS_LIBRARY_FILE = '../msdb/edb_info.json'
    with open(GNPS_LIBRARY_FILE, 'r') as f:
        GNPS_INFO = json.load(f)

    FS_GNPS_LIBRARY = []
    for key, values in tqdm(GNPS_INFO.items(), total=len(GNPS_INFO)):
        SPEC_STR, PM = values['ms2'], float(values['pepmass'])
        FS_GNPS_LIBRARY.append({
            "id": key,
            "precursor_mz": PM,
            "peaks": functions_new.spec_str2array(SPEC_STR, PM).tolist(),
            "smile": values['smiles'],
            "charge": values['charge'],
            "ion_mode": values['ion_mode']
        })

    FS_GNPS_LIBRARY_OUTPUT = '../msdb/FS_edb_info.json'
    with open(FS_GNPS_LIBRARY_OUTPUT, "w") as f:
        json.dump(FS_GNPS_LIBRARY, f)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')

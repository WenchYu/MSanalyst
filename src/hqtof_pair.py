# -*- coding: utf-8 -*-
# @Time :2025/7/31 23:36
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
# Import library
import sys,os,time,json,functools,argparse,heapq,pickle,ast,spectral_entropy
sys.path.append('../')
import numpy as np
from spectral_entropy import similarity
from tqdm import tqdm, trange
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from my_packages import functions_new
from matchms import Spectrum
from matchms.similarity import ModifiedCosine,NeutralLossesCosine,CosineGreedy

def get_spectra_from_einfo(EINFO,ID):
    '''

    :param einfo:
    :param id:
    :return:
    '''
    edb_pm = float(EINFO[ID]['PEPMASS'])
    ESMILE = EINFO[ID]['SMILE']
    temp = EINFO[ID]['MS2_SPECTRUM'] # <list>[[],[],...]
    if type(temp) == str:
        ESPEC = np.asarray(ast.literal_eval(temp)) # "[[],[],[],...]"
    else:
        ESPEC = temp
    ESPEC = spectral_entropy.clean_spectrum(ESPEC, max_mz=edb_pm+0.01)
    ESPECTRUM = Spectrum(mz=np.array(ESPEC[:, 0],dtype = float),
                  intensities=np.array(ESPEC[:, 1],dtype = float),
                  metadata={"precursor_mz": edb_pm,"smile":ESMILE})

    return ESPEC,ESPECTRUM

def test():
    # with open('../msdb/data/hqtof/compare_dict.pkl', 'rb') as f:
    #     compare_dict = pickle.load(f)
    #
    # HQTOF = functions_new.json_load('../msdb/FS_hqtof.json')
    # hqtof_search = FlashEntropySearch()
    # FS_SPECTRA = hqtof_search.build_index(HQTOF)  # Pre-clean and sorted flash spectra
    #
    # # Search
    # SEARCH_ALGOs = ['entropy']
    #
    #
    # for SEARCH_ALGO in SEARCH_ALGOs:
    #     print(SEARCH_ALGO)
    #     EMS2_RESULT = {}
    #     EMS2_SIM_MATRIX = np.zeros((len(FS_SPECTRA), len(FS_SPECTRA)), dtype=float)  # Initialize empty
    #     EMS2_PEAK_MATRIX = np.zeros((len(FS_SPECTRA), len(FS_SPECTRA)), dtype=float)
    #
    #     for query_idx, E_idxs in tqdm(compare_dict.items(), total=len(compare_dict)):
    #         query_PM = FS_SPECTRA[query_idx]['precursor_mz']
    #         query_PEAKS = FS_SPECTRA[query_idx]['peaks']
    #
    #         spectrum_1 = Spectrum(mz=np.array(query_PEAKS[:, 0], dtype=float),
    #                               intensities=np.array(query_PEAKS[:, 1], dtype=float),
    #                               metadata={"precursor_mz": query_PM})
    #
    #         for E_idx in E_idxs:  # E_idxs
    #             E_PM = FS_SPECTRA[E_idx]['precursor_mz']
    #             E_PEAKS = FS_SPECTRA[E_idx]['peaks']
    #             spectrum_2 = Spectrum(mz=np.array(E_PEAKS[:, 0], dtype=float),
    #                                   intensities=np.array(E_PEAKS[:, 1], dtype=float),
    #                                   metadata={"precursor_mz": E_PM})
    #
    #             if SEARCH_ALGO == 'modified_cosine':
    #                 modified_cosine = ModifiedCosine(tolerance=0.02)
    #                 score = modified_cosine.pair(spectrum_1, spectrum_2)
    #                 SPEC_SIM = score['score']
    #                 N_PEAK = score['matches']
    #
    #                 EMS2_SIM_MATRIX[query_idx, E_idx] = SPEC_SIM
    #                 EMS2_PEAK_MATRIX[query_idx, E_idx] = N_PEAK
    #
    #             elif SEARCH_ALGO == 'cosine':
    #                 cosine = CosineGreedy(tolerance=0.02)
    #                 score = cosine.pair(spectrum_1, spectrum_2)
    #
    #                 SPEC_SIM = score['score']
    #                 N_PEAK = score['matches']
    #
    #                 EMS2_SIM_MATRIX[query_idx, E_idx] = SPEC_SIM
    #                 EMS2_PEAK_MATRIX[query_idx, E_idx] = N_PEAK
    #
    #             elif SEARCH_ALGO == 'neutral_loss':
    #                 neutral_loss = NeutralLossesCosine(0.02)
    #                 score = neutral_loss.pair(spectrum_1, spectrum_2)
    #
    #                 SPEC_SIM = score['score']
    #                 N_PEAK = score['matches']
    #
    #                 EMS2_SIM_MATRIX[query_idx, E_idx] = SPEC_SIM
    #                 EMS2_PEAK_MATRIX[query_idx, E_idx] = N_PEAK
    #
    #                 # # if functions_new.calculate_ppm(query_PM,E_PM) <= 10:
    #             else:
    #
    #                 SPEC_SIM = similarity(query_PEAKS, E_PEAKS, method=SEARCH_ALGO, ms2_da=0.02)
    #                 EMS2_SIM_MATRIX[query_idx, E_idx] = SPEC_SIM
    #
    #     np.save(f'../msdb/data/hqtof/SpecSimMatrix/{SEARCH_ALGO}.npy', EMS2_SIM_MATRIX)
    #     np.save(f'../msdb/data/hqtof/SpecSimMatrix/{SEARCH_ALGO}_mp.npy', EMS2_PEAK_MATRIX)
    return None

if __name__ == '__main__':
    t = time.time()
    modified_cosine = ModifiedCosine(tolerance=0.02)
    cosine = CosineGreedy(tolerance=0.02)
    neutral_loss = NeutralLossesCosine(0.02)

    '''******************************************************************************************'''
    # Load
    DIR = f'../msdb/data/hqtof/'  # former: f'../msdb/data/hqtof/matrix/'
    Hqtof_CCMSIDs = list(
        np.load(os.path.join(DIR, 'idlist/H_qtof_non-redundant_CCMSIDs.npy')))  # [CCMSID1,CCMSID2,CCMSID3, ...]
    GNPS_INFO = functions_new.json_load('../msdb/GNPSLIBRARY/GNPS-LIBRARY-INFO.json')

    # Calc
    MP_MATRIX1 = np.load(os.path.join(DIR, f'SpecSimMatrix/modified_cosine0.npy'))  #
    SPEC_SIM_MATRIX1 = np.load(os.path.join(DIR, f'SpecSimMatrix/modified_cosine0.npy'))  #

    ALGO = 'modified_cosine'
    print(ALGO)
    for idx1 in trange(len(Hqtof_CCMSIDs)):
        CCMSID1 = Hqtof_CCMSIDs[idx1]
        SPEC1, SPECTRUM1 = get_spectra_from_einfo(GNPS_INFO, CCMSID1)

        for idx2 in range(idx1):
            CCMSID2 = Hqtof_CCMSIDs[idx2]
            SPEC2, SPECTRUM2 = get_spectra_from_einfo(GNPS_INFO, CCMSID2)
            if SPEC_SIM_MATRIX1[idx1, idx2] > 0:
                SPEC_SIM, MP = functions_new.clac_spec_sim(SPEC1, SPECTRUM1, SPEC2, SPECTRUM2, ALGO)
                SPEC_SIM_MATRIX1[idx1, idx2] = SPEC_SIM
                SPEC_SIM_MATRIX1[idx2, idx1] = SPEC_SIM
                MP_MATRIX1[idx1, idx2] = MP
                MP_MATRIX1[idx2, idx1] = MP

    np.save(f'../msdb/data/hqtof/SpecSimMatrix/{ALGO}.npy', SPEC_SIM_MATRIX1)
    np.save(f'../msdb/data/hqtof/SpecSimMatrix/{ALGO}_MP.npy', MP_MATRIX1)



    print(f'Finished in {(time.time() - t) / 60:.2f} min')

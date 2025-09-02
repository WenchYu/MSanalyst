# -*- coding: utf-8 -*-
# @Time :2025/6/5 18:19
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import sys
sys.path.append('../')
import time,json, functools, argparse,spectral_entropy
import numpy as np
from my_packages import functions_new, peaktools
from tqdm import tqdm,trange
from concurrent.futures import ProcessPoolExecutor



def parse_args():
    parse = argparse.ArgumentParser(description='Passsing mass spectral algorithms') # Create parameter object
    parse.add_argument('-a','--algorithm',type=str, help='Similarity algorithm of tandem mass matching used for library search'
                             , default='cosine') # add parameter to object
    args = parse.parse_args() # parse parameter object to get parse object
    return args

def spectral_sim_calc(gnps_info, CCMSIDs, ALGORITHM, args):
    '''
    Calculating spectral simialrity parallelly
    Only cosine, modified_cosine, neutral_loss print number of shared peaks
    :param gnps_info:
    :param CCMSIDs: idlist for calculating
    :param args: idx1, idx2 correspond with indices
    :return: indices, spectral simialrity and n_peaks
    '''
    idx1,idx2 = args
    spec1,spectrum1 = functions_new.GNPS_info_format(gnps_info,CCMSIDs[idx1])
    spec2, spectrum2 = functions_new.GNPS_info_format(gnps_info, CCMSIDs[idx2])

    if ALGORITHM == 'cosine':
        try:
            RESULT = peaktools.cosine(spectrum1,spectrum2, 0.05)
            SIMILARITY = round(RESULT.score,2)
            N_PEAKS = RESULT.matches
        except ZeroDivisionError:
            SIMILARITY,N_PEAKS = 0.0,0
        return idx1, idx2, SIMILARITY,N_PEAKS

    if ALGORITHM == 'modified_cosine':
        try:
            RESULT = peaktools.modified_cosine(spectrum1, spectrum2, 0.05)
            SIMILARITY = round(RESULT.score)
            N_PEAKS = RESULT.matches
        except ZeroDivisionError:
            SIMILARITY,N_PEAKS = 0,0

        return idx1, idx2, SIMILARITY,N_PEAKS
    if ALGORITHM == 'neutral_loss':
        try:
            RESULT = peaktools.neutral_loss(spectrum1, spectrum2, 0.05)
            SIMILARITY = round(RESULT.score,2)
            N_PEAKS = RESULT.matches
        except ZeroDivisionError:
            SIMILARITY, N_PEAKS = 0, 0
        return idx1, idx2, SIMILARITY, N_PEAKS

    else:
        try:
            SIMILARITY = round(spectral_entropy.similarity(spec1,spec2,method=ALGORITHM,ms2_da=0.05),4)

        except ZeroDivisionError:
            SIMILARITY = 0
        return idx1, idx2, SIMILARITY, 0

def main():
    # Loading processed json file
    GNPS_JSON = '../msdb/GNPSLIBRARY_250514/GNPS-LIBRARY-INFO.json'
    with open(GNPS_JSON, 'r') as f:
        gnps_info = json.load(f)

    # Set parameters
    CCMSIDs = np.load('../msdb/data/idlist/H_qtof_non-redundant_CCMSIDs.npy').tolist()
    indices = list(range(0,len(CCMSIDs)))
    ALGORITHM = 'cosine'
    tasks = []
    for i in trange(len(indices)):
        for j in range(i, len(indices)):
            tasks.append((indices[i], indices[j]))

    # Calculate similarity matrices
    similarity_matrix = np.zeros((len(CCMSIDs), (len(CCMSIDs))))
    peak_matrix = np.zeros((len(CCMSIDs), (len(CCMSIDs))))
    with ProcessPoolExecutor(max_workers=4) as executor:
        func = functools.partial(spectral_sim_calc, gnps_info, CCMSIDs, ALGORITHM)
        for result in tqdm(executor.map(func, tasks), total=len(tasks)):
            idx1, idx2, sim, n_peaks = result
            similarity_matrix[idx1, idx2] = sim
            similarity_matrix[idx2, idx1] = sim
            peak_matrix[idx1, idx2] = n_peaks
            peak_matrix[idx2, idx1] = n_peaks


    # Saving similarity matrices
    np.save(f'{ALGORITHM}_similarity_matrix.npy', similarity_matrix)
    np.save(f'{ALGORITHM}_peak_matrix.npy', peak_matrix)

    # Saving indices <str> : CCMSIDs
    INDEXtoCCMSID = {index: CCMSID for index, CCMSID in enumerate(CCMSIDs)}
    with open (f'{ALGORITHM}_INDEXtoCCMS.json', 'w') as f2:
        json.dump(INDEXtoCCMSID,f2)

    # # Accessing the similarity matrices
    # ALGORITHM = 'entropy'
    # with open(f'{ALGORITHM}_INDEXtoCCMS.json', 'r') as f3:
    #     INDEXtoCCMS_info = json.load(f3)
    # print(INDEXtoCCMS_info['2']) # Str type
    #
    # similarity_matrix_file = f'{ALGORITHM}_similarity_matrix.npy'
    # similarity_matrix = np.load(similarity_matrix_file)
    # SPECSIM = similarity_matrix[0, 1]
    # peak_matrix_file = f'{ALGORITHM}_peak_matrix.npy'
    # peak_matrix = np.load(peak_matrix_file)
    # NPEAKS = peak_matrix[0, 1]
    # print(type(SPECSIM), type(NPEAKS))


if __name__ == '__main__':
    t = time.time()
    # parse args
    args = parse_args()

    # Load processed json file
    GNPS_JSON = '../msdb/GNPSLIBRARY_250514/GNPS-LIBRARY-INFO.json'
    with open(GNPS_JSON, 'r') as f:
        gnps_info = json.load(f)

    # Load filtered CCMSIDs
    CCMSIDs = np.load('../msdb/data/idlist/H_qtof_non-redundant_CCMSIDs.npy').tolist()

    # Calculate matrices
    similarity_matrix = np.zeros((len(CCMSIDs), (len(CCMSIDs))))
    peak_matrix = np.zeros((len(CCMSIDs), (len(CCMSIDs))))
    ALGORITHM = args.algorithm
    for idx1 in trange(len(CCMSIDs)):
        spec1, spectrum1 = functions_new.GNPS_info_format(gnps_info, CCMSIDs[idx1])
        for idx2 in range(idx1):
            spec2, spectrum2 = functions_new.GNPS_info_format(gnps_info, CCMSIDs[idx2])
            if ALGORITHM == 'cosine':
                try:
                    RESULT = peaktools.cosine(spectrum1, spectrum2, 0.05)
                    SIMILARITY = round(RESULT.score, 2)
                    N_PEAKS = RESULT.matches
                except ZeroDivisionError:
                    SIMILARITY, N_PEAKS = 0, 0

            if ALGORITHM == 'modified_cosine':
                try:
                    RESULT = peaktools.modified_cosine(spectrum1, spectrum2, 0.05)
                    SIMILARITY = round(RESULT.score,2)
                    N_PEAKS = RESULT.matches
                except ZeroDivisionError:
                    SIMILARITY, N_PEAKS = 0, 0

            if ALGORITHM == 'neutral_loss':
                try:
                    RESULT = peaktools.neutral_loss(spectrum1, spectrum2, 0.05)
                    SIMILARITY = round(RESULT.score, 2)
                    N_PEAKS = RESULT.matches
                except ZeroDivisionError:
                    SIMILARITY, N_PEAKS = 0, 0

            else:
                try:
                    SIMILARITY = round(spectral_entropy.similarity(spec1, spec2, method=ALGORITHM, ms2_da=0.05), 2)

                except ZeroDivisionError:
                    SIMILARITY, N_PEAKS = 0, 0

            similarity_matrix[idx1, idx2] = SIMILARITY
            similarity_matrix[idx2, idx1] = SIMILARITY
            peak_matrix[idx1, idx2] = N_PEAKS
            peak_matrix[idx2, idx1] = N_PEAKS

    # Saving similarity matrices
    np.save(f'{ALGORITHM}_similarity_matrix.npy', similarity_matrix)
    np.save(f'{ALGORITHM}_peak_matrix.npy', peak_matrix)

    # Saving indices <str> : CCMSIDs
    INDEXtoCCMSID = {index: CCMSID for index, CCMSID in enumerate(CCMSIDs)}
    with open (f'{ALGORITHM}_INDEXtoCCMS.json', 'w') as f2:
        json.dump(INDEXtoCCMSID,f2)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')

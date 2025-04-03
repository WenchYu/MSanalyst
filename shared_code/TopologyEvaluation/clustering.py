# -*- coding: utf-8 -*-
# @Time :2025/3/29 22:37
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import time
from my_packages import ms2tools,config
from tqdm import tqdm

if __name__ == '__main__':
    t = time.time()
    args = config.args
    mgf = './data/std.mgf'
    csv = './std_quant.csv'

    args.output = './'
    args.quant_file = csv
    args.mgf_file = mgf
    args.self_clustering_similarity = 0.1
    args.self_clustering_peaks = 1
    meds = ['dot_product','modified_cosine','dot_product_reverse'
,'fidelity_distance','weighted_dot_product_distance','inner_product'
,'spectral_contrast_angle_distance','baroni_urbani_buser'
,'hattacharya_1','bhattacharya_2','harmonic_mean','motyka','intersection'
,'neutral_loss','unweighted_entropy','entropy','ms_for_id','ms_for_id_v1'
,'euclidean', 'squared_euclidean', 'clark', 'dice','divergence'
, 'squared_chord_distance', 'jaccard', 'matusita_distance'
, 'hellinger', 'symmetric_chi_squared', 'improved_similarity_distance'
, 'probabilistic_symmetric_chi_squared', 'vicis_symmetric_chi_squared_3'
,'pearson_correlation', 'avg_l', 'canberra', 'lorentzian', 'manhattan'
, 'chebyshev', 'mean_character_distance', 'wave_hedges_distance', 'penrose_shape', 'penrose_size', 'absolute_value'
,'whittaker_index_of_association_distance','roberts', 'ruzicka']
    for med in tqdm(meds,total=len(meds)):
        args.self_clustering_method = med
        ms2tools.self_clustering(args)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')

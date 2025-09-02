
'''
MSanalyst annotating non-targeted metabolomics data includes the following steps:
1. MS1 match conducted by 'ms1_match()'
2. 'ISDB_MS2_match()' and 'EDB_MS2_match()' will perform MS2 comparison based on the MS1 match results
3. Query features are first clustered and curated by 'self_clustering()'
4. The SpecSimNetwork will be annotated by filtered MS2 comparison results in 'molecular_generation()'
'''

import os
import ast,time,json,heapq,spectral_entropy
import pandas as pd
import numpy as np
import networkx as nx
import spectrum_utils.spectrum as sus
from rdkit import Chem
from tqdm import tqdm, trange
from joblib import Parallel, delayed
from spectral_entropy import similarity
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from matchms.similarity import ModifiedCosine,NeutralLossesCosine,CosineGreedy
from matchms import Spectrum
from matchms.exporting import save_as_mgf

from my_packages import functions,config,functions_new,evaluation,cheminfo_tools

class FlashPrecursorSearch:
    def __init__(
            self,
            max_ms2_tolerance_in_da=0.02,
            mz_index_step=0.0001,
            adduct='m+h') -> None:
        """
        Initialize the EntropySearch class.

        :param max_ms2_tolerance_in_da: The maximum MS2 tolerance used when searching the MS/MS spectra, in Dalton. Default is 0.024.
        :param mz_index_step:   The step size of the m/z index, in Dalton. Default is 0.0001.
                                The smaller the step size, the faster the search, but the larger the index size and longer the index building time.
        """
        self.mz_index_step = mz_index_step
        self.max_ms2_tolerance_in_da = max_ms2_tolerance_in_da
        self.adduct = adduct

        self.index = []  # record ions_mz

        self.index_names = [
            "all_ions_mz_idx_start",
            "all_ions_mz",
        ]
        self.index_dtypes = {
            "all_ions_mz_idx_start": np.int64,
            "all_ions_mz": np.float32,
        }

    def _generate_index_from_peak_data(self, dataframe, max_indexed_mz):
        # Sort with precursor m/z and pre-sort

        # Record the precursor m/z
        all_ions_mz = dataframe[self.adduct]

        # Build index for fast access to the ion's m/z.
        max_mz = max_indexed_mz
        search_array = np.arange(0.0, max_mz, self.mz_index_step)
        all_ions_mz_idx_start = np.searchsorted(all_ions_mz, search_array, side="left").astype(np.int64)

        ############## Step 4: Save the index. ##############
        index = [
            all_ions_mz_idx_start,
            all_ions_mz,
        ]
        return index

    def build_index(self, dataframe, max_indexed_mz: float = 1500.00005):
        """
        Build the index for the MS/MS spectra library.

        The spectra provided to this function should be a dictionary in the format of {"precursor_mz": precursor_mz, "peaks": peaks}.
        The precursor_mz is the precursor m/z value of the MS/MS spectrum;
        The peaks is a numpy array which has been processed by the function "clean_spectrum".

        :param all_spectra_list:    A list of dictionaries in the format of {"precursor_mz": precursor_mz, "peaks": peaks},
                                    the spectra in the list need to be sorted by the precursor m/z.
        :param max_indexed_mz: The maximum m/z value that will be indexed. Default is 1500.00005.
        """

        ############## Step 2: Build the index by sort with product ions. ##############
        self.index = self._generate_index_from_peak_data(dataframe, max_indexed_mz)
        return self.index

    def _find_location_from_array_with_index(self, wanted_mz, mz_array, mz_idx_start_array, side,
                                             index_number_in_one_da):
        mz_min_int = (np.floor(wanted_mz * index_number_in_one_da)).astype(int)
        mz_max_int = mz_min_int + 1

        if mz_min_int >= len(mz_idx_start_array):
            mz_idx_search_start = mz_idx_start_array[-1]
        else:
            mz_idx_search_start = mz_idx_start_array[mz_min_int].astype(int)

        if mz_max_int >= len(mz_idx_start_array):
            mz_idx_search_end = len(mz_array)
        else:
            mz_idx_search_end = mz_idx_start_array[mz_max_int].astype(int) + 1

        return mz_idx_search_start + np.searchsorted(mz_array[mz_idx_search_start:mz_idx_search_end], wanted_mz,
                                                     side=side)

    def _get_indeces(self, mz_query):
        index_number_in_one_da = int(1/self.mz_index_step)
        all_ions_mz = self.index[1]
        all_ions_mz_idx_start = self.index[0]
        product_mz_idx_min = self._find_location_from_array_with_index(
            mz_query - self.max_ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "left",
            index_number_in_one_da)

        product_mz_idx_max = self._find_location_from_array_with_index(
            mz_query + self.max_ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "right",
            index_number_in_one_da)

        return product_mz_idx_min, product_mz_idx_max

def match_mz(quant_df_row, msdb_df, mz_column='row m/z',isdb_ms1_match_threshld = 5):
    '''
    Single MS1 match against the in-silico MS1 library

    :param quant_df_row: Rows in quantification table(xxx_quant.csv) generated by MZmine
    :param msdb_df: The default database is (default: '../msdb/isdbMS1.csv')
    :param mz_column: Column name of queried MS1 (default: 'row m/z')
    :param np_ms1_match_threshold: Relative error (ppm)
    :return:tuple([...],[...],[...])
    '''
    hits_id = []
    hits_smiles = []
    for j in range(len(msdb_df.index)):
        ppm_H = functions.calculate_ppm(quant_df_row[mz_column], msdb_df.loc[j,'m+h'])

        if ppm_H < isdb_ms1_match_threshld :
            hits_id.append(msdb_df.id[j])
            hits_smiles.append(msdb_df.smiles[j])
    if not hits_id or not hits_smiles:
        hits_id.append(None)
        hits_smiles.append(None)
    return hits_id, hits_smiles

def match_edb_mz(quant_df_row,edb_df,mz_column='row m/z',edb_ms1_match_threshold = 5):
    '''
    Single MS1 match against experimental MS1 library
    '''
    hits_id = []
    hits_smiles = []
    for j in range(len(edb_df.index)):
        ppm = functions.calculate_ppm(quant_df_row[mz_column], edb_df['pepmass'][j])
        if ppm < edb_ms1_match_threshold:
            hits_id.append(edb_df.id[j])
            hits_smiles.append(str(edb_df.smiles[j]))
    if not hits_id or not hits_smiles:
        hits_id.append(None)
        hits_smiles.append(None)
    return hits_id, hits_smiles

def ms1_match(args,queryMGF = None):
    '''
    MS1 match against the entire MS1 library (experimental and in-silico)
    :param args: args.pepmass_match_tolerance,
    :return: Two MS1 match result files, prefixed with 'IS_MS1match_' and 'E_MS1match_'
    '''

    if queryMGF is None:
        SPECTRA = functions_new.load_spectra_from_file(args.spectra_file)
        FS_SPECTRA = functions_new.spectra_to_fsspectra(SPECTRA)  # Load input spectra
        quant_df = pd.DataFrame(
            {'row ID': [s['id'] for s in FS_SPECTRA], 'row m/z': [s['precursor_mz'] for s in FS_SPECTRA]})

        basename = os.path.splitext(os.path.basename(args.spectra_file))[0]
        result_dir = os.path.join(args.output, f'{basename}_result')
        os.makedirs(result_dir, exist_ok=True)


    # else: # Single MS1 search
    #     quant_df = queryMGF  # Query dataframe from single pepmass (used in '../ms1search.py')
    #     query_mz = str(queryMGF.iloc[0, 1])
    #     result_dir = os.path.join(args.output, f'{query_mz}')
    #     os.makedirs(result_dir, exist_ok=True)
    #     basename = f'{query_mz}.csv'

        # MS1 match against the entire library (experimental and in-silico)
    n_jobs = args.cpus

    '''Search against experimental library (GNPS)'''
    edb_df = functions.df_preprocess(args.edbms1_file)
    edb_results = Parallel(n_jobs=n_jobs)(
        delayed(match_edb_mz)(quant_df_row, edb_df, edb_ms1_match_threshold=args.allowed_mass_tolerance) for
        quant_df_row in
        tqdm(quant_df.to_dict('records')))

    edb_match_rows = []
    for i, (hits_id, hits_smiles) in enumerate(edb_results):
        for j in range(len(hits_id)):
            edb_match_row = {'row ID': quant_df.at[i, 'row ID'], 'row m/z': quant_df.at[i, 'row m/z'],
                             'match_id': hits_id[j], 'match_smiles': hits_smiles[j]}
            edb_match_rows.append(edb_match_row)
    edb_match_df = pd.DataFrame(edb_match_rows)
    edb_result_path = os.path.join(result_dir, f'E_MS1match_{basename}.csv')
    edb_match_df.to_csv(edb_result_path, index=False)

    quant_df['edbms1_id'], quant_df['edbms1_smiles'] = np.nan, np.nan

    for i, (hits_id, hits_smiles) in enumerate(edb_results):
        if hits_id is not None and all(isinstance(x, str) for x in hits_id):
            quant_df.at[i, 'edbms1_id'] = ';'.join([x or '' for x in hits_id])
            quant_df.at[i, 'edbms1_smiles'] = ';'.join(hits_smiles)

    ms1_result_path = os.path.join(result_dir, f'MS1match_{basename}.csv')
    quant_df.to_csv(ms1_result_path, index=False)

    '''Search against in-silico library'''
    if args.isms1_file:
        isdb_df = pd.read_csv(args.isms1_file, low_memory=False)  # Loading MS1 library
        isdb_results = Parallel(n_jobs=n_jobs)(
            delayed(match_mz)(quant_df_row, isdb_df, isdb_ms1_match_threshld=args.allowed_mass_tolerance) for quant_df_row
            in
            tqdm(quant_df.to_dict('records')))

        quant_df['isdbms1_id'],quant_df['isdbms1_smiles'] = np.nan,np.nan
        for i, (hits_id, hits_smiles) in enumerate(isdb_results):
            if hits_id is not None and all(isinstance(x, str) for x in hits_id):
                quant_df.at[i, 'isdbms1_id'] = ';'.join([x or '' for x in hits_id])
                quant_df.at[i, 'isdbms1_smiles'] = ';'.join(hits_smiles)

        isdb_match_rows = []
        for i, (hits_id, hits_smiles) in enumerate(isdb_results):
            for j in range(len(hits_id)):
                isdb_match_row = {'row ID': quant_df.at[i, 'row ID'], 'row m/z': quant_df.at[i, 'row m/z'],
                                  'match_id': hits_id[j], 'match_smiles': hits_smiles[j]}
                isdb_match_rows.append(isdb_match_row)

        isdb_match_df = pd.DataFrame(isdb_match_rows)
        isdb_result_path = os.path.join(result_dir, f'IS_MS1match_{basename}.csv')
        isdb_match_df.to_csv(isdb_result_path, index=False)

    print('MS1 matching finished!')

def ems2_match(args,queryMGF=None):
    ''' MS2 match against experimental MS2 library '''
    if queryMGF is None:

        exp_info = functions_new.spectra_process(args.spectra_file)  # Load input spectra
        basename = os.path.splitext(os.path.basename(args.spectra_file))[0]
        result_dir = os.path.join(args.output, f'{basename}_result')
        os.makedirs(result_dir, exist_ok=True)

    else:
        spectra_info = queryMGF # Query MS2 from direct input (used in '../ms2search.py')
        query_mz = spectra_info.loc[0,'pepmass']
        result_dir = os.path.join(args.output, f'{query_mz}')
        os.makedirs(result_dir, exist_ok=True)
        basename = f'{query_mz}.csv'
        exp_info = queryMGF

    edbms2_info = functions_new.json_load(args.edbms2_file) # Loading experimental MS2 library

    # Loading experimental MS1 match result (csv file generated by ms1_match())
    ems1_result_file = os.path.join(result_dir, f'E_MS1match_{basename}.csv')
    edb_ms1_df = functions.df_preprocess(ems1_result_file) # query—

    col_name = f'{args.library_matching_method}_similarity'
    if col_name in edb_ms1_df.columns: # Skip if the library match method has been applied
        pass
    else:
        for i in trange(len(edb_ms1_df)):
            row_id = str(edb_ms1_df.loc[i,'row ID'])
            match_id = edb_ms1_df.loc[i,'match_id']

            try: # Some query spectra have no ms1 hits
                SPEC, SPECTRUM = functions_new.get_spectra_from_info(exp_info, row_id)
                ESPEC,ESPECTRUM= functions_new.get_spectra_from_einfo(edbms2_info,match_id)
                SPEC_SIM, N_PEAK = functions_new.clac_spec_sim(SPEC,SPECTRUM,ESPEC,ESPECTRUM,args.library_matching_method)

                # Output
                edb_ms1_df.loc[i, f'{args.library_matching_method}_similarity'] = SPEC_SIM
                edb_ms1_df.loc[i, f'{args.library_matching_method}_mps'] = N_PEAK
                edb_ms2_path = os.path.join(result_dir, row_id, f'{match_id}.mgf')
                # save_as_mgf(ESPECTRUM, edb_ms2_path) # save matched mgf

            except:pass
    edb_ms1_df.to_csv(ems1_result_file,index = None)
    print('EMS2 matching finished!')

def isms2_match(args,queryMGF=None):
    ''' MS2 match against in-silico MS2 library '''

    if queryMGF is None:

        exp_info = functions_new.spectra_process(args.spectra_file)  # Load input spectra
        basename = os.path.splitext(os.path.basename(args.spectra_file))[0]
        result_dir = os.path.join(args.output, f'{basename}_result')
        os.makedirs(result_dir, exist_ok=True)

    else:
        spectra_info = queryMGF # Query MS2 from direct input (used in '../ms2search.py')
        query_mz = spectra_info.loc[0,'pepmass']
        result_dir = os.path.join(args.output, f'{query_mz}')
        os.makedirs(result_dir, exist_ok=True)
        basename = f'{query_mz}.csv'
        exp_info = queryMGF



    isdb_info = functions_new.json_load(args.isms2_file) # Loading in-silico MS2 library
    # Loading in-silico MS1 match result (csv file generated by ms1_match())
    is_result_path = os.path.join(result_dir, f'IS_MS1match_{basename}.csv')
    is_ms1_match_df = functions.df_preprocess(is_result_path)

    # 0, 1, 2 represent the in-silico MS2 predicted by CFM-ID at three collision energy (10, 20 and 40 eV)
    # cols_to_add = [
    #     f'{args.library_matching_method}_similarity0',
    #     f'{args.library_matching_method}_mps0',
    #     f'{args.library_matching_method}_similarity1',
    #     f'{args.library_matching_method}_mps1',
    #     f'{args.library_matching_method}_similarity2',
    #     f'{args.library_matching_method}_mps2'
    # ]
    # for c in cols_to_add:
    #     if c not in is_ms1_match_df.columns:
    #         is_ms1_match_df[c] = np.nan

    for i in trange(len(is_ms1_match_df)):
        row_id = str(is_ms1_match_df.loc[i,'row ID']) # convert int to str
        match_id = str(is_ms1_match_df.loc[i,'match_id'])
        if match_id != 'nan':
            try:  # Some features in xxx_quant.csv have no MS2
                SPEC, SPECTRUM = functions_new.get_spectra_from_info(exp_info, row_id)

                e0_SPEC,e0_SPECTRUM,e1_SPEC,e1_SPECTRUM,e2_SPEC,e2_SPECTRUM = functions_new.get_spectra_from_isinfo(isdb_info,match_id)

                SPEC_SIM0, N_PEAK0 = functions_new.clac_spec_sim(SPEC, SPECTRUM, e0_SPEC, e0_SPECTRUM,
                                                                 args.library_matching_method)
                SPEC_SIM1, N_PEAK1 = functions_new.clac_spec_sim(SPEC, SPECTRUM, e1_SPEC, e1_SPECTRUM,
                                                                 args.library_matching_method)
                SPEC_SIM2, N_PEAK2 = functions_new.clac_spec_sim(SPEC, SPECTRUM, e2_SPEC, e2_SPECTRUM,
                                                                 args.library_matching_method)

            except:
                SPEC_SIM0, SPEC_SIM1, SPEC_SIM2 = 0.0, 0.0, 0.0
                N_PEAK0, N_PEAK1, N_PEAK2 = 0, 0, 0

            # Output
            is_ms1_match_df.loc[i, f'{args.library_matching_method}_similarity0'] = SPEC_SIM0
            is_ms1_match_df.loc[i, f'{args.library_matching_method}_mps0'] = N_PEAK0
            is_ms1_match_df.loc[i, f'{args.library_matching_method}_similarity1'] = SPEC_SIM1
            is_ms1_match_df.loc[i, f'{args.library_matching_method}_mps1'] = N_PEAK1
            is_ms1_match_df.loc[i, f'{args.library_matching_method}_similarity2'] = SPEC_SIM2
            is_ms1_match_df.loc[i, f'{args.library_matching_method}_mps2'] = N_PEAK2
            is_ms2_path = os.path.join(result_dir, row_id, f'{match_id}.mgf')
            # save_as_mgf(,is_ms2_path)

    is_ms1_match_df.to_csv(is_result_path, index = None)
    print('ISMS2 matching finished!')

def mn_curating(G: nx.Graph, topk: int) -> nx.Graph:
    """
    Limit the number of neighbors of each node to no more than topk,
    retaining the topk edges with the highest 'pair_similarity'.

    """
    for node in list(G.nodes()):
        neighbors = list(G[node]) #  Current neighbor list
        if len(neighbors) <= topk:
            continue

        # Get all (similarity, neighbors) and sort
        scored_nb = [(G[node][nb].get('pair_similarity', 0), nb) for nb in neighbors]
        topk_nb   = {nb for _, nb in heapq.nlargest(topk, scored_nb)}

        # Precompute edges that need to be deleted
        edges_to_remove = [(node, nb) for nb in neighbors if nb not in topk_nb]

        G.remove_edges_from(edges_to_remove)
    return G

def self_clustering(args,exp_info):
    '''

    :param args: args.output, args.quant_file, args.spectra_file, args.self_clustering_similarity
    :return:
    '''
    G = nx.MultiGraph()  # Creating undirected graph
    for key, values in exp_info.items():
        pm1 = float(values['pepmass'])
        node_attr = {'pepmass': pm1,'level':'unmatch','class':'FEATURE'}
        G.add_node(key, **node_attr)  # add nodes and attributes

    # Self clustering
    feature_ids = [key for key, values in exp_info.items()]
    for i in trange(len(feature_ids)):
        FID1 = feature_ids[i]
        SPEC1, SPECTRUM1 = functions_new.get_spectra_from_info(exp_info, FID1)

        for j in range(i):
            FID2 = feature_ids[j]
            SPEC2, SPECTRUM2 = functions_new.get_spectra_from_info(exp_info, FID2)

            # SPEC_SIM, N_PEAK = 0.0, 0
            SPEC_SIM, N_PEAK = functions_new.clac_spec_sim(SPEC1, SPECTRUM1, SPEC2, SPECTRUM2,
                                                           args.self_clustering_method)

            if SPEC_SIM >= args.self_clustering_similarity \
                    and N_PEAK >= args.self_clustering_peaks:
                edge_attr = {f'{args.self_clustering_method}_pair_similarity': SPEC_SIM,
                             'pair_matched_peaks': N_PEAK,
                             'edge_type': args.self_clustering_method}
                G.add_edge(FID1, FID2, **edge_attr)

    G = mn_curating(G, args.top_k)
    print('Self clustering finished!')
    return G

def ismatch_filter(G):
    '''
    1. Keep IS hits share chem dice sim (≥ 0.7) with experimental hits
    2. If only IS hits were found in cluster, hits were kept at least 70% portion share chem dice sim (≥ 0.7)
    :param G:
    :return:
    '''
    clusters, cluster_nodes, cluster_map = evaluation.extract_cluster(G)
    for cluster in clusters:
        E_hits = [node for node in cluster if G.nodes[node]['class'] == 'EDB']
        IS_hits = [node for node in cluster if G.nodes[node]['class'] == 'ISDB']

        # Step 1
        if E_hits and IS_hits:  # Both E and IS hits found in cluster
            for IS_hit in IS_hits:
                IS_SMILE = G.nodes[IS_hit]['smile']
                E_SMILEs = [G.nodes[E_ID]['smile'] for E_ID in E_hits ]
                SIMs = []
                for E_SMILE in E_SMILEs:
                    try: # some experimental specra have no smiles
                        SIMs.append(cheminfo_tools.dice(IS_SMILE, E_SMILE))
                    except:
                        pass # G.remove_node(IS_hits[0])

                AVE_SIM = np.average(SIMs)
                # if np.isnan(AVE_SIM) or AVE_SIM <= 0.7:
                if AVE_SIM <= 0.7:
                    neighbor = list(G[IS_hit])
                    G.nodes[neighbor[0]]['level'] = 'unmatch'
                    G.remove_node(IS_hit) # If fail to parse edb smile, remove is match

                # elif :  # Delete low similar IS hit
                #     neighbor = list(G[IS_hit])
                #     G.nodes[neighbor[0]]['level'] = 'unmatch'
                #     G.remove_node(IS_hit)

        # Step 2
        if not E_hits and IS_hits:
            if len(IS_hits) == 1:
                G.remove_node(IS_hits[0]) # remove single is match

            elif len(IS_hits) > 1:
                G.remove_node(IS_hits[0])
            else:
                IS_SMILEs = [G.nodes[h]['smile'] for h in IS_hits]
                IDXs = cheminfo_tools.dice_pair(IS_SMILEs)
                nodes_to_remove = [IS_hits[idx] for idx in IDXs]
                for node in nodes_to_remove:
                    # neighbor = list(G[node])
                    # G.nodes[neighbor[0]]['level'] = 'unmatch'
                    G.remove_node(node)

    return G

def molecular_generation(args):
    '''
    Generating the final molecular SpecSimNetwork
    Unmatched
    Experimental matched
    In silico matched
    '''

    basename = os.path.splitext(os.path.basename(args.spectra_file))[0]
    parent_folder = f'{args.output}/{basename}_result'
    # quant_df = functions.df_preprocess(args.quant_file)
    exp_info = functions_new.spectra_process(args.spectra_file)

    G = self_clustering(args,exp_info) # Networking
    feature_ids = [str(key) for key, values in exp_info.items()]

    # Loading and preprocessing the in-silico MS1 match result generated by 'ms1_match()'
    isdbms1_result_path = os.path.join(parent_folder, f"IS_MS1match_{basename}.csv")
    isms1_match_df = functions_new.df_preprocess(isdbms1_result_path)
    isms1_match_df[f'{args.library_matching_method}_similarity'] = np.nan # Initialize an empty column to record the maximum values

    '''In silico —— load'''
    # 0, 1, 2 represent the in-silico MS2 predicted by CFM-ID at three collision energy (10, 20 and 40 eV)
    for i in range(len(isms1_match_df)):
        SIM_NP = [
            (isms1_match_df.loc[i, f'{args.library_matching_method}_similarity0'],
             isms1_match_df.loc[i, f'{args.library_matching_method}_mps0']),
            (isms1_match_df.loc[i, f'{args.library_matching_method}_similarity1'],
             isms1_match_df.loc[i, f'{args.library_matching_method}_mps1']),
            (isms1_match_df.loc[i, f'{args.library_matching_method}_similarity2'],
             isms1_match_df.loc[i, f'{args.library_matching_method}_mps2'])]


        best_SIM, best_NP = max(SIM_NP, key=lambda x: x[1]) # pick pair with maximum shared peaks
        isms1_match_df.loc[i,f'{args.library_matching_method}_similarity'] = best_SIM # Convert missing values to N/A
        isms1_match_df.loc[i,f'{args.library_matching_method}_mps'] = best_NP # Convert missing values to N/A

    index_match,index_unmatch = [], []
    for feature_id in feature_ids:  # traverse IS_MS1match_result by feature ids
        temp_df = isms1_match_df[isms1_match_df['row ID'] == int(feature_id)] # Temp df with the same feature ID


            # if temp_df is not None and not temp_df.empty: # make sure not empty
        sim_idx = temp_df[f'{args.library_matching_method}_similarity'].idxmax()  # get index of maximum spec sim
        index_match.append(sim_idx) # record idx
    index_match = [x for x in index_match if pd.notna(x)]  # Remove Nan

    '''In silico —— filter'''
    df_new_match = isms1_match_df.loc[index_match].reset_index(drop=True) # filter the dataframe by indeces
    df_new_match_well = df_new_match[(df_new_match[f'{args.library_matching_method}_similarity'] >= args.library_matching_similarity)
                                     & (df_new_match[f'{args.library_matching_method}_mps'] >= args.library_matching_peaks)].reset_index(drop=True)

    for i in range(len(df_new_match_well)):
        pair_sim = df_new_match_well.loc[i, f'{args.library_matching_method}_similarity']
        matched_peaks = int(df_new_match_well.loc[i, f'{args.library_matching_method}_mps'])
        IS_SMILE = df_new_match_well.loc[i, 'match_smiles']
        spec1_id = str(df_new_match_well.loc[i, 'row ID'])
        spec2_id = str(df_new_match_well.loc[i, 'match_id'])

        edge_attr = {f'{args.library_matching_method}_similarity': pair_sim,
                     'matched_peaks': matched_peaks,
                     'edge_type': f'{args.library_matching_method}'}
        G.add_edge(spec1_id, spec2_id, **edge_attr)
        # G.nodes[spec1_id]['level'] = 'ISmatch'
        G.nodes[spec2_id]['class'] = 'ISDB'
        G.nodes[spec2_id]['level'] = 'DB'
        G.nodes[spec2_id]['smile'] = IS_SMILE


    '''Experimental match-Load'''
    # Loading and preprocessing the experimental MS1 match result generated by 'ms1_match()'
    edbms1_result_path = os.path.join(parent_folder, f'E_MS1match_{basename}.csv')
    edbms1_match_df = functions.df_preprocess(edbms1_result_path)

    # edbms1_match_df[f'{args.library_matching_method}_similarity'] \
    #     = pd.to_numeric(edbms1_match_df[f'{args.library_matching_method}_similarity'], errors='coerce')

    edb_index_match, edb_index_unmatch = [], []
    for feature_id in feature_ids:
        temp_df = edbms1_match_df[edbms1_match_df['row ID'] == int(feature_id)]
        idx = temp_df[f'{args.library_matching_method}_mps'].idxmax()  # get index of match with maximum pair_similarity
        edb_index_match.append(idx)
    edb_index_match = [x for x in edb_index_match if pd.notna(x)] # Remove Nan

    #     elif pd.isna(temp_df['match_id']).any():  # get index of features without MS1 match
    #         edb_index_unmatch.extend(temp_df.index.values.tolist())

    '''Experimental match —— filter'''
    edb_df_new_match = edbms1_match_df.loc[edb_index_match].reset_index(drop=True) # Create a filtered new copy
    edb_quant_df_new_match_well = edb_df_new_match[
        (edb_df_new_match[f'{args.library_matching_method}_similarity'] >= args.library_matching_similarity) & (
                edb_df_new_match[f'{args.library_matching_method}_mps'] >= args.library_matching_peaks)].reset_index(drop=True)

    for i in range(len(edb_quant_df_new_match_well)):
        pair_sim = edb_quant_df_new_match_well.loc[i, f'{args.library_matching_method}_similarity']
        matched_peaks = int(edb_quant_df_new_match_well.loc[i, f'{args.library_matching_method}_mps'])
        E_SMILE = edb_quant_df_new_match_well.loc[i, 'match_smiles']
        spec1_id = str(edb_quant_df_new_match_well.loc[i, 'row ID'])
        spec2_id = str(edb_quant_df_new_match_well.loc[i, 'match_id'])
        edge_attr = {f'{args.library_matching_method}_similarity': pair_sim,
                     'matched_peaks': matched_peaks,
                     'edge_type': f'{args.library_matching_method}'}
        G.add_edge(spec1_id, spec2_id, **edge_attr)


        # G.nodes[spec1_id]['level'] = 'Ematch'
        G.nodes[spec2_id]['class'] = 'EDB'
        G.nodes[spec2_id]['level'] = 'DB'
        G.nodes[spec2_id]['smile'] = f'{E_SMILE}'

    '''IS annotation filtering'''
    G = ismatch_filter(G)

    '''Level attirbution'''
    for node in G.nodes():
        LEVEL =  G.nodes[node]['class']
        if LEVEL == 'ISDB':
            neighbors = list(G[node])
            for neighbor in neighbors:
                G.nodes[neighbor]['level'] = 'ISmatch'

        elif LEVEL == 'EDB':
            neighbors = list(G[node])
            for neighbor in neighbors:
                G.nodes[neighbor]['level'] = 'Ematch'

    MN_file = os.path.join(parent_folder,
                           f'{basename}_{args.library_matching_method}_{args.self_clustering_method}_{args.self_clustering_similarity}_{args.self_clustering_peaks}.graphml')
    nx.write_graphml(G, MN_file)

def list_graphml_files(folder_path,keywords=None):
    '''
    return the absolute path of selected graphml files
    '''
    graphml_files = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".graphml") and "—" not in filename:
            if keywords is not None:
                for keyword in keywords:
                    if filename in keyword:
                        graphml_files.append(os.path.join(folder_path, filename))
                        break
    return graphml_files

def mn_merging(mn1_G,mn2_G):
    # mn1_file = mn1_file
    # mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
    # mn1_G = nx.read_graphml(mn1_file)
    node_ids = list(mn1_G.nodes())

    # mn2_file = mn2_file
    # mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
    # mn2_G = nx.read_graphml(mn2_file)

    # node_merging
    for node in set(mn1_G.nodes()) | set(mn2_G.nodes()):
        attrs1 = mn1_G.nodes[node] if node in mn1_G else {}  # get node attributes from two graph
        attrs2 = mn2_G.nodes[node] if node in mn2_G else {}

        combined_attrs = {**attrs1, **attrs2}  # merge attibute to make sure no duplicate
        mn2_G.add_node(node, **combined_attrs)  # add node and attibutes to new graph

    # edge_merging
    for src, dst, edge_attrs in mn1_G.edges(data=True):

        if not mn2_G.has_edge(src, dst):  # check if edge in mn2 exist
            mn2_G.add_edge(src, dst, **edge_attrs)
        else:
            existing_edges = mn2_G.get_edge_data(src, dst)  # get all attributes of this node
            # Assuming mn2_G is a MultiGraph, all possible edge attributes need to be checked
            # For a normal graph, existing_edges is a dictionary, which may have only one key-value pair
            if isinstance(existing_edges, dict) and existing_edges:
                for key, attrs in existing_edges.items():
                    if attrs == edge_attrs:
                        break
                else:
                    mn2_G.add_edge(src, dst, **edge_attrs)  # if no edge attributes, find new attributes
            else:
                mn2_G.add_edge(src, dst, **edge_attrs)  # if edge attributes in mn2_G doesn't match，add new edges

    return mn2_G
    # nx.write_graphml(mn2_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')

def multimerging(gfiles):

    if len(gfiles) < 2:
        print("Need at least two molecular SpecSimNetwork")

    if len(gfiles) == 2:
        mn1_basename = os.path.basename(gfiles[0]).replace('.graphml', '')
        mn2_basename = os.path.basename(gfiles[1]).replace('.graphml', '')
        graphml_name = mn1_basename +'-'+ mn2_basename
        mn1_G = nx.read_graphml(gfiles[0])
        mn2_G = nx.read_graphml(gfiles[1])
        merged_G = mn_merging(mn1_G, mn2_G)
        nx.write_graphml(merged_G, f'{args.output}/{graphml_name}.graphml')

    else:
        mn1_basename = os.path.basename(gfiles[0]).replace('.graphml', '')
        mn2_basename = os.path.basename(gfiles[1]).replace('.graphml', '')
        mn1_G = nx.read_graphml(gfiles[0])
        mn2_G = nx.read_graphml(gfiles[1])
        merged_G = mn_merging(mn1_G, mn2_G)
        graphml_name = mn1_basename + '—' + mn2_basename
        for i in range(2, len(gfiles)):
            mni_G = nx.read_graphml(gfiles[i])
            mni_basename = os.path.basename(gfiles[i]).replace('.graphml', '')
            merged_G = mn_merging(merged_G, mni_G)  # Gradually merge with subsequent files
            graphml_name += '—' + mni_basename
        nx.write_graphml(merged_G, f'{args.output}/{graphml_name}.graphml')

if __name__ == '__main__':
    t = time.time()
    '''Param settings'''
    args = config.args
    # Input
    args.spectra_file = '../msdb/data/kutz/KutzOsmac.mgf'
    args.output =  '../msdb/data/kutz/'

    args.allowed_mass_tolerance = 10
    optimal_threhold = evaluation.get_optimal_threshold()

    args.library_matching_method = 'entropy'
    args.library_matching_similarity = optimal_threhold[args.library_matching_method]
    args.library_matching_peaks = 3

    args.self_clustering_method = 'modified_cosine'
    args.self_clustering_similarity = optimal_threhold[args.self_clustering_method]
    args.self_clustering_peaks = 3
    args.top_k = 5


    args.edbms1_file = '../msdb/edbMS1.csv'
    args.edbms2_file = '../msdb/edb_info.json'

    args.isms1_file = '../msdb/isdbMS1.csv'
    args.isms2_file = '../msdb/isdb_info.json'
    queryMGF = None

    args.merge_list = None

    '''************************************************************************************************************'''
    # ms1_match(args)
    ems2_match(args)
    isms2_match(args)
    molecular_generation(args)

    # G_FILEs = [
    #     # '../msdb/data/kutz/KutzOsmac_result_former/KutzOsmac_entropy_modified_cosine_0.72_3.graphml',
    #            '../msdb/data/kutz/KutzOsmac_result_former/KutzOsmac_modified_cosine_modified_cosine_0.72_3.graphml',
    #            '../msdb/data/kutz/KutzOsmac_result_former/KutzOsmac_peak_percentage_modified_cosine_0.72_3.graphml']
    # multimerging(G_FILEs)

    # exp_info = functions_new.spectra_process(args.spectra_file)
    # G = self_clustering(args,exp_info) # Networking
    # nx.write_graphml(G, f'{args.output}/SC_{args.self_clustering_method}_{args.self_clustering_similarity}_{args.self_clustering_peaks}.graphml')



    print(f'Finish in {(time.time() - t) / 60:.2f}min')

'''
Customized mass database
'''
import os
import json
import pandas as pd
import numpy as np
from my_packages import functions, config,functions_new
from matchms.exporting import save_as_mgf
from matchms import Spectrum

def customized_db(args):
    '''
    args.mgf_file
    args.library_info
    '''
    '''Specify output'''
    out_dir = args.output  # Create output directory or skip if exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        pass
    '''Load files for creating database'''
    mgf_csv = args.library_info  # related csv file
    mgf = args.spectra_file

    SPECTRA = functions_new.load_spectra_from_file(mgf)
    # Columns: ['compound_name', 'feature_ID', 'smiles', 'collector', 'adduct','organism', 'pepmass']
    DF = functions_new.df_preprocess(mgf_csv)

    ms1_df = pd.DataFrame(columns=['id', 'pepmass', 'smiles'])  # To create ms1 file
    ms2_data = {}  # To create ms2 file

    for index, row in DF.iterrows():
        # For ms1_library
        SPEC_INFO = SPECTRA[index]
        pm = SPEC_INFO.metadata['precursor_mz']
        ID = f"cur{row['feature_ID']}"
        SMILE = row['smiles']
        new_row = pd.Series({'id': ID, 'pepmass': pm, 'smiles': SMILE})
        ms1_df = pd.concat([ms1_df, new_row.to_frame().T], ignore_index=True)

        # For ms2 library
        charge = SPEC_INFO.metadata['charge']
        try:
            ionmode = SPEC_INFO.metadata['ionmode']
        except:
            ionmode = ''
        try:
            compound_name = SPEC_INFO.metadata['compound_name']
        except:
            compound_name = ''

        ms2_data[ID] = {
            "pepmass": pm,
            "charge": charge,
            "ion_mode": ionmode,
            "smiles": SMILE,
            "compound_name": compound_name,
            "ms2": np.column_stack((SPEC_INFO.mz, SPEC_INFO.intensities)).tolist()}

    # Save ms1_df
    base_filename = os.path.splitext(os.path.basename(mgf))[0]
    ms1_output_filename = f"{out_dir}/{base_filename}_ms1.csv"
    ms1_df.to_csv(ms1_output_filename, index=False)

    # Save ms2_json
    ms2_output_filename = f"{out_dir}/{base_filename}_ms2.json"
    with open(ms2_output_filename, 'w') as f:
        json.dump(ms2_data, f, indent=4, default=functions.arrary2list)

if __name__ == '__main__':

    args = config.args
    customized_db(args)






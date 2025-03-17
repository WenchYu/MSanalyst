# -*- coding: utf-8 -*-
# @Time :2025/2/4 00:32
# @Auther :Yuwenchao
# @Software : PyCharm
'''
customed DB moudle
'''
import os
import json
import pandas as pd
from my_packages import functions, config

def customized_db(args):
    '''
    args.mgf_file
    args.library_info
    '''
    out_dir = args.output
    if not os.path.exists(out_dir): # check or create customed_db directory
        os.makedirs(out_dir)
    else:
        pass

    '''file upload and preprocess'''
    mgf = args.mgf_file
    feature = args.library_info
    # feature = './customed_db/Xiamenmycins_quant.xlsx'
    # mgf = './customed_db/Xiamenmycins.mgf'
    ids = functions.ex_startswith(mgf, start_txt='FEATURE_ID=')
    pepmasses = functions.ex_startswith(mgf, start_txt='PEPMASS=')
    charges = functions.ex_startswith(mgf, start_txt='CHARGE=')
    ionmodes = ['positive' if charge.endswith('+') else 'negative' for charge in charges]
    ms2s = functions.ex_spectra(mgf, start_txt='MSLEVEL=2', end_txt='END IONS')

    df = functions.df_preprocess(feature)
    smiles = df.smiles.tolist()
    compound_names = df.compound_name.tolist()

    '''customed ms1 library generating'''
    ms1_df = pd.DataFrame({
        'id': ids,
        'pepmass': pepmasses,
        'smiles': smiles
    })
    base_filename = os.path.splitext(os.path.basename(mgf))[0]
    ms1_output_filename = f"./customed_db/{base_filename}_ms1.csv"
    ms1_df.to_csv(ms1_output_filename, index=False)

    '''customed ms2 library generating'''
    ms2_data = {}
    if len(ids) != len(df): raise ValueError(f"Lengthes do not match")

    for i in range(len(ids)):
        id = ids[i]
        pepmass = pepmasses[i]
        charge = charges[i]
        smile = smiles[i]
        ionmode = ionmodes[i]
        ms2 = ms2s[i]
        compound_name = compound_names[i]

        ms2_data[id] = {
            "pepmass": pepmass,
            "charge": charge,
            "ion_mode": ionmode,
            "smiles": smile,
            "compound_name": compound_name,
            "ms2": ms2
        }

    ms2_output_filename = f"./{out_dir}/{base_filename}_ms2.json"
    with open(ms2_output_filename, 'w') as f:
        json.dump(ms2_data, f, indent=4, default=functions.arrary2list)


if __name__ == '__main__':
    args = config.args
    customized_db(args)





# -*- coding: utf-8 -*-
# @Time :2023/3/25 10:58
# @Auther :Yuwenchao
# @Software : PyCharm
'''
MS2 search against entire MSanalyst MS2 library
'''
import pandas as pd
from my_packages import functions, ms2tools,config


def ms2search(args):
    '''
    args.mgf_file
    '''
    mgf_info = functions.spectra_process(args.query_ms1, args.query_ms2)
    queryms1 = float(mgf_info.loc[0, 'pepmass'])
    dict = {'row ID': '1', 'row m/z': queryms1}
    query_df = pd.DataFrame([dict])
    ms2tools.ms1_match(args, queryDF=query_df)
    ms2tools.ISDB_MS2_match(args, queryMGF=mgf_info)
    ms2tools.EDB_MS2_match(args, queryMGF=mgf_info)

if __name__=='__main__':
    args = config.args
    ms2search(args)








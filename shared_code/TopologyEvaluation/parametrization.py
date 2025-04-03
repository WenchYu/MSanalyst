# -*- coding: utf-8 -*-
# @Time :2025/4/2 23:42
# @Auther :Yuwenchao
# @Software : PyCharm
'''

Algorithms, thresholds 0.1-0.9, mps 1-5
'''
import os
import time
import pandas as pd

if __name__ == '__main__':
    t = time.time()

    root_dir = './AlgTopologyEva/'
    csv_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith('.csv'):
                csv_file_path = os.path.join(dirpath, filename)
                csv_files.append(csv_file_path)

    dfs,methods,groups,sim_ranges, n20_ranges, nacc_ranges, rccc_ranges = [],[],[],[],[],[],[]
    kwd = 'group4'
    for csv_file in csv_files:
        if kwd in csv_file:
            df = pd.read_csv(csv_file,index_col=None)
            df = df[df['NACC']!=0] # removing non-cluster netowrks
            if not df.empty:
                methods.append(csv_file.split('/')[3])
                groups.append(kwd)
                dfs.append(df)
                sim_range = [i.split('_')[-2] for i in df['name'].tolist()]
                n20_range = [round(x, 2) for x in df['N20'].tolist()]
                nacc_range = [round(x, 2) for x in df['NACC'].tolist()]
                rccc_range = [round(x, 2) for x in df['RCCC'].tolist()]

                sim_ranges.append(f'{max(sim_range)}~{min(sim_range)}')
                n20_ranges.append(f'{max(n20_range)}~{min(n20_range)}')
                nacc_ranges.append(f'{max(nacc_range)}~{min(nacc_range)}')
                rccc_ranges.append(f'{max(rccc_range)}~{min(rccc_range)}')

    sum_df = pd.DataFrame({'method':methods,'sim_range':sim_ranges,'n20_range':n20_ranges
                              ,'nacc_range':nacc_ranges,'rccc_range':rccc_ranges,'group':groups})
    sum_df.to_csv(f'{kwd}.csv', index=False)


    # combined_df = pd.concat(dfs, ignore_index=True)
    # combined_df = combined_df[combined_df['NACC'] != 0]
    # combined_df.to_excel(f'{kwd}.xlsx', index=False)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')

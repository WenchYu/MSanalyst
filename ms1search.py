
'''
MS1 search against entire MSanalyst MS1 library
'''
import sys
sys.path.append('./')
import pandas as pd
from my_packages import ms2tools,config

def ms1search(args):
    '''
    args.query_ms1
    '''
    query = float(args.query_ms1)
    dict = {'row ID': '1', 'row m/z': query}
    query_df = pd.DataFrame([dict])
    ms2tools.ms1_match(args, queryDF=query_df)

if __name__ == '__main__':
    args = config.args
    ms1search(args)






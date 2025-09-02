
'''

'''
import sys
sys.path.append('./')
from my_packages import functions

from my_packages import functions_new,ms2tools_new,config,evaluation

def main(args):
    '''Main workflow of MSanalyst'''
    functions_new.create_result_folders(args)
    ms2tools_new.ms1_match(args)
    ms2tools_new.ISDB_MS2_match(args)
    ms2tools_new.ems2_match(args)
    functions.create_subresults(args)
    ms2tools_new.molecular_generation(args)
    evaluation.connectivity_screening(args)

if __name__ == '__main__':
    args = config.args
    # main(args)





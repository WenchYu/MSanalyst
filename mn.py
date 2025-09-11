
'''

'''
import os
import sys
sys.path.append('./')
from my_packages import functions_new,ms2tools_new,config

def main(args):
    '''Main workflow of MSanalyst'''
    functions_new.create_result_folders(args)
    ms2tools_new.ms1_match(args)
    ms2tools_new.isms2_match(args)
    ms2tools_new.ems2_match(args)
    functions_new.create_subresults(args)
    ms2tools_new.molecular_generation(args)

if __name__ == '__main__':
    args = config.args
    if args.input_folder:
        functions_new.preprocess_mzml(args) # Feature ID differs between batches
        args.spectra_file = os.path.join(args.input_folder,'MS2data.mgf')
    else:
        pass
    main(args)





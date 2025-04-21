
'''
Re clustering the generated network
'''
import sys
sys.path.append('./')
from my_packages import ms2tools, config

def re_networking(args):
    MN_file = ms2tools.molecular_generation(args)
    return MN_file
    

if __name__ == '__main__':
    args = config.args
    re_networking(args)
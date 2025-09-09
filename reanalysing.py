
'''
Re clustering the generated SpecSimNetwork
'''
import sys
sys.path.append('./')
from my_packages import ms2tools_new, config

def re_networking(args):
    ms2tools_new.ems2_match(args)
    ms2tools_new.isms2_match(args)
    MN_file = ms2tools_new.molecular_generation(args)
    return MN_file


if __name__ == '__main__':
    args = config.args
    re_networking(args)
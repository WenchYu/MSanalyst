# -*- coding: utf-8 -*-
# @Time : 2025/2/7 14:27
# @Auther : Yuwenchao
# @Software : PyCharm
'''
Re clustering the generated network
'''
import time
import argparse
from my_packages import ms2tools, config

if __name__ == '__main__':
    args = config.args
    # ms2tools.self_clustering(args)
    ms2tools.molecular_generation(args)
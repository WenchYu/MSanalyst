# -*- coding: utf-8 -*-
# @Time :2022/12/11 19:35
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''

import argparse
import time
from my_packages import functions,ms2tools,config

def main(args):
    '''Main workflow of MSanalyst'''
    functions.create_result_folders(args)
    # ms2tools.spectral_entropy_calculating(args)
    ms2tools.ms1_match(args)
    ms2tools.ISDB_MS2_match(args)
    ms2tools.EDB_MS2_match(args)
    functions.create_subresults(args)
    ms2tools.molecular_generation(args)



if __name__ == '__main__':
    args = config.args
    functions.create_result_folders(args)
    # ms2tools.spectral_entropy_calculating(args)
    ms2tools.ms1_match(args)
    ms2tools.ISDB_MS2_match(args)
    ms2tools.EDB_MS2_match(args)
    functions.create_subresults(args)
    ms2tools.molecular_generation(args)

    # parser = argparse.ArgumentParser(
    #     prog="MSanalyst",
    #     description="MSanalyst designed for molecular networking and annotation",
    #     usage="python msanalyst.py main -q xxx_quant.csv -m xxx.mgf -o output_path"
    # )
    # subparsers = parser.add_subparsers(help="sub-command help")
    #
    # '''subcommand : main'''
    # parser_main = subparsers.add_parser("main", help="Default analysis workflow of MSanalyst")
    # parser_main.add_argument("-q", "--quant_file", help="Quantitative table exported by MZmine", default="./example/example_quant.csv")
    # parser_main.add_argument("-m", "--mgf_file", help="Mgf file exported by MZmine", default="./example/example.mgf")
    # parser_main.add_argument("-o", "--output", help="Output path", default="./example/")
    # parser_main.add_argument("-i1f", "--isms1_file", help="in-silico ms1 file", default="./msdb/isdbMS1.csv")
    # parser_main.add_argument("-e1f", "--edbms1_file", help="experimental ms1 file", default="./msdb/edbMS1.csv")
    # parser_main.add_argument("-i2f", "--isms2_file", help="in-silico  library", default="./msdb/isdb_info.json")
    # parser_main.add_argument("-e2f", "--edbms2_file", help="experimental ms2 library", default="./msdb/edb_info.json")
    # parser_main.add_argument('-pmt'
    #                          , '--pepmass_match_tolerance'
    #                          , help = 'Allowed ppm tolerance in MS1 matching'
    #                          , type = int
    #                          , default = 5
    #                          )
    # parser_main.add_argument('-lmm'
    #                          , '--library_matching_method'
    #                          ,help='Similarity algorithm of tandem mass matching used for library search'
    #                          ,default='modified_cosine_similarity'
    #                          )
    # parser_main.add_argument('-scm'
    #                          , '--self_clustering_method'
    #                          , help='Tandem mass self clustering methods'
    #                          ,default='modified_cosine'
    #                          )
    # parser_main.add_argument('-scs'
    #                          , '--self_clustering_similarity'
    #                          , help='Self clustering similarity threshold'
    #                          , type=float
    #                          ,default=0.7
    #                          )
    # parser_main.add_argument('-scp'
    #                          , '--self_clustering_peaks'
    #                          , help='Self clustering shared peaks threshold'
    #                          , type=int
    #                          , default=5
    #                          )
    # parser_main.add_argument('-topk'
    #                     ,'--top_k'
    #                     , help='Maximum degree of a node'
    #                     , type=int
    #                     , default=10
    #                     )
    # parser_main.add_argument('-islms'
    #                          , '--is_library_matching_similarity'
    #                          , help='In silico library matching similarity threshold'
    #                          , type=float
    #                          , default=0.7
    #                          )
    # parser_main.add_argument('-islmp'
    #                          , '--is_library_matching_peaks'
    #                          , help='In silico library matching shared peaks threshold'
    #                          , type=int
    #                          , default=5
    #                          )
    # parser_main.add_argument('-lms'
    #                          , '--library_matching_similarity'
    #                          , help='Library matching similarity threshold'
    #                          ,type=float
    #                          ,default=0.7
    #                          )
    # parser_main.add_argument('-lmp'
    #                          , '--library_matching_peaks'
    #                          , help='Library matching shared peaks threshold'
    #                          , type=int
    #                          ,default=5
    #                          )
    # parser_main.add_argument('-ppt'
    #                        , '--peak_percentage_threshold'
    #                        , help='Library matching shared peaks percentage threshold'
    #                        , type=float
    #                        , default=0.7
    #                        )
    # parser_main.set_defaults(func=main)
    #
    # '''subcommand : mn'''
    # parser_mn = subparsers.add_parser('mn', help='Re-analysis of MSanalyst results')
    # parser_mn.add_argument('-q', '--quant_file'
    #                          , help='Quantitative table exported by MZmine'
    #                          ,required=True
    #                          , default='./example/example_quant.csv'
    #                          )
    # parser_mn.add_argument('-m', '--mgf_file'
    #                          , help='Mgf file exported by MZmine'
    #                          , required=True
    #                          , default='./example/example.mgf'
    #                          )
    # parser_mn.add_argument('-o', '--output'
    #                          , help='Output path'
    #                          , required=True
    #                          , default='./example/'
    #                          )
    # parser_mn.add_argument('-pmt'
    #                          , '--pepmass_match_tolerance'
    #                          , help='Allowed ppm tolerance in MS1 matching'
    #                          , type=int
    #                          ,default=5
    #                          )
    # parser_mn.add_argument('-lmm'
    #                          , '--library_matching_method'
    #                          , help='Similarity algorithm of tandem mass matching used for library search'
    #                          , default='weighted_dot_product'
    #                          )
    # parser_mn.add_argument('-scm'
    #                          , '--self_clustering_method'
    #                          , help='Tandem mass self clustering methods'
    #                          , default='weighted_dot_product'
    #                          )
    # parser_mn.add_argument('-scs'
    #                          , '--self_clustering_similarity'
    #                          , help='Self clustering similarity threshold'
    #                          , type=float
    #                          , default=0.7
    #                          )
    # parser_mn.add_argument('-scp'
    #                          , '--self_clustering_peaks'
    #                          , help='Self clustering shared peaks threshold'
    #                          , type=int
    #                          , default=6
    #                          )
    # parser_mn.add_argument('-topk'
    #                     ,'--top_k'
    #                     , help='Maximum degree of a node'
    #                     , type=int
    #                     , default=10
    #                     )
    # parser_mn.add_argument('-islms'
    #                          , '--is_library_matching_similarity'
    #                          , help='In silico library matching similarity threshold'
    #                          , type=float
    #                          , default=0.7
    #                          )
    # parser_mn.add_argument('-islmp'
    #                          , '--is_library_matching_peaks'
    #                          , help='In silico library matching shared peaks threshold'
    #                          , type=int
    #                          , default=6
    #                          )
    # parser_mn.add_argument('-lms'
    #                          , '--library_matching_similarity'
    #                          , help='Library matching similarity threshold'
    #                          , type=float
    #                          , default=0.7
    #                          )
    # parser_mn.add_argument('-lmp'
    #                          , '--library_matching_peaks'
    #                          , help='Library matching shared peaks threshold'
    #                          , type=int
    #                          , default=6
    #                          )
    # parser_mn.add_argument('-ppt'
    #                        , '--peak_percentage_threshold'
    #                        , help='Library matching shared peaks perventage threshold'
    #                        , type=float
    #                        , default=0.7
    #                        )


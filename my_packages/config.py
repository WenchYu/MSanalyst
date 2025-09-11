
'''
Parameters used in MSanalyst
'''
import argparse
import ast

def arg_parse():
    parser = argparse.ArgumentParser(
        prog="MSanalyst",
        description="MSanalyst designed for molecular networking and annotation",
        usage="python mn.py -m xxx.mgf -o output_path"
    )

    '''In/output and database selecting'''
    parser.add_argument("-m", "--spectra_file", help="Spectra in mgf", default="./example/example.mgf")
    parser.add_argument("-f", "--input_folder", help="Folder containing at least two mzML files. UmetaFlow used for preprocessing ")
    parser.add_argument("-o", "--output", help="Output path", default="./example/")
    parser.add_argument("-e1f", "--edbms1_file", help="experimental ms1 library file", default="./msdb/edbMS1.csv")
    parser.add_argument("-e2f", "--edbms2_file", help="experimental ms2 library file", default="./msdb/edb_info.json")

    parser.add_argument("-i1f", "--isms1_file", help="in-silico ms1 library", default="./msdb/isdbMS1.csv")
    parser.add_argument("-i2f", "--isms2_file", help="in-silico ms2 library", default="./msdb/isdb_info.json")

    '''mzML,mzXML file preprocess'''
    parser.add_argument("-nl", "--noise_level", help="Removing the noise below setting threhold", default=1e4)


    '''General parameters'''
    parser.add_argument('-mt'
                             , '--allowed_mass_tolerance'
                             , help='Allowed ppm tolerance in MS1 matching'
                             , type=float
                             , default = 10.0
                             )

    '''Library searching parameters'''

    parser.add_argument('-lmm'
                             , '--library_matching_method'
                             , help='Similarity algorithm of tandem mass matching used for library search'
                             , default='modified_cosine'
                             )
    parser.add_argument('-lms'
                        , '--library_matching_similarity'
                        , help='Library matching similarity threshold'
                        , type=float
                        , default=0.7
                        )
    parser.add_argument('-lmp'
                        , '--library_matching_peaks'
                        , help='Library matching shared peaks threshold'
                        , type=int
                        , default=0
                        )



    parser.add_argument('-islms'
                        , '--is_library_matching_similarity'
                        , help='In silico library matching similarity threshold'
                        , type=float
                        , default=0.7
                        )
    parser.add_argument('-islmp'
                        , '--is_library_matching_peaks'
                        , help='In silico library matching shared peaks threshold'
                        , type=int
                        , default=5
                        )

    parser.add_argument('-ppt'
                        , '--peak_percentage_threshold'
                        , help='Library matching shared peaks percentage threshold'
                        , type=float
                        , default=0.7
                        )

    '''Self-clustering parameters'''
    parser.add_argument('-scm'
                             , '--self_clustering_method'
                             , help='Tandem mass self clustering methods'
                             , default='modified_cosine'
                             )
    parser.add_argument('-scs'
                             , '--self_clustering_similarity'
                             , help='Self clustering similarity threshold'
                             , type=float
                             , default=0.7
                             )
    parser.add_argument('-scp'
                             , '--self_clustering_peaks'
                             , help='Self clustering shared peaks threshold'
                             , type=int
                             , default=0
                             )
    parser.add_argument('-topk'
                             , '--top_k'
                             , help='Maximum degree of a node'
                             , type=int
                             , default= 100
                             )

    parser.add_argument("-qms1", "--query_ms1", help="MS1 search against entire MSanalyst library",
                        default="")
    parser.add_argument("-qms2", "--query_ms2", help="MS2 search against entire MSanalyst library",
                        default="")
    parser.add_argument("-sc", "--spectrum_clean", help="MS1 search against entire MSanalyst library",
                        type = bool , default=True)
    parser.add_argument('-li',"--library_info", type=str, help="csv file for genertating containing standard info")
    parser.add_argument('-mn1', "--mn1_file", type=str, help="Molecular SpecSimNetwork file 1")
    parser.add_argument('-mn2', "--mn2_file", type=str, help="Molecular SpecSimNetwork file 2")
    parser.add_argument('-ml', "--merge_list", type=ast.literal_eval, help="absolute path of graphml files to merge")
    parser.add_argument('-cpu', "--cpus", type=int, help="The number of CPUs allowed to be used",default=32)
    parser.add_argument('-gml', "--graphml_file", type=str, help="SpecSimNetwork .graphml file")

    return parser.parse_args()

args = arg_parse()
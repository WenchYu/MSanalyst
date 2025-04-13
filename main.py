
'''

'''

from my_packages import functions,ms2tools,config

def main(args):
    '''Main workflow of MSanalyst'''
    functions.create_result_folders(args)
    ms2tools.ms1_match(args)
    ms2tools.ISDB_MS2_match(args)
    ms2tools.EDB_MS2_match(args)
    functions.create_subresults(args)
    ms2tools.molecular_generation(args)


if __name__ == '__main__':
    args = config.args
    main(args)

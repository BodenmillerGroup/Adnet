import pandas as pd
import os
from adnet import library as lib

def convert_2_csv(folder, fn_bindat, fn_completedat, bindat_only=False):
    """
    Convert the adnet output to csv.
    """

    bin_dat = pd.read_pickle(os.path.join(folder, fn_bindat))
    dat_perbin, dat_persample = lib.bindat_2_simple_tables(bin_dat)
    dat_perbin.to_csv(os.path.join(folder, 'bindat_perbin.csv'))
    dat_persample.to_csv(os.path.join(folder, 'bindat_persample.csv'))
    if bindat_only == False:
        complete_dat = pd.read_pickle(os.path.join(fn_completedat))
        complete_dat.to_csv(os.path.join(folder, 'complete_dat.csv'))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Convert Adnet output pickle files to csv')
    parser.add_argument('folder', type=str,
                        help='The output folder of the adnet analysis')
    parser.add_argument('--bindat_name', type=str, default='bindat',
                       help='The name of the bin_dat output file')
    parser.add_argument('--completedat_name', type=str, default='complete_dat',
                       help='The name of the complete_dat output file')
    parser.add_argument('--bindat_only', action='store_true', default=False)

    # parse the arguments
    args = parser.parse_args()
    convert_2_csv(args.folder, args.bindat_name, args.completedat_name, args.bindat_only)

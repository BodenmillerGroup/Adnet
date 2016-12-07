#!/usr/bin/python
from __future__ import division
import argparse
import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import library as lib
import ast
import scipy.stats as scstats
from configparser2 import ConfigParser
import itertools


from adnet import adnet_steps

__author__ = 'vitoz'

"""
Setup
"""

# TODO: raise errors if metadata not the same for all the plates


def adnet_analysis(parser):
    """
    Perform the whole bp-R2 analysis.

    :param parser:
    :return:
    """

    """
    Begin script
    """

    base_folder = parser.get('folders', 'base_folder')
    sub_folders = ast.literal_eval(parser.get('folders', 'sub_folders'))
    plot_folder = parser.get('folders', 'plot_folder')
    name_dict = parser.get('folders', 'name_dict')

    sep = parser.get('filenames','sep')
    row_col_fn_pos = parser.getint('filenames','row_col_fn_pos')

    channel_target = parser.get('channels','channel_target')
    crap_markers = ast.literal_eval(parser.get('channels','crap_markers'))
    neg_ctrl_names = ast.literal_eval(parser.get('channels','neg_ctrl_names'))
    meta_cols_all = ast.literal_eval(parser.get('channels','meta_cols_all'))
    meta_cols_bin = ast.literal_eval(parser.get('channels','meta_cols_bin'))

    nbins = parser.getint('bpr2_options','nbins')
    bin_range = ast.literal_eval(parser.get('bpr2_options','bin_range'))
    min_cells = parser.getint('bpr2_options','min_cells')
    piecew_linear = parser.getboolean('bpr2_options', 'piecew_linear')
    bin_stat = parser.get('bpr2_options','bin_stat')

    cutoff_stat = parser.get('cuttof','cutoff_stat')
    cutoff_method = parser.get('cuttof','cutoff_method')
    cutoff_value = parser.getfloat('cuttof','cutoff_value')
    # how should the cutoff value be determined
    # 'empiric_FDR', 'gamma_FDR', 'gamma_p', 'provided'


    do_plot = parser.getboolean('plotting','do_plot')
    plot_everything = parser.getboolean('plotting','plot_everything')
    plot_pdf = parser.getboolean('plotting','plot_pdf')

    if not(os.path.exists(plot_folder)):
        os.makedirs(plot_folder)

    # Record the configuration file used
    with open(os.path.join(plot_folder, 'config.ini'), 'w') as f:
        parser.write(f)

    """
    Begin script
    """
    print('Start loading')
    complete_data, channels = adnet_steps._load_all_data(base_folder, sub_folders, row_col_fn_pos, sep)

    print('Clean the single cell data')
    complete_data, channels = adnet_steps._clean_data(complete_data, meta_cols_all, channels, crap_markers, name_dict)


    print('channels used')

    print(channels)
    print('transform data')
    complete_data.loc[:, list(channels)] = lib.transform_arcsinh(complete_data.loc[:, list(channels)])

    print('transformation finished')


    print('start bin calculations')

    bin_dat = adnet_steps._calculate_bin_stats(complete_data, meta_cols_all, meta_cols_bin, channel_target, min_cells,
                                   bin_stat, bin_range, nbins)

    print('pearson/spearman calculation started')
    bin_dat = adnet_steps._calculate_correlations(bin_dat, complete_data)

    print('pearson/spearman calculation finished')



    """
    Calculate the significance cutoff
    """

    print('find cutoff')

    bin_dat, vr_tresh  = adnet_steps._find_cutoff_from_table(bin_dat, cutoff_stat, cutoff_method, neg_ctrl_names, cutoff_value, do_plot, plot_folder)

    print('plot violins')
    if plot_everything:
        idx = bin_dat.index
    else:
        idx = bin_dat[bin_dat['bin_dat_sigfil']].index

    stats = [('stats', 'mean_var_ratio'),
             ('stats', 'corr_pearson_overall'),
             ('stats', 'corr_spearman_overall'),
             ('stats', 'corr_spearman_bin')]


    """
    save the data
    """

    # save bin data
    unibins = bin_dat[bin_stat].columns

    # get the median replicate
    bin_dat[('stats','is_median_varratio')] =bin_dat[('stats', 'mean_var_ratio')]==bin_dat[('stats', 'median_mean_var_ratio')]

    # filter to only keep the median replicate
    med_dat = bin_dat.xs(0, level='timepoint').copy()
    fil = med_dat[('stats','is_median_varratio')]
    med_dat = med_dat[fil]

    #remove all duplicated levels
    med_dat.reset_index(['row_col', 'perturbation', 'origin'], inplace=True)
    med_dat.reset_index('experiment', drop=False,inplace=True)
    med_dat['chosen_exp'] = med_dat['experiment'].groupby(level=['marker','target']).transform(lambda x: x[0])

    med_dat = med_dat[med_dat['chosen_exp'] == med_dat['experiment']]
    idx = med_dat.set_index('experiment', append=True).index
    fn = os.path.join(plot_folder, 'allbinplot_')

    # save bindat
    bin_dat.to_pickle(os.path.join(plot_folder, 'bindat'))
    bin_dat.to_csv(os.path.join(plot_folder, 'bindat.csv'))
    # save completedat
    complete_data.to_pickle(os.path.join(plot_folder, 'complete_dat'))
    complete_data.to_csv(os.path.join(plot_folder, 'complete_dat.csv'))
    print('Plot Violins')

    if plot_everything:
        idx = bin_dat.index
    else:
        idx = bin_dat[bin_dat['bin_dat_sigfil']].index

    stats = [('stats', 'mean_var_ratio'),
             ('stats', 'corr_pearson_overall'),
             ('stats', 'corr_spearman_overall'),
             ('stats', 'corr_spearman_bin')]
    bin_dat.sort_index(inplace=True)
    if do_plot:
        adnet_steps._plot_violins(complete_data, bin_dat, idx, stats, meta_cols_all, bin_range, nbins, min_cells, bin_stat,
                                  plot_pdf, plot_folder,
                                  experiment_str='experiment', marker_str='marker',
                      target_str='target', origin_str='origin', timepoint_str='timepoint')

    print('violins finished')


    print('plot trends')
    ## go through all the overexpressions and plot them versus GFP
    if plot_everything:
        idx = bin_dat.index
    else:
        idx = bin_dat[bin_dat['bin_dat_sigfil']].index

    marker_idx = idx.names.index('marker')
    target_idx = idx.names.index('target')
    origin_idx = idx.names.index('origin')
    marker_target = set((i[marker_idx], i[target_idx]) for i in idx)


    if do_plot:
        adnet_steps._plot_trends()

    print('trends finished')


if __name__ == "__main__":

    conf_parser = ConfigParser()
    # Setup the command line arguments
    parser = argparse.ArgumentParser(description='Run the adnet analysis')
    parser.add_argument('conf_filename', metavar='configuration_file', type=str,
                        help='The filename of the configuration file for the script. Look at the example to see how it '+
                             'should be formatted.')

    # parse the arguments
    args = parser.parse_args()

    # parse the configuration file
    with open(args.conf_filename, 'r') as f:
        conf_parser.read_file(f)

    # run the analysis
    adnet_analysis(conf_parser)

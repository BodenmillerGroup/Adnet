#!/usr/bin/python
from __future__ import division
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import library as lib
import scipy.stats as scstats

import itertools
__author__ = 'vitoz'


def load_all_data(base_folder, sub_folders, row_col_fn_pos, sep):
    """
    Loads the structured data into the data structure

    :param base_folder: name of the base folder
    :param sub_folders: name of the sub folder
    :param row_col_fn_pos: position of the row column (e.g. A01) in the filename
    :param sep: filename seperator
    :return: complete_data, channels: the dataframe containing all data and metadata, the names of the measured channels
    """

    complete_data, channels = lib.load_xkl_data(base_folder, sub_folders, row_col_fn_pos, sep)

    return complete_data, channels


def clean_data(complete_data, meta_columns, channels, crap_markers, name_dict):
    """

    :param complete_data:
    :param meta_columns:
    :param channels:
    :param crap_markers:
    :param name_dict:
    :return: complete_data, channels
    """
    # Remove the not wanted markers
    channels = set(channels)

    print(channels)
    channels = list(channels.difference(crap_markers))


    # clean up the names of the data table and the channel vector
    try:
        name_dict = pd.read_csv(name_dict)
        name_dict = {row['old']: row['new'] for idx, row in name_dict.iterrows()}
        complete_data.rename(columns=name_dict, inplace=True)
        for i, chan in enumerate(channels):
            if chan in name_dict.keys():
                channels[i] = name_dict[chan]
    except IOError:
        print('file %s not found!' % name_dict)

    # set the metadata columns
    complete_data = complete_data.reset_index(drop=False)
    complete_data = complete_data.set_index(meta_columns, drop=False)
    channels = set(channels)
    channels = list(channels.difference(crap_markers))
    # select the channels
    complete_data = complete_data[list(channels)]
    complete_data = complete_data.sortlevel()
    return complete_data, channels


def calculate_bin_stats(complete_data, meta_cols_all, meta_cols_bin, channel_target, min_cells,
                        bin_stat, bin_range, nbins):
    """
    Calculate the summary statistics over the bins.

    :param complete_data:
    :param meta_cols_all:
    :param meta_cols_bin:
    :param channel_target:
    :param min_cells:
    :param bin_stat:
    :return:
    """
    # group by file
    grouped = complete_data.groupby(level=meta_cols_all)

    # Choose the readout -> if no readout specified, take the one indicated in the 'readout' metadata column
    if channel_target == 'None':
        readout_idx_pos = [i for i, c in enumerate(meta_cols_all) if c == 'readout'][0]
        get_readout_channel = lambda idx: idx[readout_idx_pos]
    else:
        get_readout_channel = lambda idx: channel_target

    bin_dat = list()

    # calculate bins per overexpression and experiment
    bin_min = complete_data.groupby(level=meta_cols_bin).aggregate(lambda x: np.percentile(x, bin_range[0]))
    bin_max = complete_data.groupby(level=meta_cols_bin).aggregate(lambda x: np.percentile(x, bin_range[1]))

    # metacolumns bin
    all_idx2bin_idx = lib.generate_full2reduced(meta_cols_all, meta_cols_bin)

    # iterate over the files
    for idx, group in grouped:

        # if there is no channel target, there must be a readout column in the index
        for chan in complete_data.columns:
            if (channel_target == 'all') | (chan == get_readout_channel(idx)):

                # generate the bins by using the limits calculated over all comparable samples

                bin_idx = all_idx2bin_idx(idx)
                bins = lib.get_n_minmax_bins(bin_min.loc[bin_idx, chan],
                                             bin_max.loc[bin_idx, chan],
                                             nbins)

                # Calculate the bin statistics
                x = group[chan]
                out_dat = lib.get_binned_approx(group, x, bins=bins, minCells=min_cells, stat=bin_stat)
                out_dat = out_dat.reset_index(drop=False)
                out_dat['origin'] = chan

                for name, val in zip(meta_cols_all, idx):
                    out_dat[name] = val

                bin_dat.append(out_dat)

    # aggregate the data
    bin_dat = pd.concat(bin_dat, axis=0)

    print('finish bin calculations')

    # calculate the variance ratios per bin
    bin_dat['var_ratio'] = 1 - (bin_dat['fit_var'] / bin_dat['overall_var'])

    # fix the index
    bin_dat.set_index(meta_cols_all + ['origin', 'target', 'bin'], inplace=True)

    bin_dat = bin_dat.unstack('bin')

    # remove origin == target
    bin_dat = bin_dat[bin_dat.index.get_level_values('origin') != bin_dat.index.get_level_values('target')]
    # calculate bp R2 = mean var med ratio
    bin_dat.loc[:, ('stats', 'mean_var_ratio')] = bin_dat.loc[:, 'var_ratio'].apply(np.nanmean, axis=1)

    # calculate the variance ratio over all replicates: mean and median
    bin_dat[('stats', 'mean_mean_var_ratio')] = bin_dat[('stats', 'mean_var_ratio')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(lambda x: np.nanmean(x))
    bin_dat[('stats', 'median_mean_var_ratio')] = bin_dat[('stats', 'mean_var_ratio')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(lambda x: np.nanmedian(x))
    bin_dat[('stats', 'nr_rep')] = bin_dat[('stats', 'mean_var_ratio')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(lambda x: len(x))

    level_nr = np.array(bin_dat[bin_stat].columns)
    bin_dat[('stats', 'corr_pearson_bin')] = bin_dat[bin_stat].apply(lambda x: scstats.stats.pearsonr(level_nr, x)[0],axis=1)
    bin_dat[('stats', 'corr_spearman_bin')] = bin_dat[bin_stat].apply(lambda x: scstats.stats.spearmanr(level_nr, x)[0],axis=1)


    return bin_dat


def calculate_correlations(bin_dat, complete_data):
    # calculate correlations stats
    bin_dat[('stats', 'corr_pearson_overall')] = np.nan
    bin_dat[('stats', 'corr_spearman_overall')] = np.nan
    complete_data.sort_index(inplace=True)
    bin_dat.sort_index(inplace=True)

    meta_cols = complete_data.index.names
    for idx, grp in complete_data.groupby(level=range(len(meta_cols))):
        for idx_or, grp_or in bin_dat.xs(idx, level=meta_cols).groupby(level=['origin']):
            origin_cells = grp[idx_or]
            for idx_tar, grp_tar in grp_or.groupby(level=['target']):
                target_cells = grp[idx_tar]
                pear = scstats.stats.pearsonr(origin_cells, target_cells)[0]
                spear = scstats.stats.spearmanr(origin_cells, target_cells)[0]
                bin_dat.loc[tuple(list(idx) + [idx_or] + [idx_tar]),
                            [('stats', 'corr_pearson_overall'), ('stats', 'corr_spearman_overall')]] = pear, spear

    bin_dat[('stats', 'median_corr_pearson_overall')] = bin_dat[('stats', 'corr_pearson_overall')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(lambda x: np.nanmedian(x))
    bin_dat[('stats', 'median_corr_spearman_overall')] = bin_dat[('stats', 'corr_spearman_overall')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(lambda x: np.nanmedian(x))
    bin_dat[('stats', 'median_abs_corr_pearson_overall')] = bin_dat[('stats', 'corr_pearson_overall')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(
        lambda x: np.abs(np.nanmedian(x)))
    bin_dat[('stats', 'median_abs_corr_spearman_overall')] = bin_dat[('stats', 'corr_spearman_overall')].groupby(
        level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).transform(
        lambda x: np.abs(np.nanmedian(x)))

    return bin_dat


def _find_cutoff(values, neg_values, cutoff_method, cutoff_value, do_plot, plot_folder):
    """

    :param values:
    :param neg_values:
    :param cutoff_method:
    :param cutoff_value:
    :return:
    """

    # find the cut-off for 0% false positive rate
    if cutoff_method == 'empiric_FDR':
        vr_tresh = lib.find_fdr_cutoff(values, neg_values, fdr=cutoff_value)

    elif cutoff_method.startswith('gamma_'):
        fit_alpha, fit_loc, fit_beta = scstats.gamma.fit(neg_values)
        pred_neg = scstats.gamma.rvs(fit_alpha, loc=fit_loc, scale=fit_beta, size=1000000)
        p = plt.figure()
        plt.hist(pred_neg, normed=True, range=(0.01, 1), bins=3000, alpha=1, label="sim_ctrl")
        plt.hist(neg_values, normed=True, range=(0.01, 1), bins=300, alpha=0.8, label="ctrl")
        plt.legend(loc=1)
        if do_plot:
            p.savefig(os.path.join(plot_folder, 'Cutoff_gamma' + '.pdf'))

        if cutoff_method == 'gamma_FDR':
            vr_tresh = lib.find_fdr_cutoff(values, neg_values, fdr=cutoff_value)
        elif cutoff_method == 'gamma_p':
            vr_tresh = scstats.gamma.ppf(cutoff_value, fit_alpha, loc=fit_loc, scale=fit_beta)

            # # make the p-vlaue histograms
            # plt.figure()
            # plt.hist(1-stats.gamma.cdf(neg_values, fit_alpha, loc=fit_loc, scale=fit_beta), normed=True)
            # plt.hist(1-stats.gamma.cdf(values, fit_alpha, loc=fit_loc, scale=fit_beta), normed=True)
            #
            # # find the cut-off for 0% false positive rate
            # vr_tresh = -(lib.find_fdr_cutoff(-values, -neg_values, fdr=0.01))
            # is_tn = neg_mark_fil
    elif cutoff_method == 'provided':
        vr_tresh = cutoff_value

    else:
        raise('%s is not a valid cutoff method!' % cutoff_method)

    p = plt.figure()
    plt.hist(values, normed=True, range=(0.01, 1), bins=100, alpha=1, label="Overexpressions")
    plt.hist(neg_values, normed=True, range=(0.01, 1), bins=100, alpha=0.8, label="Control overexpressions")
    # FDR rate
    fdr = (max(sum(neg_values > vr_tresh), 1) / len(neg_values) * len(values)) / sum(values > vr_tresh)
    p.axes[0].axvline(vr_tresh, color='k', linestyle='--',
                      label='Cutoff: %.2f\n Empiric FDR: %.3f' % (vr_tresh, fdr))
    plt.legend(loc=1)
    if do_plot:
        p.savefig(os.path.join(plot_folder, 'Cutoff' + '.pdf'))

    return vr_tresh


def find_cutoff_from_table(bin_dat, cutoff_stat, cutoff_method, neg_ctrl_names, cutoff_value, do_plot, plot_folder):
    """

    :param bin_dat:
    :param cutoff_stat:
    :param cutoff_method:
    :param cutoff_value:
    :return:
    """
    neg_mark_fil = bin_dat.index.get_level_values('marker').isin(neg_ctrl_names)

    if cutoff_stat in ['median_mean_var_ratio', 'mean_mean_var_ratio',
                       'median_corr_spearman_overall',
                       'median_corr_pearson_bin_overall',
                       'median_abs_corr_spearman_overall',
                       'median_abs_corr_pearson_bin_overall']:
        # take the mean as the values should be identical over the replicates
        neg_values = bin_dat[neg_mark_fil][('stats', cutoff_stat)].groupby(
            level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).mean()
        values = bin_dat[neg_mark_fil == False][('stats', cutoff_stat)].groupby(
            level=[u'marker', u'perturbation', u'timepoint', u'origin', u'target']).mean()

    else:
        neg_values = bin_dat[neg_mark_fil][('stats', cutoff_stat)].copy()
        values = bin_dat[neg_mark_fil == False][('stats', cutoff_stat)].copy()

    if cutoff_method != 'provided':
        vr_tresh = _find_cutoff(values, neg_values, cutoff_method, cutoff_value, do_plot, plot_folder)

    else:
        vr_tresh = cutoff_value

    # filter for the median replicate
    bin_dat_sigfil = (bin_dat[('stats', cutoff_stat)] > vr_tresh)
    bin_dat['bin_dat_sigfil'] = bin_dat_sigfil
    bin_dat['bin_dat_sigfil_any'] = bin_dat['bin_dat_sigfil'].groupby(level=['marker', 'experiment', 'target']).transform(np.any)
    bin_dat['bin_dat_sigfil_any_rep'] = bin_dat['bin_dat_sigfil'].groupby(level=['marker', 'target']).transform(np.any)

    return bin_dat, vr_tresh


def plot_violins(complete_data, bin_dat, plot_idx, stats, meta_cols_all, bin_range, nbins, min_cells, bin_stat,
                  plot_pdf, plot_folder,  experiment_str='experiment', marker_str='marker',
                  target_str='target', origin_str='origin', timepoint_str='timepoint'):


    """

    :param complete_data:
    :param bin_dat:
    :param plot_idx:
    :param marker_str:
    :param target_str:
    :param origin_str:
    :return:
    """

    # calculate bins per overexpression and experiment
    marker_idx = plot_idx.names.index(marker_str)
    target_idx = plot_idx.names.index(target_str)
    origin_idx = plot_idx.names.index(origin_str)
    marker_target_origin = set((i[marker_idx], i[target_idx], i[origin_idx]) for i in plot_idx)
    bin_dat.sort_index(inplace=True)

    sns.set_style("whitegrid")
    for marker, target, origin in marker_target_origin:
        plot_dat = complete_data.xs(marker, level=marker_str)[[origin, target]].copy()
        plot_dat_meta = bin_dat.xs((marker, target, origin), level=[marker_str, target_str, origin_str])[stats]
        plot_dat_meta.name = 'mean_var_ratio'
        idx_names = plot_dat_meta.index.names

        # remove entries in plot_dat that are not in the meta data
        plot_dat = plot_dat.loc[plot_dat_meta.index]
        bin_min = plot_dat.groupby(level=['experiment', 'perturbation']).aggregate(lambda x: np.percentile(x, bin_range[0]))
        bin_max = plot_dat.groupby(level=['experiment', 'perturbation']).aggregate(lambda x: np.percentile(x, bin_range[1]))

        plot_dat_meta = plot_dat_meta.reset_index()
        plot_dat_meta['timepoint'] = plot_dat_meta['timepoint'].apply(int)

        exp_par = [c for c in meta_cols_all if c not in ['timepoint', 'marker', 'target', 'origin', 'row_col']]
        plot_dat_meta = plot_dat_meta.set_index(exp_par + ['timepoint'], drop=True)
        plot_dat.reset_index(inplace=True)
        experiment = plot_dat[exp_par].apply(lambda x: tuple(x), axis=1).unique()
        experiment.sort()

        timepoint = plot_dat['timepoint'].unique()
        timepoint.sort()

        plot_dat = plot_dat.set_index(exp_par + ['timepoint'], drop=False)
        plot_dat.sort_index(inplace=True)

        f, sub_plt = plt.subplots(nrows=len(experiment), ncols=len(timepoint),
                                  sharex=False, sharey=True, figsize=(len(timepoint) * 5 + 2, len(experiment) * 5 + 2))
        if len(experiment) * len(timepoint) > 1:
            sub_plt = sub_plt.flatten()
        else:
            sub_plt = [sub_plt]

        for i, idx in enumerate(itertools.product(experiment, timepoint)):
            idx = tuple([c for c in idx[0]] + [idx[1]])
            if idx in plot_dat.index:
                p = sub_plt[i]

                x = plot_dat.loc[idx][origin].copy()
                y = plot_dat.loc[idx][target].copy()

                # Calculate the point density
                # p.scatter(x, y, marker='.', alpha=0.2)

                # get back the original binning

                # cmap = sns.light_palette('blue', as_cmap=True)
                # p.hexbin(x, y, gridsize=30, bins='log', cmap=cmap)
                # p.scatter(x, y, c='blue', marker='.', edgecolor='', cmap=plt.get_cmap('jet'), alpha=0.2)

                idx_dict = {idxname: plot_dat.loc[idx, idxname].iloc[0] for idxname in exp_par}
                idx_dict['marker'] = marker
                # idx_dict['perturbation'] = perturbation
                idx2 = tuple([idx_dict[n] for n in ['experiment', 'perturbation']])
                cur_binmax = bin_max.loc[idx2][origin]
                x[x > cur_binmax] = cur_binmax
                cur_binmin = bin_min.loc[idx2][origin]
                x[x < cur_binmin] = cur_binmin

                bins = lib.get_n_minmax_bins(cur_binmin, cur_binmax, nbins)

                lib.violin_slice(x, y, bins=bins, alpha=0.3, ax=p, color="red", min_cells=min_cells)
                plt.plot(lib.get_bin_stat(x, x, fkt=bin_stat, bins=bins),
                         lib.get_bin_stat(y, x, fkt=bin_stat, bins=bins), color='black')
                p.scatter(x, y, c='blue', marker='.', edgecolor='', cmap=plt.get_cmap('jet'), alpha=0.2)
                ratio = plot_dat_meta.loc[idx][stats[0]]
                pears = plot_dat_meta.loc[idx][stats[1]]
                spear = plot_dat_meta.loc[idx][stats[2]]
                binspear = plot_dat_meta.loc[idx][stats[3]]
                title = "BP-R2: %1.2f\nBin Spearman: %1.2f\nPearson R: %1.2f, Spearman R: %1.2f" % (
                    ratio, binspear, pears, spear)
                p.set_title(title)
                plt.ylabel(idx[:-1])
                plt.xlabel(idx[-1])

                # fix the ylim
        for i in range(0, len(sub_plt), len(timepoint)):
            exp = experiment[int(i / len(timepoint))]
            y = plot_dat.xs(exp, level=exp_par)[target].copy()
            plt.sca(sub_plt[i])
            plt.ylim(0, np.percentile(y, 99) * 1.1)
        plt.subplots_adjust(top=0.85)
        plt.suptitle('Overexpression: ' + marker + ' /Origin: ' + origin + ' /Measured: ' + target)
        f.text(0.5, 0.04, origin + ' - ' + marker + ' [asinh(counts/5)]\n\nEGF timecourse', ha='center')
        f.text(0.04, 0.5, 'Replicates\n\n' + target + ' [asinh(counts/5)]', va='center', rotation='vertical')
        # sub_plt[np.unravel_index(idx[0], dim)].plot(d.index, d, alpha=0.3)

        f.savefig(os.path.join(plot_folder, lib.format_filename('EGF_overexpression_' + marker + '_' + origin +
                                                                '_' + target + '.png')))
        if plot_pdf:
            f.savefig(os.path.join(plot_folder, lib.format_filename('EGF_overexpression_' + marker + '_' + origin +
                                                                    '_' + target + '.pdf')))
        f.clear()
        plt.close()

def plot_trends(marker_target, bin_dat, bin_stat, nbins, plot_pdf, plot_folder):

    uni_mark = set(m for m,t in marker_target)
    plt.ioff()
    for mark in uni_mark:
        uni_targ = [t for m, t in marker_target if m == mark]
        plot_dat = bin_dat.xs(mark, level='marker').copy()
        ix = plot_dat.index.get_level_values('target').isin(uni_targ)
        plot_dat = plot_dat[ix]
        plot_dat = plot_dat[[bin_stat, ]].stack()
        plot_dat = plot_dat.reset_index(level=['target', 'experiment', 'timepoint', 'bin'])
        plot_dat['timepoint'] = plot_dat['timepoint'].apply(int)
        nTP = len(plot_dat['timepoint'].unique())
        rescaled_name = 'rescaled bin ' + bin_stat
        summary_name = bin_stat + ' [asinh(counts/5)]'
        plot_dat[summary_name] = plot_dat[bin_stat]
        plot_dat[rescaled_name] = plot_dat.groupby(['target', 'experiment'])[bin_stat].transform(
            lambda x: x / np.max(x))

        plot_dat.sort_values(by=['timepoint', 'bin'], inplace=True)
        plot_dat['bin'] = plot_dat['bin'] + 1

        p = sns.FacetGrid(plot_dat, row='experiment', col='timepoint', hue='target',
                          margin_titles=True, sharex=False, sharey=True, xlim=(0.5, nbins + 0.5), ylim=(-0.05, 1.05))
        p.map(plt.scatter, 'bin', rescaled_name)
        p.map(plt.plot, 'bin', rescaled_name)
        p.add_legend()
        p.fig.suptitle('Overexpression: ' + mark)
        p.fig.subplots_adjust(top=.9)
        p.savefig(os.path.join(plot_folder, lib.format_filename('Trends_EGF_overexpression_' + mark + '.png')),
                  figsize=(25, 15))
        if plot_pdf:
            p.savefig(os.path.join(plot_folder, lib.format_filename('Trends_EGF_overexpression_' + mark + '.pdf')),
                      figsize=(25, 15))

        plt.close()

        plot_dat.sort_values(by=['timepoint', 'bin'], inplace=True)
        p = sns.FacetGrid(plot_dat, row='experiment', col='target', hue='timepoint',
                          margin_titles=True, sharex=False, sharey=True,
                          palette=sns.cubehelix_palette(nTP + 1)[1:nTP + 1],
                          xlim=(0.5, nbins + 0.5), ylim=(-.05, 1.05 * plot_dat[bin_stat].max()))
        p.map(plt.plot, 'bin', summary_name)
        p.map(plt.scatter, 'bin', summary_name)
        p.add_legend()
        plt.suptitle('Overexpression: ' + mark)
        p.fig.subplots_adjust(top=.9)
        p.savefig(os.path.join(plot_folder,
                               lib.format_filename('Trends_EGF_overexpression_' + mark + '_TPoverBins' + '.png')),
                  figsize=(25, 15))
        if plot_pdf:
            p.savefig(os.path.join(plot_folder,
                                   lib.format_filename('Trends_EGF_overexpression_' + mark + '_TPoverBins' + '.pdf')),
                      figsize=(25, 15))

        plt.close()

        if nTP > 1:
            plot_dat.sort_values(by=['bin', 'timepoint'], inplace=True)
            p = sns.FacetGrid(plot_dat, row='experiment', col='target', hue='bin',
                              margin_titles=True, sharex=False, sharey=True,
                              palette=sns.cubehelix_palette(nbins + 1)[1:nbins + 1],
                              xlim=(np.min(plot_dat['timepoint']) - 0.05, np.max(plot_dat['timepoint']) + 0.05),
                              ylim=(-0.05, 1.05 * plot_dat[bin_stat].max()))
            p.map(plt.plot, 'timepoint', summary_name)
            p.map(plt.scatter, 'timepoint', summary_name)
            p.add_legend()
            plt.suptitle('Overexpression: ' + mark)
            p.fig.subplots_adjust(top=.9)
            p.savefig(os.path.join(plot_folder,
                                   lib.format_filename('Trends_EGF_overexpression_' + mark + '_BinsoverTP' + '.png')),
                      figsize=(25, 15))
            if plot_pdf:
                p.savefig(os.path.join(plot_folder, lib.format_filename(
                    'Trends_EGF_overexpression_' + mark + '_BinsoverTP' + '.pdf')), figsize=(25, 15))
            plt.close()
            # g.map_upper(lib.violin_slice)
            # g.map_upper(lib.plot_spline)
            # g.map_lower(lib.plot_spline)

            # g.map_lower(plt.scatter,marker='.', alpha=0.1)
            # g.map_diag(sns.kdeplot, lw=1)
            # sns.despine()
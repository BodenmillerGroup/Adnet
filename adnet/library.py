from __future__ import division
__author__ = 'vitoz'

import scipy.interpolate as interpolate
import scipy.stats as stats
import os
import fcm
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import string


def load_xkl_data(base_folder, sub_folders, row_col_fn_pos, sep='_'):
    """
    Loads all the fcs files in the base_folder and the sub_folders
    Takes the base_folder name as experiment name  and expects a a file experimentname.csv in each folder containnig
    the metadata for the files
    :param base_folder:
    :param sub_folders:
    :param row_col_fn_pos:
    :return: data, channels
    """
    data_sets = list()
    channels = set()
    get_rc = lambda x: os.path.basename(x).split(sep)[row_col_fn_pos]

    for cur_folder in sub_folders:
        fcs_files = glob.glob(os.path.join(base_folder, cur_folder, 'gated', '*.fcs'))
        dict_fn = glob.glob(os.path.join(base_folder, cur_folder, cur_folder + '*.csv'))[0]
        dict_file = pd.read_csv(dict_fn)
        dict_file['row_col'] = [r+str(c) for r,c in zip(dict_file['row'], dict_file['column'])]
        fcs_files = [f for f in fcs_files if any(rc in get_rc(f) for rc in dict_file['row_col']) ]
        dict_file.set_index('row_col',drop=False,inplace=True)

        fcs_files_id = [get_rc(fn) for fn in fcs_files]
        tmpdata = [fcm.loadFCS(fn, transform=None, auto_comp=False) for fn, id in zip(fcs_files, fcs_files_id)]

        for i,dat in enumerate(tmpdata):
            # make all channel names unique by appending -1,-2 etc
            cur_channels = make_channels_unique(dat.channels)
            channels.update(cur_channels)
            tdat = dat.view()
            tdat = tdat.byteswap().newbyteorder()
            df = pd.DataFrame(tdat,columns=cur_channels)
            df['experiment'] = cur_folder
            df['row_col'] = get_rc(dat.name)
            df.set_index('row_col', inplace=True)
            df = df.join(dict_file)
            tmpdata[i] = df
        data_sets.append(pd.concat(tmpdata))

    data = pd.concat(data_sets,axis=0, ignore_index=True)
    return data, channels



def transform_arcsinh(x, reverse=False, cofactor=5):
    if not reverse:
        return np.arcsinh(x/cofactor)
    else:
        return np.sinh(x)*cofactor

def generate_full2reduced(full_idx, reduced_idx):
    """
    Generates a function that returns a reduced index from a full index
    :param full_idx: vector of column names
    :param reduced_idx: vector of column names, must be a subset of full_idx, however can be ordered differently
    :return: a function that can convert an full index to the reduced form
    """

    assert set(reduced_idx) < set(full_idx)

    pos_reduced_in_full = [i for c_bin in reduced_idx for i, c_all in enumerate(full_idx) if c_all == c_bin]

    def full2reduced(full_idx):
        """
        Returns a reduced index from a full index
        :param full_idx: an iterable corresponding to a full length index
        :return: a tuple corresponding to the reduced index
        """

        return tuple(full_idx[i] for i in pos_reduced_in_full)

    return full2reduced

def censordat(x,perc=0.99):
    xmax = np.percentile(x, perc*100)
    x = map(lambda y: min(y,xmax), x)
    return x

def censor_0(x, minx=1):
    x = x.apply(lambda x: max(x,minx))
    return x

# plots for pairplot
def hexbin(x, y, color, gridsize=50, bins='log', **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=gridsize, bins=bins, cmap=cmap, **kwargs)

def violin_slice(x, y,  nbins=12, ax=None, alpha=0.2,bin_transform = None, bins=None, color=None, min_cells=1):
    if hasattr(x, 'values'):
        x = pd.Series(x.values)
    if hasattr(y, 'values'):
        y = pd.Series(y.values)
    xmax = x.max()
    if bin_transform is None:
        bin_transform = lambda x: x
    if bins is None:
        bins = np.linspace(0, bin_transform(xmax), nbins+1)[1:-1]
    groups = np.digitize(bin_transform(x), bins)
    dat = pd.DataFrame(data={'x':np.array(x), 'groups':groups,'y':np.array(y)}, index=groups)
    # remove single entry groups
    dat = dat[dat.groupby('groups').x.transform(len) > min_cells]
    group_widhts = list(dat.groupby('groups').apply(lambda x: 0.9*(x.max()-x.min())).x)
    group_position = dat.groupby('groups').mean().x
    group_position = list(group_position.round(2))

    datlist = list()
    uni_groups = np.unique(dat['groups'])
    uni_groups.sort()
    for grp in uni_groups:
        datlist.append(dat['y'].loc[grp])

    if ax is not None:
        plt.sca(ax)
    violins = plt.violinplot(datlist, group_position, widths=group_widhts, showextrema=False)

    for vio in violins['bodies']:
        if color is not None:
            vio.set_color(color)
            vio.set_alpha(alpha)

    # sns.violinplot(x='x', y='y',hue=None, data=dat, positions=group_position, names=None, widths=group_widhts, ax=ax,
    #                alpha=alpha, inner="quartile", **kwargs)
    xlim = np.arange(0, int(np.round(xmax))+1, 1.0)
    while len(xlim) > 5:
       xlim = xlim[np.arange(0,len(xlim),2)]
    if ax is not None:
        plt.sca(ax)
    plt.xticks(xlim, xlim)
    plt.ylim(0,max(y)+0.5)
    plt.xlim(0,xmax+1)


def plot_spline(x, y,**kwargs):
    x = pd.Series(x.values)
    y = pd.Series(y.values)
    spl = interpolate.UnivariateSpline(np.array(x),np.array(y), s=len(y) * np.var(y))
    x2 =  np.array(x.copy())
    x2.sort()
    x2 = x2[np.where(np.isfinite(spl(x2) ))]
    plt.plot(x2,spl(x2),**kwargs)

#
def pairplot_scatter(data):
    g = sns.PairGrid(data, diag_sharey=False)
    g.map_lower(lib.hexbin)
    g.map_upper(plt.scatter, marker='.', alpha=0.1)
    g.map_diag(sns.kdeplot, lw=1)
    return g

def pairpot_violin(data, bin_transform=None):
    g = sns.PairGrid(data, diag_sharey=False)
    g.map_lower(lib.violin_slice, bin_transform=lib.transform_arcsinh)
    g.map_upper(lib.violin_slice, bin_transform=lib.transform_arcsinh)
    g.map_upper(plt.scatter,marker='.', alpha=0.1)
    g.map_lower(plt.scatter,marker='.', alpha=0.1)
    g.map_diag(sns.kdeplot, lw=1)

#
def make_channels_unique(chan, sep='-'):
    "Makes the strings in a list unique by adding sep "
    uni_chan = set(chan)
    for cur_chan in uni_chan:
        n=0
        for i,t_chan in enumerate(chan):
            if cur_chan == t_chan:
                n += 1
                if n > 1:
                    chan[i] = chan[i]+sep+str(n)
    return chan

def get_n_bins(x, nbins):
    bins = np.linspace(0, x.max(), nbins+1)[1:-1]
    return bins

def get_n_minmax_bins(xmin, xmax, nbins):
    bins = np.linspace(xmin, xmax, nbins+1)[1:-1]
    return bins

def get_bin_stat(dat, x, fkt='median', nbins=20, bins=None):
    if bins is None:
        bins = get_n_bins(x, nbins)
    groups = pd.Series(np.digitize(x, bins))
    dat.index = groups
    if fkt == 'median':
        dat_summary = dat.groupby(dat.index).median()
    elif fkt == 'mean':
        dat_summary = dat.groupby(dat.index).mean()
    else:
        dat_summary = dat.groupby(dat.index).aggregate(fkt)
    return dat_summary


def get_piecewise_approx(tdat, x, nbins=20, bins=None):
    """
    This function will calculate a piecewise linear approximation of the data,


    :param tdat: a pandas dataframe
    :param x: a vector of the x values, used for binning
    :param nbins: number of bins that will be used
    :return: slopes: numpy array with the slopes
    """

    col_names = tdat.columns
    if bins is None:
        bins = get_n_bins(x, nbins)


    dbins = bins[1]-bins[0]
    ad_bins = bins-(dbins/2)
    ad_bins = list(ad_bins)
    ad_bins.append(ad_bins[-1]+dbins)
    ad_bins = np.array(ad_bins)

    groups = pd.Series(np.digitize(x, ad_bins))
    tdat.index = groups
    x.index = groups
    dat_bin = np.array(tdat.groupby(tdat.index).median())
    x_bin = x.groupby(tdat.index).median()
    y_diffs = np.diff(dat_bin, axis=0)
    x_diffs = np.diff(x_bin)
    slopes = y_diffs/x_diffs[:,None]
    intercept = dat_bin[:-1, :]-(slopes * x_bin[:-1, None])

    # predict the data for all the groups over which the data was predicted
    groups = pd.Series(np.digitize(x, bins))
    tdat.index = groups
    x.index = groups
    overall_mean= np.mean(dat_bin,0)
    overall_var = list()
    fit_var = list()
    fit_med_var = list()
    median_bin = list()
    counts = list()
    for g in range(nbins):

        dat_pred = intercept[g]+ slopes[g, :].T * x.loc[g][:, None]
        delta_pred = tdat.loc[g]-dat_pred
        delta_pred *= delta_pred
        n = delta_pred.shape[0]
        var = delta_pred.sum()/n

        delta_pred = tdat.loc[g]-tdat.loc[g].median()
        delta_pred *= delta_pred
        n = delta_pred.shape[0]
        var_med = delta_pred.sum()/n

        med = tdat.loc[g].median()

        ov_var = tdat.loc[g]-overall_mean
        ov_var *= ov_var
        ov_var = ov_var.sum()/n

        fit_var.append(var)
        fit_med_var.append(var_med)
        overall_var.append(ov_var)
        median_bin.append(med)
        counts.append(n)


    output = (fit_var, fit_med_var, overall_var, median_bin, counts)
    output = [np.vstack(l) for l in output]
    output = tuple([slopes] + [intercept] + output)
    return output


def get_binned_approx(y_val, x_val, nbins=20, bins=None, minCells=0, stat='median'):
    """
    Calculate the bin statistics

    :param y_val: a pandas series
    :param x_val:
    :param nbins:
    :param bins:
    :param minCells:
    :param stat:
    :return:
    """
    if bins is None:
        bins = get_n_bins(x_val, nbins)

    # group by bins
    groups = pd.Series(np.digitize(x_val, bins))
    y_val = y_val.set_index(groups)

    # filter out to small groups
    group_size = y_val.groupby(level=0).size()
    group_size = group_size[group_size > minCells]
    groups = group_size.index.get_level_values(0)
    y_val = y_val.loc[groups]

    # calculate bin statistics
    if stat == 'median':
        bin_med = y_val.groupby(level=0).median()
    elif stat == 'mean':
        bin_med = y_val.groupby(level=0).mean()

    # calculate overall mean/variance
    overall_mean= bin_med.mean()

    # square the differences
    tmp = (y_val.sub(overall_mean, axis=1)) ** 2
    overall_var = tmp.groupby(level=0).mean().stack()

    # square the differences
    tmp = (y_val.sub(bin_med.loc[groups])) ** 2
    fit_med_var = tmp.groupby(level=0).mean().stack()

    # save the medians
    bin_med = bin_med.stack()
    bin_med.index.names = ['bin', 'target']
    # save the number of cells per group
    counts = np.array(group_size.loc[bin_med.index.get_level_values('bin')])

    # put it together into a output data frame
    out_dat = pd.DataFrame(data={stat: bin_med, 'overall_var': overall_var, 'fit_var': fit_med_var, 'counts': counts},
                           index=bin_med.index)

    return out_dat

def column_to_hiercol(data, toplevelname, names):
    col = pd.MultiIndex.from_tuples([(toplevelname, col) for col in data.columns], names=names)
    data.columns = col




#
# def get_spline_approx(dat2, x, nbins=20):
#     """
#     This function will calculate a piecewise linear approximation of the data,
#
#
#     :param dat: a pandas dataframe
#     :param x: a vector of the x values, used for binning
#     :param nbins: number of bins that will be used
#     :return: slopes: numpy array with the slopes
#     """
#     x = dat.loc[idx][chan].copy()
#     dat2 = dat.loc[idx].copy()
#     col_names = dat2.columns
#
#     bins = np.linspace(0, x.max(), nbins+1)[1:]
#     dbins = bins[1]-bins[0]
#     bins -= dbins/2
#
#     groups = pd.Series(np.digitize(x, bins))
#     dat.index = groups
#     x.index = groups
#
#
#     col = dat2.columns[32]
#     #spl = interpolate.LSQUnivariateSpline(np.array(x),np.array(dat[col]),t=bins)
#
#     for col in dat2.columns:
#         plt.figure()
#         y = dat2[col]
#         spl = interpolate.UnivariateSpline(np.array(x),np.array(y), s=len(y) * np.var(y))
#         x2 =  np.array(x.copy())
#         x2.sort()
#         x2 = x2[np.where(np.isfinite(spl(x2) ))]
#         plt.plot(x2,spl(x2))
#         plt.scatter(np.array(x), np.array(dat2[col]), alpha=0.3)
#
#         spl_der1 = spl.derivative()
#         spl_der1.integral(x.min(),x.max())
#         spl_der2 = spl.derivative(n=2)
#         spl_der2.integral(x.min(),x.max())

def prefix_channels(dat, prefix, sel_channels=None, sep='_', inplace=True):
    if sel_channels is None:
        sel_channels = dat.columns
    col_dict = {chan: prefix+sep+chan for chan in sel_channels}
    if inplace:
        dat.rename(columns=col_dict, inplace=inplace)
    else:
        return dat.rename(columns=col_dict, inplace=inplace)

def get_var_to_ref(x, ref):
    res = (((x-ref)*(x-ref)).sum(0)).sqrt(0)
    return res

def find_fdr_cutoff(values, neg_values, fdr=0.05):
    """
    calculates the empirical FDR cuttof using a series of values values and a true negative label is_tn
    :param value:
    :param is_tn:
    :return:
    """
    univals = list(values.tolist())
    univals.extend(neg_values.tolist())
    univals.sort(reverse=True)
    n_neg = len(neg_values)
    n_vals = len(values)
    for i, co in enumerate(univals):
        if n_neg > 0:
            fpr = np.sum(neg_values >= co)/n_neg
        else:
            n_neg
        if fpr > 0:
            e_fdr = (fpr*n_vals)/np.sum(values >= co)
            if e_fdr > fdr:
                return univals[i-1]
    return co

def format_filename(s):
    """
    From a gist:
    Take a string and return a valid filename constructed from the string.
Uses a whitelist approach: any characters not present in valid_chars are
removed. Also spaces are replaced with underscores.

Note: this method may produce invalid filenames such as ``, `.` or `..`
When I use this method I prepend a date string like '2009_01_15_19_46_32_'
and append a file extension like '.txt', so I avoid the potential of using
an invalid filename.

"""
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    filename = ''.join(c for c in s if c in valid_chars)
    filename = filename.replace(' ','-') # I don't like spaces in filenames.
    return filename


# def sort_linkage(Z, order):
#     """
#     Sorts the linkage matrix such that it should be correctly ordered when plotted
#     without reordering.
#     :param Z: linkage matrix
#     :param order: order of the leaves
#     :return: ordered linkage matrix Z
#     """
#
#     Z_ord = Z.copy()
#     for idx in order:
#


def part_cor_coef(corr_mat, ref_name):
    # converts the correlation matrix to a partial correlation matrix in respect to the column ref_name
    out_mat = corr_mat.copy()
    for col in corr_mat.columns:
        for row in corr_mat.index:
            xy = corr_mat.loc[row, col]
            xu = corr_mat.loc[row, ref_name]
            yu =corr_mat.loc[col, ref_name]
            out_mat.loc[row, col] = (xy-(xu*yu))/np.sqrt((1-xu**2)*(1-yu**2))


    out_mat.drop(ref_name, axis=0, inplace=True)
    out_mat.drop(ref_name, axis=1, inplace=True)
    return(out_mat)


def hoeffings_d(x,y):
    x = np.array(x)
    y = np.array(y)
    N = len(x);
    r_x = stats.rankdata(x, method='average')
    r_y = stats.rankdata(y)

    assert(len(x) == len(y))
    Q = np.zeros(x.size)
    for i in range(N):
        Q[i] = 1 + np.sum( (r_x < r_x[i]) & (r_y < r_y[i]) )
        # and deal with cases where one or both values are ties, which contribute less
        Q[i] = Q[i] + 1/4 * (np.sum( (r_x == r_x[i]) & (r_y == r_y[i]) ) - 1) # both indices tie.  -1 because we know point i matches
        Q[i] = Q[i] + 1/2 * np.sum( (r_x == r_x[i]) & (r_y < r_y[i]) ) # one index ties.
        Q[i] = Q[i] + 1/2 * np.sum( (r_x < r_x[i]) & (r_y == r_y[i])) # one index ties.


    D1 = np.sum( (Q-1)*(Q-2) )
    D2 = np.sum( (r_x-1)*(r_x-2)*(r_y-1)*(r_y-2) )
    D3 = np.sum( (r_x-2)*(r_y-2)*(Q-1) )

    D = 30*((N-2)*(N-3)*D1 + D2 - 2*(N-2)*D3) / (N*(N-1)*(N-2)*(N-3)*(N-4))


def dremi(x,y):
    pass

def simulate_perfect_distribution(n, dist_type=None):
    """
    A function to produce different distributions
    
    :param n: number of points of the distribution
    :param dist_type: The name for the distribution to be simulated
    :return: a dict with the distribution and a numpy array with n values
    """

    
    if 'independent' == dist_type:
        # Independent
        x = np.random.normal(size=n)
        y= np.array(np.arange(n))

    elif 'linear' == dist_type:
        x =  np.array(np.arange(-2, 2+4/(n-1), 4/(n-1))*0.5)
        y = np.array(np.arange(n))

    elif 'exponential' == dist_type:
        x = np.exp(np.arange(-30, 20+50/(n-1), 50/(n-1))) * 0.01
        y = np.array(np.arange(n))

    elif 'quadratic' == dist_type:
        x = np.array(np.arange(-3, 3+6/(n-1), 6/(n-1)))
        x = np.multiply(x,x)
        y = np.array(np.arange(n))

    elif 'sine' == dist_type:
        x = 2*np.sin(np.arange(0, 10+10/(n-1), 10/(n-1)))
        y = np.array(np.arange(n))
        
    elif 'circumference' == dist_type:
        m = np.mod(n,2) + n
        x = np.zeros(m)
        y = np.zeros(m)
        x[:(m/2)] = np.arange(-5, 5+10/(m/2-1), 10/(m/2-1))
        x[(m/2 ):m] = np.arange(5, -5+-10/(m/2 - 1), -10/(m/2 - 1))
        y[:(m/2)] = np.sqrt(25 - x[:(m/2)]*x[:(m/2)])
        y[(m/2 ):m] = -np.sqrt(25 - x[(m/2 ):m]*x[(m/2 ):m])
        idx = np.random.choice(range(m), n, replace=False)
        x = x[idx]
        y = y[idx]


    elif 'cross' == dist_type:
        m =  np.mod(n,2) + n
        x = np.zeros(m)
        y = np.zeros(m)
        x[:(m/2)-1] = np.arange(-5, 5, 10/(m/2-1))
        x[(m/2 ):m-1] = np.arange(5, -5, -10/(m/2 - 1))
        y[:(m/2)] = x[0:(m/2)]
        y[(m/2 ):m] = -x[(m/2 ):m]

        idx = np.random.choice(range(m), n, replace=False)
        x = x[idx]
        y = y[idx]

    elif 'square' == dist_type:
        m =  4-np.mod(n,4) + n

        y = np.arange(m)
        x = np.arange(m)
        x[:(m/4)] = np.arange(6, 9+3/(m/4-1), 3/(m/4-1))
        x[(m/4 ):(m/2)] = np.repeat(9, m/4)
        x[(m/2):(3*m/4)] = np.arange(9, 6+3/(m/4), -3/(m/4))
        x[(3*m/4):] = np.repeat(6, m/4)

        y[:(m/4)] = np.repeat(6, m/4)
        y[(m/4 ):(m/2)] = np.arange(6, 9+3/(m/4-1), 3/(m/4-1))
        y[(m/2):(3*m/4)] = np.repeat(9, m/4)
        y[(3*m/4):] = np.arange(9, 6+3/(m/4), -3/(m/4))

        idx = np.random.choice(range(m), n, replace=False)
        x = x[idx]
        y = y[idx]

    else:
        raise NameError('No valid distribution type!')

    return x, y


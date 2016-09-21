
# coding: utf-8

# # BpR2 vs number of bins
# 
# This should show how bpR2 changes when changing the number of bins

# In[61]:

from  __future__ import division
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import numpy as np
from matplotlib_venn import venn3, venn3_circles
import bokeh.plotting as bplt
import bokeh.models as bmod
import bokeh as bk
from scipy import stats
bk.io.output_notebook()
get_ipython().magic('matplotlib inline')


# In[10]:

folder_dict ={
    4:'nbin_experiment_4_2.5perc_bpr2_median_25_final',
    5:'nbin5_2.5perc_bpr2_median_25_test',
    8:'nbin_experiment_8_2.5perc_bpr2_median_25_final',
    10: 'nbin10_2.5perc_bpr2_median_25_final',
    12: 'nbin_experiment_12_2.5perc_bpr2_median_25_final',
    15:'nbin_experiment_15_2.5perc_bpr2_median_25_final',
    20: 'nbin20_2.5perc_bpr2_median_25_test',
    25: 'nbin_experiment_25_2.5perc_bpr2_median_25_final',
    30: 'nbin_experiment_30_2.5perc_bpr2_median_25_final',
    40: 'nbin40_2.5perc_bpr2_median_25_test',
    60: 'nbin_experiment_60_2.5perc_bpr2_median_25_final',
    80: 'nbin80_2.5perc_bpr2_median_25_test',
    100:  'nbin_experiment_100_2.5perc_bpr2_median_25_final'}
    
base_folder = '/mnt/imls-bod/Xiao-Kang/EGF transfection/plots/'

file_name = 't_bindat'  

        


# In[11]:

tdat = pd.read_pickle(os.path.join(base_folder,folder_dict[5],file_name))

bp_dict = dict((key,pd.read_pickle(os.path.join(base_folder,fol,file_name))[('stats', 'mean_mean_var_ratio')]) 
               for key, fol in folder_dict.items())

bp_dat = pd.DataFrame.from_dict(bp_dict)
bp_dat = bp_dat.stack().reset_index(-1,drop=False)
bp_dat.columns = ['nbins', 'bpR2']
#bp_dat = bp_dat.reset_index(drop=False)


# In[12]:

grouped = bp_dat.groupby(level=['experiment','marker','perturbation','row_col','timepoint','origin', 'target'])

for name, group in grouped:
    plt.plot(group['nbins'], group['bpR2'], alpha=0.1)

plt.show()


# Try to divide it by tp=10 to show the deviations

# In[39]:

bp_dat['bpR2_vs_10'] = np.array(bp_dat.groupby(level=['experiment','marker','perturbation','row_col','timepoint','origin', 'target']
                                   ).apply(lambda g: g['bpR2']/np.float(g['bpR2'][g['nbins'] == 10])))

bp_dat['bpR2_at_10'] = np.array(bp_dat.groupby(level=['experiment','marker','perturbation','row_col','timepoint','origin', 'target']
                                   ).apply(lambda g: g['bpR2']/g['bpR2']*np.float(g['bpR2'][g['nbins'] == 10])))
                                
bp_dat['bpR2_at_20'] = np.array(bp_dat.groupby(level=['experiment','marker','perturbation','row_col','timepoint','origin', 'target']
                                   ).apply(lambda g: g['bpR2']/g['bpR2']* np.float(g['bpR2'][g['nbins'] == 20])))


# In[43]:

#bp_dat.hist(column='bpR2_vs_10',by='nbins')

bp_dat.groupby('nbins').describe()


# Only look at bpR2 that are over 0.05 in at least one binning

# In[67]:

thresh = 0.05
fil = np.array(bp_dat['bpR2'].groupby(level=['experiment','marker','perturbation','row_col','timepoint','origin', 'target']
                                   ).transform(lambda g: np.any(g > thresh)))

bp_dat_large = bp_dat.loc[fil == 1]


# In[49]:

bp_dat_large.query('bpR2_vs_10 > 1.5')
print('upper 99')
print(bp_dat_large.groupby('nbins')['bpR2_vs_10'].apply(lambda x: np.percentile(x, 99)))
print('upper 95')
print(bp_dat_large.groupby('nbins')['bpR2_vs_10'].apply(lambda x: np.percentile(x, 95)))

print('\nlower 05')
print(bp_dat_large.groupby('nbins')['bpR2_vs_10'].apply(lambda x: np.percentile(x, 5)))
print('\nlower 01')
print(bp_dat_large.groupby('nbins')['bpR2_vs_10'].apply(lambda x: np.percentile(x, 1)))

g = sns.FacetGrid(bp_dat_large, col='nbins')
g.map(plt.scatter, 'bpR2_at_20' , 'bpR2', size=1)
#sns.plt.ylim(0,1)


# In[68]:

bp_mat = bp_dat[['nbins','bpR2']].pivot(columns='nbins')


def corrfunc(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f}".format(r),
                xy=(.1, .8), xycoords=ax.transAxes)

g = sns.PairGrid(bp_mat)
g.map(plt.scatter)
g.map(corrfunc)


# In[ ]:




# In[131]:


folder_dict ={
            1:'mincell_experiment_10_2.5perc_bpr2_median_1_final',
        5:'mincell_experiment_10_2.5perc_bpr2_median_5_final',
    10:'mincell_experiment_10_2.5perc_bpr2_median_10_final',
        15:'mincell_experiment_10_2.5perc_bpr2_median_15_final',
        20:'mincell_experiment_10_2.5perc_bpr2_median_20_final',
    25:'mincell_experiment_10_2.5perc_bpr2_median_25_final',
        30:'mincell_experiment_10_2.5perc_bpr2_median_30_final',

        40:'mincell_experiment_10_2.5perc_bpr2_median_40_final',

        60:'mincell_experiment_10_2.5perc_bpr2_median_60_final',
    80:'mincell_experiment_10_2.5perc_bpr2_median_80_final',
    100:'mincell_experiment_10_2.5perc_bpr2_median_100_final',
        150:'mincell_experiment_10_2.5perc_bpr2_median_150_final'
}
    
base_folder = '/mnt/imls-bod/Xiao-Kang/EGF transfection/plots/'

file_name = 'bindat'  


#tdat = pd.read_pickle(os.path.join(base_folder,folder_dict[5],file_name))

bp_dict = dict((key,pd.read_pickle(os.path.join(base_folder,fol,file_name))[('stats', 'mean_mean_var_ratio')]) 
               for key, fol in folder_dict.items())

bp_dat = pd.DataFrame.from_dict(bp_dict)
bp_dat = bp_dat.stack().reset_index(-1,drop=False)
bp_dat.columns = ['mincells', 'bpR2']
#bp_dat = bp_dat.reset_index(drop=False)


# In[133]:

bp_mat = bp_dat[['mincells','bpR2']].pivot(columns='mincells')

def corrfunc(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f}".format(r),
                xy=(.1, .8), xycoords=ax.transAxes)

g = sns.PairGrid(bp_mat, diag_sharey=False)
g.map(plt.scatter)
g.map(corrfunc)
#g.set(ylim=(0, 0.5))
#g.set(xlim=(0, 0.5))


# In[134]:

thresh = 0.3
fil = np.array(bp_dat['bpR2'].groupby(level=['experiment','marker','perturbation','row_col','timepoint','origin', 'target']
                                   ).transform(lambda g: np.any(g < thresh)))

bp_dat_large = bp_dat.loc[fil == 1]

bp_mat = bp_dat_large[['mincells','bpR2']].pivot(columns='mincells')

def corrfunc(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.3f}".format(r),
                xy=(.1, .8), xycoords=ax.transAxes)

g = sns.PairGrid(bp_mat, diag_sharey=False)
g.map(plt.scatter)
g.map(corrfunc)
g.set(ylim=(0, 0.5))
g.set(xlim=(0, 0.5))


# In[97]:

bp_dat.loc[np.isnan(bp_dat['bpR2'])]


# In[98]:

np.all(np.isfinite(bp_dat['bpR2']))


# In[99]:

bp_dat.xs('Beta-catenin', level='target').query('mincells >5')


# In[88]:

bp_mat


# In[ ]:




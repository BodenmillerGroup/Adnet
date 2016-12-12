
# coding: utf-8

# # Compare Readouts vs bpR2
# 
# by Vito Zanotelli
# vito.zanotelli@gmail.com
# 
# Look at the overlap between using DREMI/Spearman correlatio/Pearson correlation as a cutoff.
# 
# The bpR2 was calculated using the adnet/Spearman/Pearson correlatio were calulated using the Adnet script.
# The DREMI value was calculated using the reference matlab GUI implementation from the DREMI paper (Conditional density-based analysis of T cell signaling in single-cell data, Science 2014) modified to write out the results as a .csv file.
# 
# 

# In[1]:

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
bk.io.output_notebook()
get_ipython().magic('matplotlib inline')


# Set the file paths

# In[2]:

# the bin_dat pickle file, outputed by the adnet analysis script
bin_dat_fn = '/mnt/imls-bod/XiaoKangL/EGF transfection/plots/nbin10_2.5perc_bpr2_median_25_final/t_bindat'
# exported from the matlab dremi gui
dremi_fn = '/home/vitoz/imls-bod/XiaoKangL/EGF transfection/benchmark/20160314_dremi_values_all_overexpressions.csv'
# a dictionary with the old and new names
name_dict = '/home/vitoz/imls-bod/XiaoKangL/EGF transfection/name_dict.csv'
out_folder = '/home/vitoz/Data/Analysis/XKL'
neg_ctrl_names = ['empty-1','empty-2','GFP-FLAG-1','GFP-FLAG-2']

#neg_ctrl_names = ['empty-1','empty-2']


crap_names = ['cleaved PARP-cleaved caspase3', 'cyclin B1', 'p-4EBP1', 'p-HH3', 'p-RB', 'beads']


# Load the binned data file

# In[3]:

bin_dat = pd.read_pickle(bin_dat_fn)
bin_dat.index.get_level_values('target').unique()


# Prepare the name dict
# 

# In[4]:

name_dict = pd.read_csv(name_dict)
name_dict = {row['old']: row['new'] for idx, row in name_dict.iterrows()}


# Read the DREMI data

# In[5]:

dat_dremi = pd.read_csv(dremi_fn, sep=',',index_col=False)
dat_dremi.head()


# make the name dict compatible with the DREMI names (only alphanumeric, lower case)

# In[6]:

name_dict_dremi = dict((filter(str.isalnum, oldname.lower()), nicename) for oldname, nicename in name_dict.iteritems())


# In[7]:

dat_dremi['target'] = dat_dremi['target'].map(lambda x: name_dict_dremi[x])


# Merge the dremi data with the bpr2 data

# In[8]:

dat_dremi_stacked = dat_dremi.copy()
dat_dremi_stacked['origin'] = dat_dremi_stacked['origin'].map(lambda x: x.upper())
dat_dremi_stacked['origin'].map(lambda x: x.upper())
dat_dremi_stacked = dat_dremi_stacked.set_index(['origin','target'])

dat_dremi_stacked =  pd.DataFrame(dat_dremi_stacked.stack(), columns=['dremi'])
dat_dremi_stacked.index.names = ['origin', 'target', 'filename']

# extract experiment and row-col from the filename

get_experiment = lambda x: x.split('/')[5]
get_rowcol = lambda x: x.split('/')[5]
dat_dremi_stacked['experiment'] = [x.split('/')[5] for x in dat_dremi_stacked.index.get_level_values('filename')]
dat_dremi_stacked['row_col'] = [x.split('_')[1] for x in dat_dremi_stacked.index.get_level_values('filename')]

dat_dremi_stacked = dat_dremi_stacked.reset_index('filename', drop=True)
dat_dremi_stacked = dat_dremi_stacked.reset_index()

#pd.merge(bin_dat, dat_dremi_stacked, left_index=['origin', 'target', 'row_col', 'experiment'], right_index=['origin', 'target', 'row_col', 'experiment'])

t_dat = bin_dat.reset_index()
t_dat = t_dat[bin_dat.index.names]
t_dat.columns = t_dat.columns.get_level_values(0)
dat_dremi_stacked = pd.merge(t_dat, dat_dremi_stacked, how='outer')

# due to the way the merging was performed, both bindat and dat_dremi are now aligned

dat_dremi_stacked = dat_dremi_stacked.dropna(subset=['marker'])
bin_dat[('stats', 'dremi')] = dat_dremi_stacked['dremi'].tolist()



# Calculate the dremi median over the replicates

# In[9]:

bin_dat[('stats', 'dremi_median')] = bin_dat[('stats', 'dremi')].groupby(level=['marker', 'origin', 'target', 'timepoint','perturbation']).transform(np.median)


# Prepare the negative control filter

# In[10]:

bin_dat = bin_dat.loc[bin_dat.index.get_level_values('target').isin(crap_names) == False, :]
neg_mark_fil = bin_dat.index.get_level_values('marker').isin(neg_ctrl_names)


# Filter out the non wanted markers

# ### Start the comparison
# 

# Look only at timepoint 0

# In[11]:

#fil = [ x not in filtarget for x in bin_dat.index.get_level_values('target')]
fil = [ x ==0  for x in bin_dat.index.get_level_values('timepoint')]
bin_dat = bin_dat.loc[fil]
neg_mark_fil = np.array([negmark for negmark, f in zip(neg_mark_fil, fil) if f])
bin_dat.loc[neg_mark_fil == False].index.get_level_values('marker').unique()


# Calculate the cutoffs using the negative control for dremi, bpr2 and spearman

# In[12]:

bins = np.arange(0,1,0.025)

bin_dat.loc[neg_mark_fil == False, ('stats', 'dremi_median')].hist(normed=1,bins=bins)
bin_dat.loc[neg_mark_fil, ('stats', 'dremi_median')].hist(normed=1, alpha=0.4,bins=bins)

maxneg_dremi = bin_dat.loc[neg_mark_fil, ('stats', 'dremi_median')].max()
maxneg_dremi_99 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'dremi_median')].dropna(),99,)
maxneg_dremi_90 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'dremi_median')].dropna(),90,)
print(bin_dat.loc[neg_mark_fil == False, ('stats', 'dremi_median')].max())
print(bin_dat.loc[neg_mark_fil, ('stats', 'dremi_median')].max())

print(maxneg_dremi)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'dremi_median')] > maxneg_dremi)/3)
print(maxneg_dremi_99)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'dremi_median')] > maxneg_dremi_99)/3)


# In[13]:

bin_dat['stats']
bins = np.arange(0,1,0.025)

bin_dat.loc[neg_mark_fil == False, ('stats', 'median_mean_var_ratio')].hist(normed=1,bins=bins)
bin_dat.loc[neg_mark_fil, ('stats', 'median_mean_var_ratio')].hist(normed=1, alpha=0.4,bins=bins)

maxneg_bp = bin_dat.loc[neg_mark_fil, ('stats', 'median_mean_var_ratio')].max()
maxneg_bp_99 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'median_mean_var_ratio')].dropna(),99)
maxneg_bp_90 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'median_mean_var_ratio')].dropna(),90)

print(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_mean_var_ratio')].max())
print(bin_dat.loc[neg_mark_fil, ('stats', 'median_mean_var_ratio')].max())

print(maxneg_bp)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_mean_var_ratio')] > maxneg_bp)/3)
print(maxneg_bp_99)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_mean_var_ratio')] > maxneg_bp_99)/3)


# In[14]:

maxneg_bp_99


# In[15]:

bin_dat['stats']
bins = np.arange(0,1,0.025)

bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_spearman_overall')].hist(normed=1,bins=bins)
bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_spearman_overall')].hist(normed=1, alpha=0.4,bins=bins)

maxneg_sp = bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_spearman_overall')].max()
maxneg_sp_99 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_spearman_overall')].dropna(),99,)
maxneg_sp_90 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_spearman_overall')].dropna(),90,)
print(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_spearman_overall')].max())
print(bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_spearman_overall')].max())

print(maxneg_sp)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_spearman_overall')] > maxneg_sp)/3)
print(maxneg_bp_99)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_spearman_overall')] > maxneg_sp_99)/3)


# In[16]:

bin_dat['stats']
bins = np.arange(0,1,0.025)

bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_pearson_overall')].hist(normed=1,bins=bins)
bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_pearson_overall')].hist(normed=1, alpha=0.4,bins=bins)

maxneg_pc = bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_pearson_overall')].max()
maxneg_pc_99 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_pearson_overall')].dropna(),99,)
maxneg_pc_90 = np.percentile(bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_pearson_overall')].dropna(),90,)
print(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_pearson_overall')].max())
print(bin_dat.loc[neg_mark_fil, ('stats', 'median_abs_corr_pearson_overall')].max())

print(maxneg_pc)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_pearson_overall')] > maxneg_sp)/3)
print(maxneg_pc_99)
print(np.sum(bin_dat.loc[neg_mark_fil == False, ('stats', 'median_abs_corr_pearson_overall')] > maxneg_sp_99)/3)


# In[17]:

#fil = bin_dat[('stats', 'dremi_median')] > 0.24
#fildat = bin_dat.loc[(fil) & (neg_mark_fil)]
#hits_dremi_GFP = set('_'.join([m, t, str(tp)]) for m, t, tp in
#             zip(fildat.reset_index()['marker'],fildat.reset_index()['target'],fildat.reset_index()['timepoint']))
#hits_dremi_GFP


# Compare the overlap of hits between 

# In[18]:

fil = bin_dat[('stats', 'dremi_median')] > maxneg_dremi_99
fildat = bin_dat.loc[(fil) & (neg_mark_fil == False)]

hits_dremi = set('_'.join([m, t, str(tp)]) for m, t, tp in
             zip(fildat.reset_index()['marker'],fildat.reset_index()['target'],fildat.reset_index()['timepoint']))

fil = bin_dat[('stats', 'median_mean_var_ratio')] > maxneg_bp_99
fildat = bin_dat.loc[(fil) & (neg_mark_fil == False)]

hits_bp = set('_'.join([m, t, str(tp)]) for m, t, tp in
             zip(fildat.reset_index()['marker'],fildat.reset_index()['target'],fildat.reset_index()['timepoint']))

fil = bin_dat[('stats', 'median_abs_corr_spearman_overall')] > maxneg_sp_99

fildat = bin_dat.loc[(fil) & (neg_mark_fil == False)]

hits_sp = set('_'.join([m, t, str(tp)]) for m, t, tp in
             zip(fildat.reset_index()['marker'],fildat.reset_index()['target'],fildat.reset_index()['timepoint']))

fil = bin_dat[('stats', 'median_abs_corr_pearson_overall')] > maxneg_pc_99

fildat = bin_dat.loc[(fil) & (neg_mark_fil == False)]

hits_pc = set('_'.join([m, t, str(tp)]) for m, t, tp in
             zip(fildat.reset_index()['marker'],fildat.reset_index()['target'],fildat.reset_index()['timepoint']))


# Look at the different hits as a venn diagramm

# In[19]:


venn3([hits_bp,hits_sp,hits_dremi],set_labels=['bp-R2', 'Spearman', 'DREMI'])

#plt.savefig(os.path.join(out_folder,'20160429_readout_comparison_venn_wocc.pdf'))


# In[20]:

venn3([hits_bp,hits_pc,hits_dremi],set_labels=['bp-R2', 'Pearson', 'DREMI'])


# ## Look how exactly the hits differ

# In[21]:

hits_bp


# In[22]:

hits_bp.difference(hits_dremi)


# In[23]:



hits_bp.difference(hits_dremi.union(hits_sp))


# In[24]:

hits_bp.difference(hits_sp)


# In[25]:

hits_dremi.difference(hits_bp)


# In[26]:

hits_dremi.difference(hits_bp.union(hits_sp))


# In[27]:

hits_dremi.difference(hits_sp)


# In[28]:

hits_sp.difference(hits_bp)


# In[29]:

hits_sp.difference(hits_bp.union(hits_dremi))


# In[30]:

hits_bp.union(hits_sp).difference(hits_dremi)


# In[31]:

hits_bp.intersection(hits_dremi).difference(hits_sp)


# ## Make a table of the hits and save it 

# In[32]:

hit_dict = dict()

hit_dict['all_hits'] = list(set(list(hits_sp)+ list(hits_dremi)+ list(hits_bp)))
hit_dict['is_bp'] = [h in hits_bp for h in hit_dict['all_hits'] ]
hit_dict['is_sp'] = [h in hits_sp for h in hit_dict['all_hits'] ]
hit_dict['is_dremi'] = [h in hits_dremi for h in hit_dict['all_hits'] ]


hit_tab = pd.DataFrame.from_dict(hit_dict)
#hit_tab = hit_tab.set_index('all_hits')

hit_tab = hit_tab.sort_values(by=['is_bp','is_dremi', 'is_sp'],ascending=False)

hit_tab.to_csv(os.path.join(out_folder,'20160429_readout_comparison_vsemptygfp_woCC.csv'),index=False)


# In[33]:

hit_tab


# In[34]:

#bin_dat.plot(kind='scatter', x=('stats', 'mean_var_ratio'), y=('stats', 'dremi'))



def draw_tooltip_scatter(df, x, y, plot_width=400, plot_height=400, pointsize=8, title='Title'):
    source = bplt.ColumnDataSource(data=dict([
        ('idx', bin_dat.index.tolist()),
        ('x',df[x]),
        ('y',df[y])]
        ))

    hover = bmod.HoverTool(
            tooltips=[
                ("idx", "@idx"),
            ("(x,y)", "($x, $y)")
            ]
        )

    p = bplt.figure(plot_width=plot_width, plot_height=plot_height, 
               title=title,
                   x_axis_label=str(x), y_axis_label=str(y))
    p.add_tools(hover)

    p.circle('x','y', size=pointsize, source=source)
    bplt.show(p)
    return p

xcol =('stats', 'mean_var_ratio')
ycol = ('stats', 'dremi')

draw_tooltip_scatter(bin_dat,xcol, ycol )


# In[35]:

#bin_dat.plot(kind='scatter', x=('stats', 'median_mean_var_ratio'), y=('stats', 'median_abs_corr_spearman_overall'))

xcol =('stats', 'mean_var_ratio')
ycol = ('stats', 'corr_spearman_overall')
draw_tooltip_scatter(bin_dat,xcol, ycol )


# In[36]:

#bin_dat.plot(kind='scatter', x=('stats', 'median_mean_var_ratio'), y=('stats', 'median_abs_corr_spearman_overall'))

xcol =('stats', 'mean_var_ratio')
ycol = ('stats', 'corr_pearson_overall')
draw_tooltip_scatter(bin_dat,xcol, ycol )


# Compare different bins

# In[37]:

bin_dat.plot(kind='scatter', x=('stats', 'mean_var_ratio'), y=('stats', 'corr_spearman_overall'))


# In[38]:

fil = (bin_dat[('stats', 'median_mean_var_ratio')] < 0.5) &     (bin_dat[('stats', 'median_mean_var_ratio')] > 0.4) & (bin_dat[('stats', 'median_abs_corr_spearman_overall')] <0.4)
    
bin_dat.loc[fil]


# In[39]:

bin_dat.xs('p-GSK3-Beta', level='target').xs('MAP2K6', level='marker')['stats'].columns


# In[40]:

np.median(1-(bin_dat.xs('p-GSK3-Beta', level='target').xs('MAP2K6', level='marker')['fit_var']/
   bin_dat.xs('p-GSK3-Beta', level='target').xs('MAP2K6', level='marker')[u'overall_var']), axis=1)


# In[41]:

print(bin_dat.index.get_level_values('target').unique())
print(bin_dat.index.get_level_values('marker').unique())


# In[ ]:




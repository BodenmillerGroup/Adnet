
# coding: utf-8

# ## Aim
# 
# Compare the coverage of the Signor network with the markers measured and the overexpressions performed
# 
# This version adapted to the reviewer comment:
# 
# "Reviewer #3:
# 
# Remarks to the Author:
# 
#    The revised manuscript submitted by Bodenmiller and colleagues addresses in a satisfactory manner most of the referees’ comments and criticisms with new controls, experiments and analyses.
# I only have the following additional remark I could not find a description of the algorithm used to obtain the shortest “directed” path in SIGNOR to be correlated with the experimental correlation.
# I am wondering whether in their analysis the authors took in consideration the sign associated to the path that should correlate with the sign of the experimental correlation. To make a simple example, let’s consider two pairs of proteins that correlate experimentally with a positive sign.  According to the signed directed network protein A directly activates the phosphorylation of B while protein C directly dephosphorylates protein D, In both cases the network distance is 1 but only the correlation of A with B is in accord with the experimental results
# "

# In[132]:

import networkx as nx
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from IPython.display import display, HTML
import os
get_ipython().magic('matplotlib inline')


# In[3]:

# The complete signor network
signor_fn = '/mnt/imls-bod/Xiao-Kang/EGF transfection/20160110_signor_all_data.tsv'


# In[4]:

# The antibody names in uniprot
antibody_fn = '/mnt/imls-bod/Xiao-Kang/EGF transfection/201610_antibody_uniprot.csv'


# In[5]:

# The overexpressions in uniprot
overexp_fn = '/mnt/imls-bod/Xiao-Kang/EGF transfection/201602_gene_uniprot.csv'


# In[159]:

outfolder ='/home/vitoz/imls-bod/Xiao-Kang/EGF transfection/benchmark'


# In[163]:

# ignored effects due to spillove
ignored_fn ='/home/vitoz/imls-bod/Xiao-Kang/EGF transfection/benchmark/20161020_ignored_relations.csv'


# In[218]:

# Effect mechanism filters for signor
not_used_effects = pd.read_csv(ignored_fn)


# In[7]:

signor_tab = pd.read_csv(signor_fn, sep='\t')


# In[20]:

signor_tab


# Read the network

# Make a dict for the names

# In[9]:

signor_dict = {k: v for k,v in zip(signor_tab['IdA'],signor_tab['EntityA'])}
signor_dict.update({k: v for k,v in zip(signor_tab['IdB'],signor_tab['EntityB'])})


#     Read the antibody and gene information

# In[10]:

overexp_tab = pd.read_csv(overexp_fn)
overexp_tab.head(2)

overexp_dict = {k: v for k,v in zip(overexp_tab['gene_xkl'],overexp_tab['UniProt'])}


# In[11]:

antibody_tab = pd.read_csv(antibody_fn)
antibody_tab.head(2)

antibody_dict= dict()
for ab, tab in antibody_tab.groupby('Antibody'):
    if not(tab['UniProt'].isnull().values.any()):
        antibody_dict[ab] = set(tab['UniProt'])
        
    


# In[40]:

antibody_tab


# In[60]:

antibody_dict_extended= dict()
for ab, tab in antibody_tab.groupby('Antibody'):
    if not(tab['UniProt'].isnull().values.any()):
       cur_entry = [(row['UniProt'],row['Residue'].split('/'), row['Phospho'])  for idx, row in tab.iterrows()]
       antibody_dict_extended.update({ab: cur_entry})


# In[61]:

antibody_dict_extended


# # Generate the graph

# It needs to be a multidigraph to allow for parallel edges

# In[12]:

g = nx.MultiDiGraph()

for idx, row in signor_tab.iterrows():
    g.add_edge(row['IdA'], row['IdB'], attr_dict=dict((col, row[col]) for col in ['Effect', 'EffectMechanism',
       'MechanismResidues', 'MechanismSequences','Direct', 'Sentence']))


# In[ ]:




# 
# 
# 
# g= nx.from_pandas_dataframe(signor_tab,'IdA','IdB',edge_attr=[u'Effect', u'EffectMechanism',
#        u'MechanismResidues', u'MechanismSequences',u'Direct'], create_using=nx.MultiDiGraph())

# Check if the graph contains all measured antibodies/overexpressions

# In[13]:

xkl_nodes = np.unique(overexp_tab['UniProt'].dropna().tolist() + antibody_tab['UniProt'].dropna().tolist())

all_nodes = g.nodes()
signor_xkl_nodes = [n for n in all_nodes if n in xkl_nodes]
len(signor_xkl_nodes) == len(xkl_nodes)


# -> All relevant species are present

# ## Check the average shortest path length between overexpressions and phosphoproteins

# Calculate pairwise shortest path distance, taking the minimum distanc if the antibody is not specific

# In[14]:

col = antibody_dict.keys()
row = overexp_tab['gene_xkl']

pair_sp = np.empty((len(row), len(col)))
pair_sp[:] = np.NAN
for i, over in enumerate(row):
    for j, ab in enumerate(col):
        sp = list()
        for ab_prot in antibody_dict[ab]:
            ov_prot = overexp_dict[over]
            try:
                sp.append(nx.shortest_path_length(g, ov_prot, ab_prot))
            except nx.NetworkXNoPath:
                pass
        if len(sp) > 0:
            pair_sp[i,j] = np.min(sp)
        else:
            pair_sp[i,j] = 10
 



# In[ ]:




# In[15]:

pair_dat = pd.DataFrame(pair_sp,columns=col, index=row)
pair_dat = pair_dat.stack()
pair_dat.index = pair_dat.index.rename(['marker','target'])


# Look at the distribution

# In[16]:

bins = np.arange(int(np.ceil(pair_dat.max()))+1)-0.5
pair_dat.hist(bins=bins)
pair_dat.describe()


# In[ ]:




# In[ ]:




# # As proposed by the reviewer: make a signed version
# 
# Check for each shortest path which sign is predicted.
# 
# Save two versions: positive and negative
# 
# Algorithm: 
# - Get all shortest path with  nx.all_shortest_paths(g, source, target)
# - calculate the sign based on the annotation of the interaction
# - If any 

# set(entry[1]['Effect'] for e in g.edge.items() for egroup in e[1].items() for entry in egroup[1].items())

# set((entry[1]['EffectMechanism'], entry[1]['Effect'])  for e in g.edge.items() for egroup in e[1].items() for entry in egroup[1].items())

# The problem here is that it is not trivial to assign positive and negative interactions.
# 
# Lets simplify and make a new network only having entities that are indicated to down or upregulate

# In[31]:

g_updown = g.copy()


# In[32]:

g_updown.edges()


# In[33]:

for e in set(g_updown.edges()):
    for key, item in g_updown.edge[e[0]][e[1]].items():
        e_dat = item
        if (e_dat['Effect'].startswith('down-regulates')):
            e_dat.update(({'mode':-1}))
        elif (e_dat['Effect'].startswith('up-regulates')):
            e_dat.update(({'mode':1}))
        else:
            e_dat.update(({'mode':0}))


# In[ ]:




# In[ ]:




# In[34]:

g_updown.edge[e[0]][e[1]]


# How many edges do we lose:

# In[35]:

dif = len(g_updown.edges())-len(g.edges())
print(dif)
print(abs(dif)/float(len(g.edges())))


# In[ ]:




# In[ ]:




# Which is quite low compared to the over 8000 edges in the network.

# Make a function that gives the shortest path given the network and a sign

# In[64]:

def get_shortest_dir_path(graph, origin, target, mode, maxlen=6, target_residues=None, phospho_dir=1):

    cur_len = nx.shortest_path_length(graph, origin, target)
    
    if cur_len == 0:
        return [origin]
    while cur_len <= maxlen:
        pathgen = nx.all_simple_paths(graph, origin, target, cutoff=cur_len)
        
        prev_paths = set()
        for path in pathgen:
            path = tuple(path)
            if path not in prev_paths:
                curmode = get_overal_mode(graph, path, target_residues=target_residues, phospho_dir=phospho_dir)
                if (curmode == mode) | (mode == 0):
                    return path
                else:
                    prev_paths.add(path)
            
        cur_len += 1
            
    raise(nx.NetworkXNoPath)
    
def get_overal_mode(graph, nodes, target_residues=None, phospho_dir=1):
    """
    """
    if target_residues is None:
        target_residues = []
    modes = list()
    for i in range(len(nodes)-1):
        
        idx, edges = zip(*graph.edge[nodes[i]][nodes[i+1]].items())
        if (i == len(nodes)-2):
            ## for the last node, check if it is a phosphorylation and if any of the target sites is in there
            ## if it is, check if it is 'dephosphorylation' or 'phosphorylation' and set sign accordingly.
            ## default back to 'mode'
            curmod = set()
            for e in edges:
                if e['MechanismResidues'] in target_residues:
                    if e['EffectMechanism'] == 'phosphorylation':
                        curmod.add(phospho_dir)
                    elif e['EffectMechanism'] == 'dephosphorylation':
                        curmod.add(phospho_dir*-1)
                    else:
                        curmod.add(0)
                else:
                    curmod.add(e['mode']*phospho_dir)        
        else:
            # take the annotation mode as mode
            curmod = set(e['mode'] for e in edges)
        if len(curmod) > 1:
            return 0
        else:
            modes.append(curmod.pop())
            
    overal_mode = np.prod(modes)
    return overal_mode


# In[ ]:




# In[37]:

#get_shortest_dir_path(g_updown, origin='P68400', target='Q9H2G2', mode=1)
g_updown.edge['P68400']['Q9H2G2']


# In[38]:

nx.shortest_path_length(g_updown, 'P68400', 'Q9H2G2')
pathgen = nx.all_simple_paths(g_updown, 'P68400', 'Q9H2G2', cutoff=1)
for p in pathgen:
    print(get_overal_mode(g_updown, p))


# In[76]:

def calc_disttable_signed(g, antibody, overexpression, direction, antibody_dict, overexp_dict):
    col = antibody
    row = overexpression

    pair_tab = np.empty((len(row), len(col)))
    pair_tab[:] = np.NAN
    for i, over in enumerate(row):
        for j, ab in enumerate(col):
            sp = list()
            for ab_prot in antibody_dict[ab]:
                ov_prot = overexp_dict[over]

                try:
                    sp.append(len(get_shortest_dir_path(g, ov_prot, ab_prot, direction))-1)
                except nx.NetworkXNoPath:
                    1
            if len(sp) > 0:
                pair_tab[i,j] = np.min(sp)
            else:
                pair_tab[i,j] = np.NaN
    pair_dat = pd.DataFrame(pair_tab,columns=col, index=row)
    pair_dat = pair_dat.stack(dropna=False)
    pair_dat.index = pair_dat.index.rename(['marker','target'])
    return pair_dat
    


pair_sp_pos = calc_disttable_signed(g=g_updown,
                                    antibody=antibody_dict.keys(), 
                      overexpression=overexp_tab['gene_xkl'],
                      direction=1,
                      antibody_dict=antibody_dict,
                      overexp_dict=overexp_dict
                     )
pair_sp_neg = calc_disttable_signed(g=g_updown, antibody=antibody_dict.keys(), 
                      overexpression=overexp_tab['gene_xkl'],
                      direction=-1,
                      antibody_dict=antibody_dict,
                      overexp_dict=overexp_dict
                     )




# In[77]:

def calc_disttable_signed_v2(g, antibody, overexpression, direction, antibody_dict_ext, overexp_dict):
    col = antibody
    row = overexpression

    pair_tab = np.empty((len(row), len(col)))
    pair_tab[:] = np.NAN
    for i, over in enumerate(row):
        for j, ab in enumerate(col):
            sp = list()
            for ab_prot, target_resid, phospho_dir in antibody_dict_ext[ab]:
                ov_prot = overexp_dict[over]
                

                try:
                    sp.append(len(get_shortest_dir_path(g, ov_prot, ab_prot, direction,target_residues=target_resid,
                                                        phospho_dir=phospho_dir))-1)
                except nx.NetworkXNoPath:
                    1
            if len(sp) > 0:
                pair_tab[i,j] = np.min(sp)
            else:
                pair_tab[i,j] = np.NaN
    pair_dat = pd.DataFrame(pair_tab,columns=col, index=row)
    pair_dat = pair_dat.stack(dropna=False)
    pair_dat.index = pair_dat.index.rename(['marker','target'])
    return pair_dat
    


pair_dat_pos_2 = calc_disttable_signed_v2(g=g_updown,
                                    antibody=antibody_dict_extended.keys(), 
                      overexpression=overexp_tab['gene_xkl'],
                      direction=1,
                      antibody_dict_ext=antibody_dict_extended,
                      overexp_dict=overexp_dict
                     )
pair_dat_neg_2 = calc_disttable_signed_v2(g=g_updown, antibody=antibody_dict_extended.keys(), 
                      overexpression=overexp_tab['gene_xkl'],
                      direction=-1,
                      antibody_dict_ext=antibody_dict_extended,
                      overexp_dict=overexp_dict
                     )


# In[73]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[78]:

bins = np.arange(int(np.ceil(pair_dat_pos.max()))+1)-0.5
pair_dat_pos.hist(bins=bins)
pair_dat_pos.describe()


# In[79]:

bins = np.arange(int(np.ceil(pair_dat_pos_2.max()))+1)-0.5
pair_dat_pos_2.hist(bins=bins)
pair_dat_pos_2.describe()


# In[80]:

bins = np.arange(int(np.ceil(pair_dat_neg_2.max()))+1)-0.5
pair_dat_neg_2.hist(bins=bins)
pair_dat_neg_2.describe()


# In[81]:

bins = np.arange(int(np.ceil(pair_dat_neg.max()))+1)-0.5
pair_dat_neg.hist(bins=bins)
pair_dat_neg.describe()


# In[ ]:




# ## Look at the minium number of measured nodes at the shortest path

# In[36]:

def get_min_measured_nodes_shortest_path(g, source, target, measured_nodes):
    """
    
    """

    
    try:
        sp = nx.shortest_path(g, source, target)
    except nx.NetworkXNoPath:
        return 10

    sp = nx.all_shortest_paths(g, source, target)
        
    n_measured_between = max(sum(1 for n in nl if n in measured_nodes) for nl in sp)
    
    return n_measured_between

col = antibody_dict.keys()
row = overexp_tab['gene_xkl']

pair_sp_nbetween = np.empty((len(row), len(col)))
pair_sp_nbetween[:] = np.NAN
for i, over in enumerate(row):
    for j, ab in enumerate(col):
        sp = list()
        for ab_prot in antibody_dict[ab]:
            ov_prot = overexp_dict[over]
            try:
                sp.append(get_min_measured_nodes_shortest_path(g, ov_prot, ab_prot, set(xkl_nodes)))
            except nx.NetworkXNoPath:
                pass
        if len(sp) > 0:
            pair_sp_nbetween[i,j] = np.min(sp)-1
        else:
            pair_sp_nbetween[i,j] = 10

pair_dat_nbetween = pd.DataFrame(pair_sp_nbetween,columns=col, index=row)
pair_dat_nbetween = pair_dat_nbetween.stack()
pair_dat_nbetween.index = pair_dat_nbetween.index.rename(['marker','target'])


# In[37]:

bins = np.arange(int(np.ceil(pair_dat_nbetween.max()))+1)-0.5
pair_dat_nbetween.hist(bins=bins)
pair_dat_nbetween.describe()


# ## Check consistency with bp-R2 analysis

# In[85]:

bp_file = '/mnt/imls-bod/Xiao-Kang/EGF transfection/plots/nbin10_2.5perc_bpr2_median_25_final/t_bindat'


# In[86]:

bin_dat = pd.read_pickle(bp_file)


# In[87]:

print(bin_dat.index.get_level_values('target').unique())
print(bin_dat.index.get_level_values('marker').unique())

#pair_dat = pair_dat_nbetween


# In[88]:

sig_dat = bin_dat.loc[bin_dat['bin_dat_sigfil_any']].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']
plt.figure()
bins = np.arange(int(np.ceil(pair_dat.max()))+1)-0.5
plt.hist(pair_dat, bins=bins,normed=True)

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))


sp_dist_hits = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits, bins=bins,normed=True, alpha=0.5)
print(sp_dist_hits.mean())
print(sp_dist_hits.median())
plt.figure()
plt.hist(sp_dist_hits, bins=bins, alpha=0.5)


# In[ ]:




# In[ ]:




# -> slight enrichment when considering all pairs significant in at least one timepoint

# In[75]:

sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 0) & bin_dat['bin_dat_sigfil'] & (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']
plt.figure()

plt.hist(pair_dat, bins=bins,normed=True)



sp_dist_hits_ss = sig_dat.apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits_ss, bins=bins,normed=True, alpha=0.5)
sig_dat_0 = sig_dat.copy()
print(sp_dist_hits_ss.mean())
print(sp_dist_hits_ss.median())
plt.figure()
plt.hist(sp_dist_hits_ss, bins=bins, alpha=0.5)


# -> stronger enrichment when considering  only significant pairs at steady state

# In[76]:


sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 5) & bin_dat['bin_dat_sigfil']].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']

plt.figure()
plt.hist(pair_dat, bins=bins,normed=True)


sig_dat_0['temp'] = 1
sig_dat = pd.merge(sig_dat, sig_dat_0, how='outer')
sig_dat = sig_dat[sig_dat['temp'] != 1]

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sp_dist_hits = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits, bins=bins,normed=True, alpha=0.5)
print(sp_dist_hits.mean())
print(sp_dist_hits.median())
plt.figure()
plt.hist(sp_dist_hits, bins=bins, alpha=0.5)


# -> weaker enrichement when considering only after 5 min

# In[77]:


sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 15) & bin_dat['bin_dat_sigfil']].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']

plt.figure()

plt.hist(pair_dat, bins=bins,normed=True)


sig_dat_0['temp'] = 1
sig_dat = pd.merge(sig_dat, sig_dat_0, how='outer')
sig_dat = sig_dat[sig_dat['temp'] != 1]

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sp_dist_hits = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits, bins=bins,normed=True, alpha=0.5)

print(sp_dist_hits.mean())
print(sp_dist_hits.median())
plt.figure()
plt.hist(sp_dist_hits, bins=bins, alpha=0.5)


# -> even weaker after only 15

# In[78]:


sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 30) & bin_dat['bin_dat_sigfil']].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']

plt.figure()
plt.hist(pair_dat, bins=bins,normed=True)


sig_dat_0['temp'] = 1
sig_dat = pd.merge(sig_dat, sig_dat_0, how='outer')
sig_dat = sig_dat[sig_dat['temp'] != 1]

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sp_dist_hits = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits, bins=bins,normed=True, alpha=0.5)
print(sp_dist_hits.mean())
print(sp_dist_hits.median())
plt.figure()
plt.hist(sp_dist_hits, bins=bins, alpha=0.5)


# -> even less direct after 30

# In[79]:


sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 60) & bin_dat['bin_dat_sigfil']].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']

plt.figure()
plt.hist(pair_dat, bins=bins,normed=True)


sig_dat_0['temp'] = 1
sig_dat = pd.merge(sig_dat, sig_dat_0, how='outer')
sig_dat = sig_dat[sig_dat['temp'] != 1]

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sp_dist_hits = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits, bins=bins,normed=True, alpha=0.5)
print(sp_dist_hits.mean())
print(sp_dist_hits.median())
plt.figure()
plt.hist(sp_dist_hits, bins=bins, alpha=0.5)


# -> enriched for distant after 60 min

# In[80]:

fil2 = (bin_dat.index.get_level_values('timepoint') == 0) & (bin_dat['bin_dat_sigfil'] == False)

sig_dat = bin_dat.loc[fil2 & (bin_dat.index.get_level_values('timepoint') == 5) & bin_dat['bin_dat_sigfil']].reset_index(drop=False)
sig_dat.columns


# In[ ]:




# In[ ]:




# ### Look at the hists at different distances

# Add the distance to the bindat

# In[84]:

sig_dat = bin_dat.loc[bin_dat['bin_dat_sigfil_any']].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']
plt.figure()
bins = np.arange(int(np.ceil(pair_dat.max()))+1)-0.5
plt.hist(pair_dat, bins=bins,normed=True)

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = sig_dat.apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
sig_dat = sig_dat.sort_values('dist', ascending=False)
plt.hist(sig_dat['dist'], bins=bins,normed=True, alpha=0.5)


# In[83]:

sig_dat = bin_dat.loc[ (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target','timepoint']].drop_duplicates(['marker', 'target','timepoint'])
sig_dat.columns = ['marker', 'target','timepoint']
plt.figure()
bins = np.arange(int(np.ceil(pair_dat.max()))+1)-0.5
plt.hist(pair_dat, bins=bins,normed=True)

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = sig_dat.apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
sig_dat = sig_dat.sort_values(['dist','marker','target','timepoint'], ascending=False)
plt.hist(sig_dat['dist'], bins=bins,normed=True, alpha=0.5)

sig_dat = sig_dat.sort_values(['dist','marker','target','timepoint'], ascending=True)
sig_dat = sig_dat.drop_duplicates(['marker', 'target'])
sig_dat = sig_dat.sort_values(['dist','marker','target','timepoint'], ascending=False)
HTML(sig_dat.to_html())


# Use the in 'measured nodes in between' distance
# 

# In[84]:

sig_dat = bin_dat.loc[ (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target','timepoint']].drop_duplicates(['marker', 'target','timepoint'])
sig_dat.columns = ['marker', 'target','timepoint']
plt.figure()
bins = np.arange(int(np.ceil(pair_dat_nbetween.max()))+1)-0.5
plt.hist(pair_dat_nbetween, bins=bins,normed=True)

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = sig_dat.apply(lambda x: pair_dat_nbetween.loc[x['marker'], x['target']],axis=1)
sig_dat = sig_dat.sort_values(['dist','marker','target','timepoint'], ascending=False)
plt.hist(sig_dat['dist'], bins=bins,normed=True, alpha=0.5)
sig_dat.loc[sig_dat['timepoint'] == 0]


# In[85]:

sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 0) & bin_dat['bin_dat_sigfil'] & (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat.drop_duplicates(['marker', 'target'])
sig_dat =sig_dat[['marker', 'target','stats']]
sig_dat = sig_dat.drop_duplicates(['marker', 'target'])
fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))


sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = [pair_dat.loc[m,t] for m, t in zip(sig_dat['marker'], sig_dat['target'])]


# In[86]:

sig_dat = sig_dat.sort_values('dist', ascending=False)
sig_dat.head(50)


# In[ ]:




# ## investigate how many of the 'close' relationships where recovered

# In[87]:

sig_dat = bin_dat.reset_index(drop=False)
#sig_dat =sig_dat.drop_duplicates(['marker', 'target'])
fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))


sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = [pair_dat.loc[m,t] for m, t in zip(sig_dat['marker'], sig_dat['target'])]


# In[88]:

fil = (sig_dat['dist'] <2)& (sig_dat['timepoint'] == 5)
#
sig_dat.loc[fil]['bin_dat_sigfil_any'].mean()


fil = sig_dat['timepoint'] == 0
[sig_dat.loc[fil & (sig_dat['dist'] == i)]['bin_dat_sigfil'].mean() for i in range(10)]


# In[ ]:




# In[89]:

fil = (sig_dat['bin_dat_sigfil_any'])
#

print((sig_dat.loc[fil]['dist'] >1).mean())
print((sig_dat.loc[fil]['dist'] <2).mean())
print((sig_dat.loc[fil]['dist'] >5).mean())



# In[ ]:




# # With the signed directed network

# In[ ]:




# In[89]:

sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 0) & bin_dat['bin_dat_sigfil'] & (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat.groupby(['marker', 'target']).apply(lambda x: np.sign(x[('stats', 'corr_spearman_bin')].mean()))
sig_dat.name = 'sign'
sig_dat = sig_dat.reset_index(drop=False)


# In[220]:

plt.hist(pair_dat_pos.append(pair_dat_neg).dropna(), bins=bins,normed=True)


fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))

def _get_data(x):
    if x['sign'] > 0:
        return(pair_dat_pos.loc[x['marker'], x['target']])
    else:
        return(pair_dat_neg.loc[x['marker'], x['target']])

sp_dist_hits_sign = sig_dat.loc[fil].apply(lambda x: _get_data(x),axis=1)

sp_dist_hits_nonsign = sig_dat.apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)

#plt.hist(sp_dist_hits_nonsign.dropna(), bins=bins,normed=True, alpha=0.5)
plt.hist(sp_dist_hits_sign.dropna(), bins=bins,normed=True, alpha=0.5)

sig_dat_0 = sig_dat.copy()
print(sp_dist_hits_sign.mean())
print(sp_dist_hits_sign.median())
plt.figure()
plt.hist(sp_dist_hits_sign, bins=bins, alpha=0.5)


# In[221]:

from IPython.display import display, HTML

overexp_tab = overexp_tab.set_index('gene_xkl',drop=False)
full_dat = pd.concat([sig_dat.loc[fil].reset_index(drop=True),
                      pd.DataFrame({'dist nonsign':sp_dist_hits_nonsign}).reset_index(drop=True),
                      pd.DataFrame({'dist sign':sp_dist_hits_sign}).reset_index(drop=True),
                      pd.DataFrame({'protein': overexp_tab.loc[sig_dat.loc[fil,'marker'],'Overexpressed proteins'].reset_index(drop=True)})], axis=1)

HTML(full_dat.sort_values('dist sign').to_html())


# In[ ]:




# In[223]:

fil = (pair_dat.index.get_level_values('marker').isin(not_used_effects['marker']) &
       pair_dat.index.get_level_values('target').isin(not_used_effects['target'])) == False 


# In[225]:

pair_dat_pos_neg = pair_dat_pos_2.append(pair_dat_neg_2).dropna()

fil = (pair_dat_pos_neg.index.get_level_values('marker').isin(not_used_effects['marker']) &
       pair_dat_pos_neg.index.get_level_values('target').isin(not_used_effects['target'])) == False 

plt.hist(pair_dat_pos_neg.loc[fil], bins=bins,normed=True)


fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys())) & (sig_dat['marker'].isin(not_used_effects['marker']) &
       sig_dat['target'].isin(not_used_effects['target'])) == False 


def _get_data(x):
    if x['sign'] > 0:
        return(pair_dat_pos_2.loc[x['marker'], x['target']])
    else:
        return(pair_dat_neg_2.loc[x['marker'], x['target']])

sp_dist_hits_sign = sig_dat.loc[fil].apply(lambda x: _get_data(x),axis=1)

sp_dist_hits_nonsign = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)

#plt.hist(sp_dist_hits_nonsign.dropna(), bins=bins,normed=True, alpha=0.5)
plt.hist(sp_dist_hits_sign.dropna(), bins=bins,normed=True, alpha=0.5)

sig_dat_0 = sig_dat.copy()
print(sp_dist_hits_sign.mean())
print(sp_dist_hits_sign.median())
plt.figure()
plt.hist(sp_dist_hits_sign, bins=bins, alpha=0.5)


# In[226]:

from IPython.display import display, HTML

overexp_tab = overexp_tab.set_index('gene_xkl',drop=False)
full_dat = pd.concat([sig_dat.loc[fil].reset_index(drop=True),
                      pd.DataFrame({'dist nonsign':sp_dist_hits_nonsign}).reset_index(drop=True),
                      pd.DataFrame({'dist sign':sp_dist_hits_sign}).reset_index(drop=True),
                      pd.DataFrame({'protein': overexp_tab.loc[sig_dat.loc[fil,'marker'],'Overexpressed proteins'].reset_index(drop=True)})], axis=1)

HTML(full_dat.sort_values('dist sign').to_html())


# In[ ]:




# In[ ]:




# In[228]:

def get_path(origin, target, sign):
    origin_prot = overexp_dict[origin]
    paths = list()
    
    for target_prot, target_sites, phospho_dir in antibody_dict_extended[target]:
        try:
            paths.append(get_shortest_dir_path(g_updown, origin_prot, target_prot, sign,
                                   target_residues=target_sites,phospho_dir=phospho_dir))
        except:
            pass
    return paths


full_dat['path']=full_dat.apply(lambda x: ' _ '.join([' '.join(path) for path in get_path(x['marker'], x['target'],x['sign'])]), axis=1)


# In[229]:

HTML(full_dat.sort_values('dist sign').to_html())

full_dat.sort_values('dist sign').to_csv(os.path.join(outfolder, '20161020_shortest_paths_revised.csv'))


# Find for each antiboy if its phosphosite is activing or inhibitory as in SIGNOR

# ## Do the readout comparsion based on the networks
# 

# In[230]:

# from 20160312_readout_comparison, accessed 20161018

fn_hits = os.path.join('/home/vitoz/imls-bod/Xiao-Kang/EGF transfection/benchmark',
                       '20160429_readout_comparison_vsemptygfp_woCC.csv')


# In[231]:

hittab = pd.read_csv(fn_hits)
hittab['marker'] =hittab['all_hits'].map(lambda x: x.split('_')[0])
hittab['target'] =hittab['all_hits'].map(lambda x: x.split('_')[1])


# In[ ]:




# In[232]:

fil = (pair_dat.index.get_level_values('marker').isin(not_used_effects['marker']) &
       pair_dat.index.get_level_values('target').isin(not_used_effects['target'])) == False 

pair_dat_wo = pair_dat.loc[fil,:]


# In[ ]:




# In[ ]:




# In[233]:

plt.hist(pair_dat_wo, bins=bins,normed=True)
plt.figure()
plt.hist(pair_dat_wo, bins=bins,normed=True)
print(np.mean(pair_dat_wo))
 
hittab_fil = (hittab['marker'].isin(not_used_effects['marker']) &
       hittab['target'].isin(not_used_effects['target'])) == False 



#plt.hist(sp_dist_hits_nonsign.dropna(), bins=bins,normed=True, alpha=0.5)
cols= ['is_bp', 'is_dremi','is_sp']
for c in cols:
    fil = (hittab[c] == True) & hittab_fil
    tdist = pair_dat_wo.loc[[(row['marker'],row['target']) for idx, row in hittab.loc[fil].iterrows()]]
    plt.hist(tdist, bins=bins,normed=True, alpha=0.5)
    print([c, tdist.mean()])


for c in cols:
    plt.figure()
    fil = (hittab[c] == True) & hittab_fil
    tdist = pair_dat_wo.loc[[(row['marker'],row['target']) for idx, row in hittab.loc[fil].iterrows()]]

    plt.hist(tdist, bins=bins, alpha=0.5)
    plt.title(c)
    
plt.figure()
for c in cols:
    fil = (hittab[c] == True) & hittab_fil
    tdist = pair_dat_wo.loc[[(row['marker'],row['target']) for idx, row in hittab.loc[fil].iterrows()]]
    plt.hist(tdist, bins=bins, alpha=0.5)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




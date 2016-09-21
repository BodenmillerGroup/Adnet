
# coding: utf-8

# ## Aim
# 
# Compare the coverage of the Signor network with the markers measured and the overexpressions performed

# In[1]:

import networkx as nx
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from IPython.display import display, HTML
get_ipython().magic('matplotlib inline')


# In[2]:

# The complete signor network
signor_fn = '/mnt/imls-bod/Xiao-Kang/EGF transfection/20160110_signor_all_data.tsv'


# In[3]:

# The antibody names in uniprot
antibody_fn = '/mnt/imls-bod/Xiao-Kang/EGF transfection/201602_antibody_uniprot.csv'


# In[4]:

# The overexpressions in uniprot
overexp_fn = '/mnt/imls-bod/Xiao-Kang/EGF transfection/201602_gene_uniprot.csv'


# In[5]:

# Effect mechanism filters for signor
not_used_effects = ['']


# Read the network

# In[6]:

signor_tab = pd.read_csv(signor_fn, sep='\t')
print(signor_tab.columns)
signor_tab



# In[7]:

# Unique  interaction types
signor_tab[u'EffectMechanism'].unique()



# In[8]:

g= nx.from_pandas_dataframe(signor_tab,'IdA','IdB',edge_attr=[u'Effect', u'EffectMechanism',
       u'MechanismResidues', u'MechanismSequences',u'Direct'], create_using=nx.DiGraph())


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


# ## Make a minimum spanning tree network containing all the nodes measured by antibodies or overexpressions

# In[12]:

xkl_nodes = np.unique(overexp_tab['UniProt'].dropna().tolist() + antibody_tab['UniProt'].dropna().tolist())



# In[13]:

all_nodes = g.nodes()
signor_xkl_nodes = [n for n in all_nodes if n in xkl_nodes]
len(signor_xkl_nodes) == len(xkl_nodes)


# -> All XKL nodes are actually in Signor! Cool!

# Prune the graph to contain only nodes contained in all shortest paths between all the xkl_nodes

# In[14]:

def shortest_path_nodes(graph, source, target, max_cutoff = 10, all_cutoff=5):
    """
    Returns all nodes between the shortest path between source and target
    max_cutoff: Return empty set if the shortest path exceeds this cutoff.
    all_cutoff: Consider all possible shortest paths, for paths shorter than this cuttoff.
    """
    
    assert(max_cutoff >= all_cutoff)
    
    try:
        sp = nx.shortest_path(g, source, target)
    except nx.NetworkXNoPath:
        sp =set()
    sp_len = len(sp)
    if 0 < sp_len <= all_cutoff:
        sp = nx.all_shortest_paths(g, source, target)
        sp = set(n for nl in sp for n in nl)
    if sp_len <= max_cutoff:
        return(set(sp))
    
    else:
        return(set())


# In[15]:

#nx.single_source_shortest_path(g,signor_xkl_nodes[0])
#list(nx.all_simple_paths(g,signor_xkl_nodes[0],signor_xkl_nodes[1],5))


n=len(signor_xkl_nodes)
all_xkl_spnodes = set()
for i in range(n):
    for j in range(n):
        nodes = shortest_path_nodes(g, signor_xkl_nodes[i], signor_xkl_nodes[j],
                                           max_cutoff = 10, all_cutoff=5)
        all_xkl_spnodes.update(nodes)

g_xkl = g.copy()
for n in g_xkl.nodes():
    if n not in all_xkl_spnodes:
        g_xkl.remove_node(n)
        
print(len(all_xkl_spnodes))
#print(all_xkl_spnodes)


# In[ ]:




# In[16]:

plt.figure(figsize=(20,20))
nx.draw_spring(g_xkl)


# The graph is really to big for visualization
# 
# Sanity check: check for hubs -> might be artificial shortcuts

# In[17]:

centr = pd.DataFrame.from_dict([nx.degree_centrality(g_xkl)])
centr = centr.swapaxes(0,1)

centr = centr.reset_index(drop=False)
centr.columns = ['id', 'centrality']


centr['name'] = [signor_dict[id] for id in centr['id']]
centr.sort_values('centrality', ascending=False, inplace=True)


# In[18]:

centr.head(10)


# How does this look for direct interactions only?

# In[19]:

g_dir = g.copy()


# In[20]:

for e in g_dir.edges():
    e_dat = g_dir.get_edge_data( e[0], e[1])
    if e_dat['Direct'] != 'YES':
        g_dir.remove_edge( e[0], e[1])



# In[21]:

print(len(g_dir.edges()))

print(len(g.edges()))


# In[22]:

n=len(signor_xkl_nodes)
all_dir_xkl_spnodes = set()
for i in range(n):
    for j in range(n):
        nodes = shortest_path_nodes(g_dir, signor_xkl_nodes[i], signor_xkl_nodes[j],
                                           max_cutoff = 10, all_cutoff=5)
        all_dir_xkl_spnodes.update(nodes)

g_dir_xkl = g_dir.copy()
for n in g_dir_xkl.nodes():
    if n not in all_dir_xkl_spnodes:
        g_dir_xkl.remove_node(n)
        
print(len(all_dir_xkl_spnodes))


# -> In terms of total nodes considering only direct edges does not change anything...

# ## Check the average shortest path length between overexpressions and phosphoproteins

# In[ ]:




# In[ ]:




# Calculate pairwise shortest path distance, taking the minimum distanc if the antibody is not specific

# In[23]:

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
 



# In[24]:

pair_dat = pd.DataFrame(pair_sp,columns=col, index=row)
pair_dat = pair_dat.stack()
pair_dat.index = pair_dat.index.rename(['marker','target'])


# Look at the distribution

# In[25]:

bins = np.arange(int(np.ceil(pair_dat.max()))+1)-0.5
pair_dat.hist(bins=bins)
pair_dat.describe()


# 
# Look at the minium number of measured nodes at the shortest path

# In[26]:

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


# In[27]:

bins = np.arange(int(np.ceil(pair_dat_nbetween.max()))+1)-0.5
pair_dat_nbetween.hist(bins=bins)
pair_dat_nbetween.describe()


# In[28]:

(sum(1 for n in range(10) if n >10))


# In[53]:

# From comparison analysis
# BPR2 only

x = [(s.split('_')[0],s.split('_')[1]) for s in {'AKT1_Beta-catenin_0',
 'AKT1_p-4EBP1_0',
 'AKT1_p-GSK3-Beta_0',
 'AKT1_p-SMAD1-5_0',
 'AKT1_p-SMAD2-3_0',
 'BRAF_p-MKK3-6_0',
 'GSK3B_Beta-catenin_0',
 'GSK3B_p-FAK_0',
 'GSK3B_p-MEK1-2_0',
 'GSK3B_p-PDPK1_0',
 'GSK3B_p-STAT5_0',
 'HRAS_p-SMAD1-5_0',
 'MAP2K1_p-4EBP1_0',
 'MAP2K1_p-BTK_0',
 'MAP2K1_p-GSK3-Beta_0',
 'MAP2K6_p-GSK3-Beta_0',
 'PDPK1_p-MAPKAPK2_0',
 'PDPK1_p-S6_0',
 'PIK3CA_p-SMAD1-5_0',
 'PIK3CA_p-p38_0',
 'PTPN11_p-JNK_0',
 'SRC_p-JNK_0'}


]


print(pair_dat.loc[x])
print(np.nanmean(pair_dat.loc[x]))


# In[54]:

# bpR2 but not Dremi
x = [(s.split('_')[0],s.split('_')[1]) for s in {'AKT1_Beta-catenin_0',
 'AKT1_p-4EBP1_0',
 'AKT1_p-AKT_0',
 'AKT1_p-GSK3-Beta_0',
 'AKT1_p-MKK3-6_0',
 'AKT1_p-SMAD1-5_0',
 'AKT1_p-SMAD2-3_0',
 'AKT1_p-p70S6K_0',
 'BRAF_p-MKK3-6_0',
 'CRAF_p-p90RSK_0',
 'GSK3B_Beta-catenin_0',
 'GSK3B_p-FAK_0',
 'GSK3B_p-MEK1-2_0',
 'GSK3B_p-PDPK1_0',
 'GSK3B_p-SHP2_0',
 'GSK3B_p-STAT5_0',
 'GSK3B_p-p70S6K_0',
 'HRAS_p-JNK_0',
 'HRAS_p-SMAD1-5_0',
 'HRAS_p-SMAD2-3_0',
 'MAP2K1_p-4EBP1_0',
 'MAP2K1_p-BTK_0',
 'MAP2K1_p-GSK3-Beta_0',
 'MAP2K6_p-GSK3-Beta_0',
 'MAP2K6_p-STAT1_0',
 'MAP2K7_p-SMAD1-5_0',
 'MAP3K5_p-GSK3-Beta_0',
 'MAP3K5_p-HH3_0',
 'MAP3K5_p-MKK3_0',
 'MAP3K5_p-PDPK1_0',
 'MAP3K5_p-SMAD2-3_0',
 'MAP3K5_p-STAT1_0',
 'MAPK1_p-ERK1-2_0',
 'PDPK1_p-BTK_0',
 'PDPK1_p-MAPKAPK2_0',
 'PDPK1_p-S6_0',
 'PIK3CA_p-JNK_0',
 'PIK3CA_p-MKK3-6_0',
 'PIK3CA_p-SMAD1-5_0',
 'PIK3CA_p-p38_0',
 'PTPN11_p-JNK_0',
 'RPS6KA1_p-PDPK1_0',
 'SRC_p-JNK_0'}
]


print(pair_dat.loc[x])
print(np.nanmean(pair_dat.loc[x]))
pair_dat.loc[x].hist()


# In[55]:

# dremi only
x = [(s.split('_')[0],s.split('_')[1]) for s in {'AKT1_p-RB_0',
 'AKT1_p-S6_0',
 'BRAF_p-S6_0',
 'GFP-FLAG-1_cyclin B1_0',
 'GFP-FLAG-1_p-RB_0',
 'GFP-FLAG-2_cyclin B1_0',
 'GFP-FLAG-2_p-RB_0',
 'HRAS_p-RB_0',
 'MAP2K1_cyclin B1_0',
 'MAP2K1_p-RB_0',
 'MAP2K1_p-S6_0',
 'MAP2K6_cyclin B1_0',
 'MAP2K6_p-RB_0',
 'MAP3K5_Beta-catenin_0',
 'MAP3K5_p-RB_0',
 'MAP8_cyclin B1_0',
 'MAP8_p-PLCg2_0',
 'MAPK14_cyclin B1_0',
 'MAPK14_p-RB_0',
 'PDPK1_p-AKT_0',
 'PTPN11_p-MARCKS_0',
 'PTPN11_p-RB_0',
 'RPS6KA1_p-RB_0',
 'RPS6KB1_cyclin B1_0',
 'RPS6KB1_p-RB_0',
 'SRC_Beta-catenin_0'}
]


print(pair_dat.loc[x])
print(np.nanmean(pair_dat.loc[x]))
pair_dat.loc[x].hist()


# ## Check consistency with bp-R2 analysis

# In[32]:

bp_file = '/mnt/imls-bod/Xiao-Kang/EGF transfection/plots/nbin10_2.5perc_bpr2_median_25_final/t_bindat'


# In[33]:

bin_dat = pd.read_pickle(bp_file)


# In[34]:

print(bin_dat.index.get_level_values('target').unique())
print(bin_dat.index.get_level_values('marker').unique())

#pair_dat = pair_dat_nbetween


# In[35]:

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




# -> slight enrichment when considering all pairs significant in at least one timepoint

# In[36]:

sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 0) & bin_dat['bin_dat_sigfil'] & (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat[['marker', 'target']].drop_duplicates(['marker', 'target'])
sig_dat.columns = ['marker', 'target']
plt.figure()

plt.hist(pair_dat, bins=bins,normed=True)

fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))


sp_dist_hits = sig_dat.loc[fil].apply(lambda x: pair_dat.loc[x['marker'], x['target']],axis=1)
plt.hist(sp_dist_hits, bins=bins,normed=True, alpha=0.5)
sig_dat_0 = sig_dat.copy()
print(sp_dist_hits.mean())
print(sp_dist_hits.median())
plt.figure()
plt.hist(sp_dist_hits, bins=bins, alpha=0.5)


# -> stronger enrichment when considering  only significant pairs at steady state

# In[37]:


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

# In[38]:


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

# In[39]:


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

# In[40]:


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

# In[41]:

fil2 = (bin_dat.index.get_level_values('timepoint') == 0) & (bin_dat['bin_dat_sigfil'] == False)

sig_dat = bin_dat.loc[fil2 & (bin_dat.index.get_level_values('timepoint') == 5) & bin_dat['bin_dat_sigfil']].reset_index(drop=False)
sig_dat.columns


# In[ ]:




# In[42]:

fil2& fil


# ### Look at the hists at different distances

# Add the distance to the bindat

# In[43]:

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


# In[44]:

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

# In[45]:

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


# In[46]:

sig_dat = bin_dat.loc[(bin_dat.index.get_level_values('timepoint') == 0) & bin_dat['bin_dat_sigfil'] & (bin_dat['bin_dat_sigfil'])].reset_index(drop=False)
sig_dat =sig_dat.drop_duplicates(['marker', 'target'])
sig_dat =sig_dat[['marker', 'target','stats']]
sig_dat = sig_dat.drop_duplicates(['marker', 'target'])
fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))


sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = [pair_dat.loc[m,t] for m, t in zip(sig_dat['marker'], sig_dat['target'])]


# In[47]:

sig_dat = sig_dat.sort_values('dist', ascending=False)
sig_dat.head(50)


# In[ ]:




# ## investigate how many of the 'close' relationships where recovered

# In[48]:

sig_dat = bin_dat.reset_index(drop=False)
#sig_dat =sig_dat.drop_duplicates(['marker', 'target'])
fil = (sig_dat['marker'].apply(lambda x: x in overexp_dict.keys())) & (
    sig_dat['target'].apply(lambda x: x in antibody_dict.keys()))


sig_dat = sig_dat.loc[fil]
sig_dat['dist'] = [pair_dat.loc[m,t] for m, t in zip(sig_dat['marker'], sig_dat['target'])]


# In[49]:

fil = (sig_dat['dist'] <2)& (sig_dat['timepoint'] == 5)
#
sig_dat.loc[fil]['bin_dat_sigfil_any'].mean()


fil = sig_dat['timepoint'] == 0
[sig_dat.loc[fil & (sig_dat['dist'] == i)]['bin_dat_sigfil'].mean() for i in range(10)]


# In[50]:

fil = sig_dat['timepoint'] == 0
sig_dat.loc[fil & (sig_dat['dist'] == 0) & (sig_dat['bin_dat_sigfil'] == False)]
#sig_dat.loc[fil & (sig_dat['dist'] == 0) & (sig_dat['bin_dat_sigfil'] == True)]


# In[51]:

fil = (sig_dat['bin_dat_sigfil_any'])
#

print((sig_dat.loc[fil]['dist'] >1).mean())
print((sig_dat.loc[fil]['dist'] <2).mean())
print((sig_dat.loc[fil]['dist'] >5).mean())



# In[ ]:




# In[12]:

np.round(2**np.linspace(2,7,10))


# In[ ]:




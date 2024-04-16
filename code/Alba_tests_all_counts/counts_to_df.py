# import libraries
import pandas as pd
import scipy.stats
import statsmodels.stats.multitest
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# project specific libs
import os
import matplotlib.pyplot as plt

# set working directory file path
path = '/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/code/'

# read in initial df
df_raw = pd.read_csv(path + 'All_counts.tsv', sep='\t')#, index_col=0)

# get rid of all vars that are irrelevant, also optionally do log transform
# get hard copy
df_counts = df_raw.copy()

# take log
# df_counts['logLive'] = np.log(df_counts['Live.cells.cm2'])
df_counts['Total'] = df_counts['Total.cells.cm2']

# remove unnecessary columns
df_counts = df_counts.drop(['Live_sorted_pct','Dead_sorted_pct','Total.cells.cm2','Dead.cells.cm2','Live.cells.cm2'],axis=1)

# grab subset dfs on site
df_cs = df_counts[df_counts['DeliveryMode'] == 'Csection']
df_vg = df_counts[df_counts['DeliveryMode'] == 'Vaginal']

# drop duplication?
df_cs = df_cs.drop(395) # Gauze duplication, PID 13706
df_cs = df_cs.drop(387) # Gauze duplication, PID 10805

site_to_ogdf = {
    'Csection': df_cs,
    'Vaginal': df_vg
}

site_to_newdf = {}
site_to_newdf_wide = {}

# for each site
for s in ['Csection', 'Vaginal']:
    # get relevant df
    df = site_to_ogdf[s]
    
    # drop DeliveryMode
    df = df.drop('DeliveryMode',axis=1)
    
    # get PIDs and Speciments
    # for each PID and Specimen pairing, if logLive is not a float, make it NaN
    PIDs = list(set(df['PID'].unique())) # 80 unique PIDs
    types = list(set(df['Specimen'].unique())) # 4 unique types

    for p in PIDs:
        for t in types:
            ll = df[(df['PID'] == p) & (df['Specimen'] == t)]['Total'].values # will be empty list if DNE
            if ll.size == 0:
                new_row = pd.DataFrame(data={'PID':p, 'Specimen': t, 'Total': np.nan}, index=[p]) 
                df = pd.concat([df, new_row])#, ignore_index=True)
    site_to_newdf[s] = df
    df.to_csv(path + s + '_long_df.csv')
    site_to_newdf_wide[s] = df.pivot(index='PID', columns='Specimen', values='Total')
    site_to_newdf_wide[s].to_csv(path + s + '_wide_df.csv')
    
df.head()
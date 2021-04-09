#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import numpy as np

# In[2]:

parser=argparse.ArgumentParser()
parser.add_argument('outfile')
parser.add_argument('rsemfiles',nargs='+')
opts = parser.parse_args()
files = opts.rsemfiles

gene_data_list = []
for file in files:
    print(file)
    data = pd.read_table(file,index_col=0)
    gene_data = data[['expected_count']]
    gene_data.columns = [file]
    gene_data_list.append(gene_data)
    
all_samples = pd.concat(gene_data_list,axis=1,join='inner')

# In[41]:


def normalize(samples_dataframe):
    samples_data = samples_dataframe.values
    N = len(samples_dataframe.columns)
    print(samples_dataframe.columns)
    ratios = np.zeros((N,))
    sample_ratios = {'sample':[],'ratio':[]}
    for u,sample in enumerate(samples_dataframe.columns):
        print(sample)
        g_ratio = np.zeros((N,))
        for v,_ in enumerate(samples_dataframe.columns):
            g_ratio_sample = []
            for j,gene in enumerate(samples_dataframe.index):
                if samples_data[j,v] != 0 and samples_data[j,u] != 0:
                    #ignore any genes with zero in either sample
                    g_ratio_sample.append(samples_data[j,u]/samples_data[j,v])
            g_ratio[v] = np.median(g_ratio_sample)
        ratios[u] = np.prod(g_ratio)**(1.0/N)
        sample_ratios['sample'].append(sample)
        sample_ratios['ratio'].append(ratios[u])
    print(pd.DataFrame(sample_ratios))
    return samples_dataframe/ratios,pd.DataFrame(sample_ratios)


normalized_gene_expression, norm_ratios =normalize(all_samples)

normalized_gene_expression.to_csv(opts.outfile+'.normalized.expr')

norm_ratios.to_csv(opts.outfile+'.ratios')

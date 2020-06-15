#!/usr/bin/env python
# coding: utf-8

# # Statistical tests for AI and AI within ASE

# Import libraries and data

# In[1]:


import pandas as pd
import numpy as np
import math
from scipy import stats
import scipy as sc
import re #regex


# Let's work out the file containing the SNPs and the allele frequencies created in the script "ap bias 2maps"

# In[2]:


df=pd.read_csv('SNPs2maps_GF3KF6.csv')


# In[3]:


df.head()


# Rename the column names to addapt to the script

# In[4]:


# get the column names
for col in df.columns: 
    print(col) 


# In[5]:


df.rename(columns={'Gene.refGene_x': 'Gene.refGene',
                'Func.refGene_x': 'Func.refGene',
                'ExonicFunc.refGene_x': 'ExonicFunc.refGene',
                '1GF_R_GF3KF6.GT':'1GF_R_final.GT',
                '1GF_A_GF3KF6.GT':'1GF_A_final.GT',
                '1GF_R_GF3KF6.AD':'1GF_R_final.AD',
                '1GF_A_GF3KF6.AD':'1GF_A_final.AD',
                '2GF_R_GF3KF6.GT':'2GF_R_final.GT',
                '2GF_A_GF3KF6.GT':'2GF_A_final.GT',
                '2GF_R_GF3KF6.AD':'2GF_R_final.AD',
                '2GF_A_GF3KF6.AD':'2GF_A_final.AD',
                '3GF_R_GF3KF6.GT':'3GF_R_final.GT',
                '3GF_A_GF3KF6.GT':'3GF_A_final.GT',
                '3GF_R_GF3KF6.AD':'3GF_R_final.AD',
                '3GF_A_GF3KF6.AD':'3GF_A_final.AD',
                '4GF_R_GF3KF6.GT':'4GF_R_final.GT',
                '4GF_A_GF3KF6.GT':'4GF_A_final.GT',
                '4GF_R_GF3KF6.AD':'4GF_R_final.AD',
                '4GF_A_GF3KF6.AD':'4GF_A_final.AD',
                '5GF_R_GF3KF6.GT':'5GF_R_final.GT',
                '5GF_A_GF3KF6.GT':'5GF_A_final.GT',
                '5GF_R_GF3KF6.AD':'5GF_R_final.AD',
                '5GF_A_GF3KF6.AD':'5GF_A_final.AD',
                '6GF_R_GF3KF6.GT':'6GF_R_final.GT',
                '6GF_A_GF3KF6.GT':'6GF_A_final.GT',
                '6GF_R_GF3KF6.AD':'6GF_R_final.AD',
                '6GF_A_GF3KF6.AD':'6GF_A_final.AD',
                '1GS_R_GF3KF6.GT':'1GS_R_final.GT',
                '1GS_A_GF3KF6.GT':'1GS_A_final.GT',
                '1GS_R_GF3KF6.AD':'1GS_R_final.AD',
                '1GS_A_GF3KF6.AD':'1GS_A_final.AD',
                '2GS_R_GF3KF6.GT':'2GS_R_final.GT',
                '2GS_A_GF3KF6.GT':'2GS_A_final.GT',
                '2GS_R_GF3KF6.AD':'2GS_R_final.AD',
                '2GS_A_GF3KF6.AD':'2GS_A_final.AD',
                '3GS_R_GF3KF6.GT':'3GS_R_final.GT',
                '3GS_A_GF3KF6.GT':'3GS_A_final.GT',
                '3GS_R_GF3KF6.AD':'3GS_R_final.AD',
                '3GS_A_GF3KF6.AD':'3GS_A_final.AD',
                '4GS_R_GF3KF6.GT':'4GS_R_final.GT',
                '4GS_A_GF3KF6.GT':'4GS_A_final.GT',
                '4GS_R_GF3KF6.AD':'4GS_R_final.AD',
                '4GS_A_GF3KF6.AD':'4GS_A_final.AD',
                '5GS_R_GF3KF6.GT':'5GS_R_final.GT',
                '5GS_A_GF3KF6.GT':'5GS_A_final.GT',
                '5GS_R_GF3KF6.AD':'5GS_R_final.AD',
                '5GS_A_GF3KF6.AD':'5GS_A_final.AD',
                '6GS_R_GF3KF6.GT':'6GS_R_final.GT',
                '6GS_A_GF3KF6.GT':'6GS_A_final.GT',
                '6GS_R_GF3KF6.AD':'6GS_R_final.AD',
                '6GS_A_GF3KF6.AD':'6GS_A_final.AD',
                '1KF_R_GF3KF6.GT':'1KF_R_final.GT',
                '1KF_A_GF3KF6.GT':'1KF_A_final.GT',
                '1KF_R_GF3KF6.AD':'1KF_R_final.AD',
                '1KF_A_GF3KF6.AD':'1KF_A_final.AD',
                '2KF_R_GF3KF6.GT':'2KF_R_final.GT',
                '2KF_A_GF3KF6.GT':'2KF_A_final.GT',
                '2KF_R_GF3KF6.AD':'2KF_R_final.AD',
                '2KF_A_GF3KF6.AD':'2KF_A_final.AD',
                '3KF_R_GF3KF6.GT':'3KF_R_final.GT',
                '3KF_A_GF3KF6.GT':'3KF_A_final.GT',
                '3KF_R_GF3KF6.AD':'3KF_R_final.AD',
                '3KF_A_GF3KF6.AD':'3KF_A_final.AD',
                '4KF_R_GF3KF6.GT':'4KF_R_final.GT',
                '4KF_A_GF3KF6.GT':'4KF_A_final.GT',
                '4KF_R_GF3KF6.AD':'4KF_R_final.AD',
                '4KF_A_GF3KF6.AD':'4KF_A_final.AD',
                '5KF_R_GF3KF6.GT':'5KF_R_final.GT',
                '5KF_A_GF3KF6.GT':'5KF_A_final.GT',
                '5KF_R_GF3KF6.AD':'5KF_R_final.AD',
                '5KF_A_GF3KF6.AD':'5KF_A_final.AD',
                '6KF_R_GF3KF6.GT':'6KF_R_final.GT',
                '6KF_A_GF3KF6.GT':'6KF_A_final.GT',
                '6KF_R_GF3KF6.AD':'6KF_R_final.AD',
                '6KF_A_GF3KF6.AD':'6KF_A_final.AD',
                '2KS_R_GF3KF6.GT':'2KS_R_final.GT',
                '2KS_A_GF3KF6.GT':'2KS_A_final.GT',
                '2KS_R_GF3KF6.AD':'2KS_R_final.AD',
                '2KS_A_GF3KF6.AD':'2KS_A_final.AD',
                '3KS_R_GF3KF6.GT':'3KS_R_final.GT',
                '3KS_A_GF3KF6.GT':'3KS_A_final.GT',
                '3KS_R_GF3KF6.AD':'3KS_R_final.AD',
                '3KS_A_GF3KF6.AD':'3KS_A_final.AD',
                '4KS_R_GF3KF6.GT':'4KS_R_final.GT',
                '4KS_A_GF3KF6.GT':'4KS_A_final.GT',
                '4KS_R_GF3KF6.AD':'4KS_R_final.AD',
                '4KS_A_GF3KF6.AD':'4KS_A_final.AD',
                '5KS_R_GF3KF6.GT':'5KS_R_final.GT',
                '5KS_A_GF3KF6.GT':'5KS_A_final.GT',
                '5KS_R_GF3KF6.AD':'5KS_R_final.AD',
                '5KS_A_GF3KF6.AD':'5KS_A_final.AD',
                '6KS_R_GF3KF6.GT':'6KS_R_final.GT',
                '6KS_A_GF3KF6.GT':'6KS_A_final.GT',
                '6KS_R_GF3KF6.AD':'6KS_R_final.AD',
                '6KS_A_GF3KF6.AD':'6KS_A_final.AD'},inplace=True)


# Collect also a SNP dictionary for the translation of SNPs based in the SNP ID with vlookup or similar

# In[6]:


SNP_dicc = df[['SNP_ID','CHROM','POS','Gene.refGene','Func.refGene','ExonicFunc.refGene',
               '1GF_R_final.GT',
               '1GF_A_final.GT',
               '2GF_R_final.GT',
               '2GF_A_final.GT',
               '3GF_R_final.GT',
               '3GF_A_final.GT',
               '4GF_R_final.GT',
               '4GF_A_final.GT',
               '5GF_R_final.GT',
               '5GF_A_final.GT',
               '6GF_R_final.GT',
               '6GF_A_final.GT',
               '1KF_R_final.GT',
               '1KF_A_final.GT',
               '2KF_R_final.GT',
               '2KF_A_final.GT',
               '3KF_R_final.GT',
               '3KF_A_final.GT',
               '4KF_R_final.GT',
               '4KF_A_final.GT',
               '5KF_R_final.GT',
               '5KF_A_final.GT',
               '6KF_R_final.GT',
               '6KF_A_final.GT',
               'AF_1GF',
               'AF_2GF',
               'AF_3GF',
               'AF_4GF',
               'AF_5GF',
               'AF_6GF',
               'AF_1GS',
               'AF_2GS',
               'AF_3GS',
               'AF_4GS',
               'AF_5GS',
               'AF_6GS',
               'AF_1KF',
               'AF_2KF',
               'AF_3KF',
               'AF_4KF',
               'AF_5KF',
               'AF_6KF',
               'AF_2KS',
               'AF_3KS',
               'AF_4KS',
               'AF_5KS',
               'AF_6KS',
               'AF',
               'AF_GF',
               'AF_GS',
               'AF_KF',
               'AF_KS']]
SNP_dicc.to_csv('SNP_dictionary_Weismann_2maps_GF3KF6.csv')


# In[7]:


df['AF_GF'].head()


# Make the total reads, including alternative and reference allele, per group

# In[8]:


df['GFT'] = df['1GF_R_final.AD']+df['1GF_A_final.AD']+df['2GF_R_final.AD']+df['2GF_A_final.AD']+df['3GF_R_final.AD']+df['3GF_A_final.AD']+df['4GF_R_final.AD']+df['4GF_A_final.AD']+df['5GF_R_final.AD']+df['5GF_A_final.AD']+df['6GF_R_final.AD']+df['6GF_A_final.AD']
df['GFT'] = df['GFT'].astype(int)


# In[9]:


df['GST'] = df['1GS_R_final.AD']+df['1GS_A_final.AD']+df['2GS_R_final.AD']+df['2GS_A_final.AD']+df['3GS_R_final.AD']+df['3GS_A_final.AD']+df['4GS_R_final.AD']+df['4GS_A_final.AD']+df['5GS_R_final.AD']+df['5GS_A_final.AD']+df['6GS_R_final.AD']+df['6GS_A_final.AD']
df['GST'] = df['GST'].astype(int)


# In[10]:


df['KFT'] = df['1KF_R_final.AD']+df['1KF_A_final.AD']+df['2KF_R_final.AD']+df['2KF_A_final.AD']+df['3KF_R_final.AD']+df['3KF_A_final.AD']+df['4KF_R_final.AD']+df['4KF_A_final.AD']+df['5KF_R_final.AD']+df['5KF_A_final.AD']+df['6KF_R_final.AD']+df['6KF_A_final.AD']
df['KFT'] = df['KFT'].astype(int)


# In[11]:


df['KST'] = df['2KS_R_final.AD']+df['2KS_A_final.AD']+df['3KS_R_final.AD']+df['3KS_A_final.AD']+df['4KS_R_final.AD']+df['4KS_A_final.AD']+df['5KS_R_final.AD']+df['5KS_A_final.AD']+df['6KS_R_final.AD']+df['6KS_A_final.AD']
df['KST'] = df['KST'].astype(int)


# Clean the rows containing alleles that have allele counts = 0 in O niloticus to prevent the collapse of the test. From now on, the matrix will be called df_binom.

# In[12]:


index_binom = df[((df['GFT'] == 0) & (df['GST'] == 0)&(df['KFT'] == 0) & (df['KST'] == 0))].index
index_binom


# In[13]:


df_binom = df.drop(index = index_binom)


# In[14]:


df_binom.shape


# In[15]:


for col in df_binom.columns: 
    print(col) 


# ## Heatmap

# In[16]:


import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mat


# Prepare the matrix of allele frequencies per group and Chromosome position. 

# Fist, divide the dataframe by genomic regions, spearating the linkage groups from the scaffolds

# In[17]:


df_binom.tail()


# In[18]:


index_scf = df_binom[df_binom['CHROM'].str.match(r'^NC_')].index
df_scf = df_binom.drop(index_scf)
df_scf.head()


# In[19]:


index_chrom = df_binom[df_binom['CHROM'].str.match(r'^NW_')].index
df_chrom = df_binom.drop(index_chrom)
df_chrom.head()


# Now make a general comparison between treatments. Most of the colors in these graphs are or 0 or 1, since they correspond to the allelic imbalance with significance.

# In[20]:


htmap = df_binom[['AF_GF','AF_GS','AF_KF','AF_KS']]
#htmap = df_binom[['AF_GF','AF_KF']]


# In[21]:


ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[22]:


htmap = df_chrom[['AF_GF','AF_GS','AF_KF','AF_KS']]
#htmap = df_chrom[['AF_GF','AF_KF']]


# In[23]:


ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[24]:


htmap = df_scf[['AF_GF','AF_GS','AF_KF','AF_KS']]
#htmap = df_scf[['AF_GF','AF_KF']]


# In[25]:


ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# Now make a graph for each treatment group

# In[26]:


htmap = df_chrom[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_GF')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[27]:


htmap = df_chrom[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_GS')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[28]:


htmap = df_chrom[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_KF')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[29]:


htmap = df_chrom[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_KS')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[30]:


htmap = df_scf[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_GF')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[31]:


htmap = df_scf[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_GS')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[32]:


htmap = df_scf[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_KF')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# In[33]:


htmap = df_scf[['CHROM','POS','AF_GF','AF_GS','AF_KF','AF_KS']]
htmap = htmap.pivot('POS','CHROM','AF_KS')
ax = sns.heatmap(htmap,cmap= 'YlGnBu')


# # Statistical tests

# We implement:
#     Chi² test for groups based in allele frequencies 
#     Binomial test for allelic imbalance based in allele counts
#     Fisher exact test for allele specific expression based in allele counts

# A Bonferroni one-step correction will be applied to the p-values of the chi-square test and a Benjamini-Hochberg non-negative test will be applied to the Fisher and Binomial tests respectively.

# In[34]:


import scipy.stats as stats
from scipy.stats import chisquare
import statsmodels.api as sm
import statsmodels.stats as smt


# ## Chi-square test for conditions based on allele frequencies

# In this test we compare the allelic imbalance in the gills and in the liver for a single individual. The list of genes that present AI will be compared between the two individuals. The chosen individual in fresh water is the 6 and the chosen individual in sea water is the 28.

#  scipy.stats.chisquare(f_obs, f_exp=None, ddof=0, axis=0)[source]¶
# 
#     Calculate a one-way chi-square test.
# 
#     The chi-square test tests the null hypothesis: the categorical data has the expected frequencies.
# 
#     Parameters
# 
#         f_obsarray_like
#             Observed frequencies in each category.
#         
#         f_exparray_like, optional
#             Expected frequencies in each category. By default the categories are assumed to be equally likely.
#         
#         ddofint, optional
#             “Delta degrees of freedom”: adjustment to the degrees of freedom for the p-value. The p-value is computed using a chi-squared distribution with k - 1 - ddof degrees of freedom, where k is the number of observed frequencies. The default value of ddof is 0.
#         
#         axisint or None, optional
#             The axis of the broadcast result of f_obs and f_exp along which to apply the test. If axis is None, all values in f_obs are treated as a single data set. Default is 0.
# 
# 
#     Returns
#     
#         chisqfloat or ndarray
#             The chi-squared test statistic. The value is a float if axis is None or f_obs and f_exp are 1-D.
#         pfloat or ndarray
#             The p-value of the test. The value is a float if ddof and the return value chisq are scalars.
# 
# 

# ## Chi-square for treatment group in each tissue

# This test will compare the differences between gills and liver in the O. niloticus. For this, we need to clean the dataset of the rows containing allele frequences equal to 0 for O. niloticus.
# The formula needs the next variables: sampling size, expected probability (reference allele frequency from control group) and observed probability (reference allele frequency from the experimental group)
# First test: gills and liver in the fresh water. Gills will be the control group. The second test will use the data of the sea water, also for gills and liver.
# Let's sum the values of the reference allele in liver and gills for the three individuals fresh and sea water and add it in a new column.

# We can use the data clean previously by droping the values where AF_NF is equal to 0. The dataset is called df_binom.

# Now let's perform the test.

# In[35]:


print(df_binom[['AF_GF','AF_GS','AF_KF','AF_KS']].head(10))


# The control group is the fresh water allele frequency, and the test will be in both tissues separately.

# In[36]:


chi_GS = []
p_chi_GS = []
stat_chi_GS = []
fdr_pval_GS = []
fdr_statistic_GS = []
for index, row in df_binom.iterrows() : 
    chi_GS.append(chisquare([row['AF_1GS'],row['AF_2GS'],row['AF_3GS'],row['AF_4GS'],row['AF_5GS'],row['AF_6GS']],
                              f_exp=[row['AF_1GF'],row['AF_2GF'],row['AF_3GF'],row['AF_4GF'],row['AF_5GF'],row['AF_6GF']]))
    p_chi_GS.append(chi_GS[index].pvalue)
    stat_chi_GS.append(chi_GS[index].statistic)
    #add the multitest also
    multi = smt.multitest.multipletests(chi_GS[index].pvalue, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    fdr_pval_GS.append(multi[1][0])
    fdr_statistic_GS.append(multi[0][0])


# In[37]:


df_binom['Chi_Gsalinity_pvalue'] = p_chi_GS


# In[38]:


df_binom['Chi_Gsalinity_p-fdr'] = fdr_pval_GS


# In[39]:


chi_KS = []
p_chi_KS = []
stat_chi_KS = []
fdr_pval_KS = []
fdr_statistic_KS = []

for index, row in df_binom.iterrows() : 
    chi_KS.append(chisquare([row['AF_2KS'],row['AF_3KS'],row['AF_4KS'],row['AF_5KS'],row['AF_6KS']],
                              f_exp=[row['AF_2KF'],row['AF_3KF'],row['AF_4KF'],row['AF_5KF'],row['AF_6KF']]))
    p_chi_KS.append(chi_KS[index].pvalue)
    stat_chi_KS.append(chi_KS[index].statistic)
    #add the multitest also
    multi = smt.multitest.multipletests(chi_KS[index].pvalue, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    fdr_pval_KS.append(multi[1][0])
    fdr_statistic_KS.append(multi[0][0])


# In[40]:


df_binom['Chi_Ksalinity_pvalue'] = p_chi_KS


# In[41]:


df_binom['Chi_Ksalinity_p-fdr'] = fdr_pval_KS


# Another test will be the different expression of the tissues: gills against kidney. The control group will be the gills and this test wil be performed only in the fresh water groups

# In[42]:


chi_KF = []
p_chi_KF = []
stat_chi_KF = []
fdr_pval_KF = []
fdr_statistic_KF = []

for index, row in df_binom.iterrows() : 
    chi_KF.append(chisquare([row['AF_1KF'],row['AF_2KF'],row['AF_3KF'],row['AF_4KF'],row['AF_5KF'],row['AF_6KF']],
                              f_exp=[row['AF_1GF'],row['AF_2GF'],row['AF_3GF'],row['AF_4GF'],row['AF_5GF'],row['AF_6GF']]))
    p_chi_KF.append(chi_KF[index].pvalue)
    stat_chi_KF.append(chi_KF[index].statistic)
    #add the multitest also
    multi = smt.multitest.multipletests(chi_KF[index].pvalue, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    fdr_pval_KF.append(multi[1][0])
    fdr_statistic_KF.append(multi[0][0])


# In[43]:


df_binom['Chi_tissue_pvalue'] = p_chi_KF


# In[44]:


df_binom['Chi_tissue_p-fdr'] = fdr_pval_KF


# In[45]:


df_binom.columns


# # Fisher exact test for ASE in tissues

# This test is based on the read counts of each allele, including alternative and reference alleles. It is valid to test one tissue to the other for the same individual.

# In[46]:


p_Fisher_KF_1 = []
stat_Fisher_KF_1 = []
fdr_pval_KF1 = []
fdr_statistic_KF1 = []
for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['1KF_R_final.AD'],
                                              row['1KF_A_final.AD']], 
                                             [row['1GF_R_final.AD'],
                                              row['1GF_A_final.AD']]])
    stat_Fisher_KF_1.append(oddsratio)
    p_Fisher_KF_1.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KF_1[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KF1.append(multi[1][0])
    fdr_statistic_KF1.append(multi[0][0])
    


# In[47]:


df_binom['F1_Fisher_tissue_pvalue'] = p_Fisher_KF_1


# In[48]:


df_binom['F1_Fisher_tissue_p-fdr'] = fdr_pval_KF1


# In[49]:


p_Fisher_KF_2 = []
stat_Fisher_KF_2 = []
fdr_pval_KF2 = []
fdr_statistic_KF2 = []
for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['2KF_R_final.AD'],
                                              row['2KF_A_final.AD']], 
                                             [row['2GF_R_final.AD'],
                                              row['2GF_A_final.AD']]])
    stat_Fisher_KF_2.append(oddsratio)
    p_Fisher_KF_2.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KF_2[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KF2.append(multi[1][0])
    fdr_statistic_KF2.append(multi[0][0])


# In[50]:


df_binom['F2_Fisher_tissue_pvalue'] = p_Fisher_KF_2


# In[51]:


df_binom['F2_Fisher_tissue_p-fdr'] = fdr_pval_KF2


# In[52]:


p_Fisher_KF_3 = []
stat_Fisher_KF_3 = []
fdr_pval_KF3 = []
fdr_statistic_KF3 = []
for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['3KF_R_final.AD'],
                                              row['3KF_A_final.AD']], 
                                             [row['3GF_R_final.AD'],
                                              row['3GF_A_final.AD']]])
    stat_Fisher_KF_3.append(oddsratio)
    p_Fisher_KF_3.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KF_3[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KF3.append(multi[1][0])
    fdr_statistic_KF3.append(multi[0][0])


# In[53]:


df_binom['F3_Fisher_tissue_pvalue'] = p_Fisher_KF_3


# In[54]:


df_binom['F3_Fisher_tissue_p-fdr'] = fdr_pval_KF3


# In[55]:


p_Fisher_KF_4 = []
stat_Fisher_KF_4 = []
fdr_pval_KF4 = []
fdr_statistic_KF4 = []

for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['4KF_R_final.AD'],
                                              row['4KF_A_final.AD']], 
                                             [row['4GF_R_final.AD'],
                                              row['4GF_A_final.AD']]])
    stat_Fisher_KF_4.append(oddsratio)
    p_Fisher_KF_4.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KF_4[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KF4.append(multi[1][0])
    fdr_statistic_KF4.append(multi[0][0])


# In[56]:


df_binom['F4_Fisher_tissue_pvalue'] = p_Fisher_KF_4


# In[57]:


df_binom['F4_Fisher_tissue_p-fdr'] = fdr_pval_KF4


# In[58]:


p_Fisher_KF_5 = []
stat_Fisher_KF_5 = []
fdr_pval_KF5 = []
fdr_statistic_KF5 = []

for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['5KF_R_final.AD'],
                                              row['5KF_A_final.AD']], 
                                             [row['5GF_R_final.AD'],
                                              row['5GF_A_final.AD']]])
    stat_Fisher_KF_5.append(oddsratio)
    p_Fisher_KF_5.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KF_5[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KF5.append(multi[1][0])
    fdr_statistic_KF5.append(multi[0][0])


# In[59]:


df_binom['F5_Fisher_tissue_pvalue'] = p_Fisher_KF_5


# In[60]:


df_binom['F5_Fisher_tissue_p-fdr'] = fdr_pval_KF5


# In[61]:


p_Fisher_KF_6 = []
stat_Fisher_KF_6 = []
fdr_pval_KF6 = []
fdr_statistic_KF6 = []
for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['6KF_R_final.AD'],
                                              row['6KF_A_final.AD']], 
                                             [row['6GF_R_final.AD'],
                                              row['6GF_A_final.AD']]])
    stat_Fisher_KF_6.append(oddsratio)
    p_Fisher_KF_6.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KF_6[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KF6.append(multi[1][0])
    fdr_statistic_KF6.append(multi[0][0])


# In[62]:


df_binom['F6_Fisher_tissue_pvalue'] = p_Fisher_KF_6


# In[63]:


df_binom['F6_Fisher_tissue_p-fdr'] = fdr_pval_KF6


# In[64]:


p_Fisher_KS_2 = []
stat_Fisher_KS_2 = []
fdr_pval_KS2 = []
fdr_statistic_KS2 = []
for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['2KS_R_final.AD'],
                                              row['2KS_A_final.AD']], 
                                             [row['2GS_R_final.AD'],
                                              row['2GS_A_final.AD']]])
    stat_Fisher_KS_2.append(oddsratio)
    p_Fisher_KS_2.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KS_2[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KS2.append(multi[1][0])
    fdr_statistic_KS2.append(multi[0][0])


# In[65]:


df_binom['S2_Fisher_tissue_pvalue'] = p_Fisher_KS_2


# In[66]:


df_binom['S2_Fisher_tissue_p-fdr'] = fdr_pval_KS2


# In[67]:


p_Fisher_KS_3 = []
stat_Fisher_KS_3 = []
fdr_pval_KS3 = []
fdr_statistic_KS3 = []

for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['3KS_R_final.AD'],
                                              row['3KS_A_final.AD']], 
                                             [row['3GS_R_final.AD'],
                                              row['3GS_A_final.AD']]])
    stat_Fisher_KS_3.append(oddsratio)
    p_Fisher_KS_3.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KS_3[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KS3.append(multi[1][0])
    fdr_statistic_KS3.append(multi[0][0])


# In[68]:


df_binom['S3_Fisher_tissue_pvalue'] = p_Fisher_KS_3


# In[69]:


df_binom['S3_Fisher_tissue_p-fdr'] = fdr_pval_KS3


# In[70]:


p_Fisher_KS_4 = []
stat_Fisher_KS_4 = []
fdr_pval_KS4 = []
fdr_statistic_KS4 = []

for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['4KS_R_final.AD'],
                                              row['4KS_A_final.AD']], 
                                             [row['4GS_R_final.AD'],
                                              row['4GS_A_final.AD']]])
    stat_Fisher_KS_4.append(oddsratio)
    p_Fisher_KS_4.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KS_4[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KS4.append(multi[1][0])
    fdr_statistic_KS4.append(multi[0][0])


# In[71]:


df_binom['S4_Fisher_tissue_pvalue'] = p_Fisher_KS_4


# In[72]:


df_binom['S4_Fisher_tissue_p-fdr'] = fdr_pval_KS4


# In[73]:


p_Fisher_KS_5 = []
stat_Fisher_KS_5 = []
fdr_pval_KS5 = []
fdr_statistic_KS5 = []
for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['5KS_R_final.AD'],
                                              row['5KS_A_final.AD']], 
                                             [row['5GS_R_final.AD'],
                                              row['5GS_A_final.AD']]])
    stat_Fisher_KS_5.append(oddsratio)
    p_Fisher_KS_5.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KS_5[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KS5.append(multi[1][0])
    fdr_statistic_KS5.append(multi[0][0])


# In[74]:


df_binom['S5_Fisher_tissue_pvalue'] = p_Fisher_KS_5


# In[75]:


df_binom['S5_Fisher_tissue_p-fdr'] = fdr_pval_KS5


# In[76]:


p_Fisher_KS_6 = []
stat_Fisher_KS_6 = []
fdr_pval_KS6 = []
fdr_statistic_KS6 = []

for index, row in df_binom.iterrows() :
    
    oddsratio, pvalue = stats.fisher_exact([[row['6KS_R_final.AD'],
                                              row['6KS_A_final.AD']], 
                                             [row['6GS_R_final.AD'],
                                              row['6GS_A_final.AD']]])
    stat_Fisher_KS_6.append(oddsratio)
    p_Fisher_KS_6.append(pvalue)
    #add the multitest also
    multi = smt.multitest.multipletests(p_Fisher_KS_6[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_KS6.append(multi[1][0])
    fdr_statistic_KS6.append(multi[0][0])


# In[77]:


df_binom['S6_Fisher_tissue_pvalue'] = p_Fisher_KS_6


# In[78]:


df_binom['S6_Fisher_tissue_p-fdr'] = fdr_pval_KS6


# # Binomial test for AI in samples

# The binomial test will determine the allelic imbalance in each sample by comparing the counts of the reference allele with the total counts of the SNP.

# ## Analysis for group data

# Perform the binomial tests in the read counts of the reference allele (it can be done by experimental group or by individual)(x) over the total number of reads in this group (as before, it can be done by experimental group or by individual) (n). For this analysis, the expected frequency is 0.5 and the observed allele frequency is determined for the reference allele in each individual or for the control group (p). This test will provide the p-values for the significant differences of the groups. 
# The formula is extracted from wikipedia and adapted to each test. 
# In python the formula is scipy.stats.binom_test (x,n,p, alternative "two-sided").

# In this test we compare the allelic imbalance in the gills and in the liver for each individual of the Oreochromis niloticus. The list of genes that present AI will be compared between the tissues and the treatments in an analysis group manner: outer and inner genes within each group.

# First prepare the matrix for the binomial analysis and the subset of the GeneID

# In[79]:


df_binom.rename(columns={"Gene.refGene": "GeneID"}, inplace = True)


# ## Binomial test for all the samples

# Calculate the total reads for reference and alternative allele in each sample, each tissue of each individual

# In[80]:


df_binom['1GFT'] = df_binom['1GF_R_final.AD']+df_binom['1GF_A_final.AD']
df_binom['1GFT'] = df_binom['1GFT'].astype(int)

df_binom['2GFT'] = df_binom['2GF_R_final.AD']+df_binom['2GF_A_final.AD']
df_binom['2GFT'] = df_binom['2GFT'].astype(int)

df_binom['3GFT'] = df_binom['3GF_R_final.AD']+df_binom['3GF_A_final.AD']
df_binom['3GFT'] = df_binom['3GFT'].astype(int)

df_binom['4GFT'] = df_binom['4GF_R_final.AD']+df_binom['4GF_A_final.AD']
df_binom['4GFT'] = df_binom['4GFT'].astype(int)

df_binom['5GFT'] = df_binom['5GF_R_final.AD']+df_binom['5GF_A_final.AD']
df_binom['5GFT'] = df_binom['5GFT'].astype(int)

df_binom['6GFT'] = df_binom['6GF_R_final.AD']+df_binom['6GF_A_final.AD']
df_binom['6GFT'] = df_binom['6GFT'].astype(int)


# In[81]:


df_binom['1GST'] = df_binom['1GS_R_final.AD']+df_binom['1GS_A_final.AD']
df_binom['1GST'] = df_binom['1GST'].astype(int)

df_binom['2GST'] = df_binom['2GS_R_final.AD']+df_binom['2GS_A_final.AD']
df_binom['2GST'] = df_binom['2GST'].astype(int)

df_binom['3GST'] = df_binom['3GS_R_final.AD']+df_binom['3GS_A_final.AD']
df_binom['3GST'] = df_binom['3GST'].astype(int)

df_binom['4GST'] = df_binom['4GS_R_final.AD']+df_binom['4GS_A_final.AD']
df_binom['4GST'] = df_binom['4GST'].astype(int)

df_binom['5GST'] = df_binom['5GS_R_final.AD']+df_binom['5GS_A_final.AD']
df_binom['5GST'] = df_binom['5GST'].astype(int)

df_binom['6GST'] = df_binom['6GS_R_final.AD']+df_binom['6GS_A_final.AD']
df_binom['6GST'] = df_binom['6GST'].astype(int)


# In[82]:


df_binom['1KFT'] = df_binom['1KF_R_final.AD']+df_binom['1KF_A_final.AD']
df_binom['1KFT'] = df_binom['1KFT'].astype(int)

df_binom['2KFT'] = df_binom['2KF_R_final.AD']+df_binom['2KF_A_final.AD']
df_binom['2KFT'] = df_binom['2KFT'].astype(int)

df_binom['3KFT'] = df_binom['3KF_R_final.AD']+df_binom['3KF_A_final.AD']
df_binom['3KFT'] = df_binom['3KFT'].astype(int)

df_binom['4KFT'] = df_binom['4KF_R_final.AD']+df_binom['4KF_A_final.AD']
df_binom['4KFT'] = df_binom['4KFT'].astype(int)

df_binom['5KFT'] = df_binom['5KF_R_final.AD']+df_binom['5KF_A_final.AD']
df_binom['5KFT'] = df_binom['5KFT'].astype(int)

df_binom['6KFT'] = df_binom['6KF_R_final.AD']+df_binom['6KF_A_final.AD']
df_binom['6KFT'] = df_binom['6KFT'].astype(int)


# In[83]:


df_binom['2KST'] = df_binom['2KS_R_final.AD']+df_binom['2KS_A_final.AD']
df_binom['2KST'] = df_binom['2KST'].astype(int)

df_binom['3KST'] = df_binom['3KS_R_final.AD']+df_binom['3KS_A_final.AD']
df_binom['3KST'] = df_binom['3KST'].astype(int)

df_binom['4KST'] = df_binom['4KS_R_final.AD']+df_binom['4KS_A_final.AD']
df_binom['4KST'] = df_binom['4KST'].astype(int)

df_binom['5KST'] = df_binom['5KS_R_final.AD']+df_binom['5KS_A_final.AD']
df_binom['5KST'] = df_binom['5KST'].astype(int)

df_binom['6KST'] = df_binom['6KS_R_final.AD']+df_binom['6KS_A_final.AD']
df_binom['6KST'] = df_binom['6KST'].astype(int)


# Make a for loop to calculate the p-value of the binomial test per each row and create a new column for it:

# In[84]:


p_binom_1GF=[]
p_binom_2GF=[]
p_binom_3GF=[]
p_binom_4GF=[]
p_binom_5GF=[]
p_binom_6GF=[]


# In[85]:


p_binom_1GS=[]
p_binom_2GS=[]
p_binom_3GS=[]
p_binom_4GS=[]
p_binom_5GS=[]
p_binom_6GS=[]


# In[86]:


p_binom_1KF=[]
p_binom_2KF=[]
p_binom_3KF=[]
p_binom_4KF=[]
p_binom_5KF=[]
p_binom_6KF=[]


# In[87]:


p_binom_2KS=[]
p_binom_3KS=[]
p_binom_4KS=[]
p_binom_5KS=[]
p_binom_6KS=[]


# In[88]:


fdr_pval_bi1GF = []
fdr_pval_bi2GF = []
fdr_pval_bi3GF = []
fdr_pval_bi4GF = []
fdr_pval_bi5GF = []
fdr_pval_bi6GF = []
fdr_pval_bi1GS = []
fdr_pval_bi2GS = []
fdr_pval_bi3GS = []
fdr_pval_bi4GS = []
fdr_pval_bi5GS = []
fdr_pval_bi6GS = []
fdr_pval_bi1KF = []
fdr_pval_bi2KF = []
fdr_pval_bi3KF = []
fdr_pval_bi4KF = []
fdr_pval_bi5KF = []
fdr_pval_bi6KF = []
fdr_pval_bi2KS = []
fdr_pval_bi3KS = []
fdr_pval_bi4KS = []
fdr_pval_bi5KS = []
fdr_pval_bi6KS = []


# Perform the test

# In[89]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['1GF_R_final.AD'], row['1GFT'], 0.5, alternative = 'two-sided'))
    p_binom_1GF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_1GF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi1GF.append(multi[1][0])


# In[90]:


df_binom['binom_1GF'] = p_binom_1GF    
df_binom['binom_1GF_p-fdr'] = fdr_pval_bi1GF


# In[91]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['2GF_R_final.AD'], row['2GFT'], 0.5, alternative = 'two-sided'))
    p_binom_2GF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_2GF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi2GF.append(multi[1][0])


# In[92]:


df_binom['binom_2GF'] = p_binom_2GF     
df_binom['binom_2GF_p-fdr'] = fdr_pval_bi2GF


# In[93]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['3GF_R_final.AD'], row['3GFT'], 0.5, alternative = 'two-sided'))
    p_binom_3GF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_3GF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi3GF.append(multi[1][0])


# In[94]:


df_binom['binom_3GF'] = p_binom_3GF     
df_binom['binom_3GF_p-fdr'] = fdr_pval_bi3GF


# In[95]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['4GF_R_final.AD'], row['4GFT'], 0.5, alternative = 'two-sided'))
    p_binom_4GF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_4GF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi4GF.append(multi[1][0])


# In[96]:


df_binom['binom_4GF'] = p_binom_4GF    
df_binom['binom_4GF_p-fdr'] = fdr_pval_bi4GF


# In[97]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['5GF_R_final.AD'], row['5GFT'], 0.5, alternative = 'two-sided'))
    p_binom_5GF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_5GF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi5GF.append(multi[1][0])


# In[98]:


df_binom['binom_5GF'] = p_binom_5GF     
df_binom['binom_5GF_p-fdr'] = fdr_pval_bi5GF


# In[99]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['6GF_R_final.AD'], row['6GFT'], 0.5, alternative = 'two-sided'))
    p_binom_6GF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_6GF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi6GF.append(multi[1][0])


# In[100]:


df_binom['binom_6GF'] = p_binom_6GF     
df_binom['binom_6GF_p-fdr'] = fdr_pval_bi6GF


# In[101]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['1GS_R_final.AD'], row['1GST'], 0.5, alternative = 'two-sided'))
    p_binom_1GS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_1GS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi1GS.append(multi[1][0])


# In[102]:


df_binom['binom_1GS'] = p_binom_1GS     
df_binom['binom_1GS_p-fdr'] = fdr_pval_bi1GS


# In[103]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['2GS_R_final.AD'], row['2GST'], 0.5, alternative = 'two-sided'))
    p_binom_2GS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_2GS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi2GS.append(multi[1][0])


# In[104]:


df_binom['binom_2GS'] = p_binom_2GS     
df_binom['binom_2GS_p-fdr'] = fdr_pval_bi2GS


# In[105]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['3GS_R_final.AD'], row['3GST'], 0.5, alternative = 'two-sided'))
    p_binom_3GS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_3GS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi3GS.append(multi[1][0])


# In[106]:


df_binom['binom_3GS'] = p_binom_3GS     
df_binom['binom_3GS_p-fdr'] = fdr_pval_bi3GS


# In[107]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['4GS_R_final.AD'], row['4GST'], 0.5, alternative = 'two-sided'))
    p_binom_4GS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_4GS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi4GS.append(multi[1][0])


# In[108]:


df_binom['binom_4GS'] = p_binom_4GS     
df_binom['binom_4GS_p-fdr'] = fdr_pval_bi4GS


# In[109]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['5GS_R_final.AD'], row['5GST'], 0.5, alternative = 'two-sided'))
    p_binom_5GS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_5GS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi5GS.append(multi[1][0])


# In[110]:


df_binom['binom_5GS'] = p_binom_5GS     
df_binom['binom_5GS_p-fdr'] = fdr_pval_bi5GS


# In[111]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['6GS_R_final.AD'], row['6GST'], 0.5, alternative = 'two-sided'))
    p_binom_6GS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_6GS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi6GS.append(multi[1][0])


# In[112]:


df_binom['binom_6GS'] = p_binom_6GS     
df_binom['binom_6GS_p-fdr'] = fdr_pval_bi6GS


# In[113]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['1KF_R_final.AD'], row['1KFT'], 0.5, alternative = 'two-sided'))
    p_binom_1KF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_1KF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi1KF.append(multi[1][0])


# In[114]:


df_binom['binom_1KF'] = p_binom_1KF     
df_binom['binom_1KF_p-fdr'] = fdr_pval_bi1KF


# In[115]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['2KF_R_final.AD'], row['2KFT'], 0.5, alternative = 'two-sided'))
    p_binom_2KF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_2KF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi2KF.append(multi[1][0])


# In[116]:


df_binom['binom_2KF'] = p_binom_2KF     
df_binom['binom_2KF_p-fdr'] = fdr_pval_bi2KF


# In[117]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['3KF_R_final.AD'], row['3KFT'], 0.5, alternative = 'two-sided'))
    p_binom_3KF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_3KF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi3KF.append(multi[1][0])


# In[118]:


df_binom['binom_3KF'] = p_binom_3KF     
df_binom['binom_3KF_p-fdr'] = fdr_pval_bi3KF


# In[119]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['4KF_R_final.AD'], row['4KFT'], 0.5, alternative = 'two-sided'))
    p_binom_4KF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_4KF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi4KF.append(multi[1][0])


# In[120]:


df_binom['binom_4KF'] = p_binom_4KF     
df_binom['binom_4KF_p-fdr'] = fdr_pval_bi4KF


# In[121]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['5KF_R_final.AD'], row['5KFT'], 0.5, alternative = 'two-sided'))
    p_binom_5KF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_5KF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi5KF.append(multi[1][0])


# In[122]:


df_binom['binom_5KF'] = p_binom_5KF     
df_binom['binom_5KF_p-fdr'] = fdr_pval_bi5KF


# In[123]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['6KF_R_final.AD'], row['6KFT'], 0.5, alternative = 'two-sided'))
    p_binom_6KF.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_6KF[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi6KF.append(multi[1][0])


# In[124]:


df_binom['binom_6KF'] = p_binom_6KF     
df_binom['binom_6KF_p-fdr'] = fdr_pval_bi6KF


# In[125]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['2KS_R_final.AD'], row['2KST'], 0.5, alternative = 'two-sided'))
    p_binom_2KS.append(binom)
#add the multitest also
    multi = smt.multitest.multipletests(p_binom_2KS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi2KS.append(multi[1][0])


# In[126]:


df_binom['binom_2KS'] = p_binom_2KS     
df_binom['binom_2KS_p-fdr'] = fdr_pval_bi2KS


# In[127]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['3KS_R_final.AD'], row['3KST'], 0.5, alternative = 'two-sided'))
    p_binom_3KS.append(binom)
#add the multitest also
    multi = smt.multitest.multipletests(p_binom_3KS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi3KS.append(multi[1][0])


# In[128]:


df_binom['binom_3KS'] = p_binom_3KS     
df_binom['binom_3KS_p-fdr'] = fdr_pval_bi3KS


# In[129]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['4KS_R_final.AD'], row['4KST'], 0.5, alternative = 'two-sided'))
    p_binom_4KS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_4KS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi4KS.append(multi[1][0])


# In[130]:


df_binom['binom_4KS'] = p_binom_4KS     
df_binom['binom_4KS_p-fdr'] = fdr_pval_bi4KS


# In[131]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['5KS_R_final.AD'], row['5KST'], 0.5, alternative = 'two-sided'))
    p_binom_5KS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_5KS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi5KS.append(multi[1][0])


# In[132]:


df_binom['binom_5KS'] = p_binom_5KS     
df_binom['binom_5KS_p-fdr'] = fdr_pval_bi5KS


# In[133]:


for index, row in df_binom.iterrows() : 
    binom=(sc.stats.binom_test(row['6KS_R_final.AD'], row['6KST'], 0.5, alternative = 'two-sided'))
    p_binom_6KS.append(binom)
    #add the multitest also
    multi = smt.multitest.multipletests(p_binom_6KS[index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    fdr_pval_bi6KS.append(multi[1][0])


# In[134]:


df_binom['binom_6KS'] = p_binom_6KS     
df_binom['binom_6KS_p-fdr'] = fdr_pval_bi6KS


# And more data summaries for gene function and whitin the exonic function

# In[135]:


df_binom['Func.refGene'].value_counts()


# In[136]:


df_binom['ExonicFunc.refGene'].value_counts()


# The datafframe df_binom includes p-values for the binomial tests and all the allele frequencies from all the samples that is necessary to create all the next reports. The significant p-values have not been filtered here in order to be used in the analysis of the SNPs in AI in one organ and not AI for the other organ withing the same individual/experimental group. The section title in this script is called Differential imbalance depeding on tissue.

# In[137]:


df_binom.columns


# In[138]:


df_binom.drop(['Unnamed: 0', 
               'Unnamed: 0.1',
               'Unnamed: 0.1.1',
               'Unnamed: 0.1.1.1'], axis=1, inplace=True)


# In[139]:


df_binom.to_csv('Total_SNPs_pval_GF3KF6.csv')


# # Differential imbalance depeding on tissue

# In this analysis we will identify the SNPs that are balanced in one tissue and not in the other, for gills and liver. AI in this case is independent of ASE, since the gene will be expressed in both tissues. The objective of this analysis is to find SNPs in a situation of opposite AI between salinities in each organ or in both organs. The final table will be a list of SNPs that can be input in a Venn diagram.

# First prepare for the total SNPs before filtering for the AI significance in each sample.

# In[140]:


df_total=pd.read_csv('Total_SNPs_pval_GF3KF6.csv')


# In[141]:


for col in df_total.columns: 
    print(col) 


# Clean the dataframe to get only the p-values and the allele frequencies of the groups, as key data for tracking the SNPs such as SNP_ID, chromosome and position.

# In[142]:


df_total = df_total[['CHROM',
                    'POS',
                    'SNP_ID',
                    'GeneID',
                    'Func.refGene',
                    'ExonicFunc.refGene',
                    'AF',
                    'AF_GF',
                    'AF_GS',
                    'AF_KF',
                    'AF_KS',
                    'GFT',
                    'GST',
                    'KFT',
                    'KST',
                    'binom_1GF',
                    'binom_2GF',
                    'binom_3GF',
                    'binom_4GF',
                    'binom_5GF',
                    'binom_6GF',
                    'binom_1KF',
                    'binom_2KF',
                    'binom_3KF',
                    'binom_4KF',
                    'binom_5KF',
                    'binom_6KF',
                    'binom_1GS',
                    'binom_2GS',
                    'binom_3GS',
                    'binom_4GS',
                    'binom_5GS',
                    'binom_6GS',
                    'binom_2KS',
                    'binom_3KS',
                    'binom_4KS',
                    'binom_5KS',
                    'binom_6KS']]


# Make the test for the different AI between tissues in the same individual. To do so, compare which indexes (which rows) contain a significant p-value for the binomial test performed in the liver and at the same time, this index is not present in the significant p-values of the binomial test performed in the gills of the same individual.

# In[143]:


AI1GF = df_total[(df_total['binom_1KF'] < 0.05)&(df_total['binom_1GF'] >= 0.05)]
AI1KF = df_total[(df_total['binom_1GF'] < 0.05)&(df_total['binom_1KF'] >= 0.05)]

AI2GF = df_total[(df_total['binom_2KF'] < 0.05)&(df_total['binom_2GF'] >= 0.05)]
AI2KF = df_total[(df_total['binom_2GF'] < 0.05)&(df_total['binom_2KF'] >= 0.05)]

AI3GF = df_total[(df_total['binom_3KF'] < 0.05)&(df_total['binom_3GF'] >= 0.05)]
AI3KF = df_total[(df_total['binom_3GF'] < 0.05)&(df_total['binom_3KF'] >= 0.05)]

AI4GF = df_total[(df_total['binom_4KF'] < 0.05)&(df_total['binom_4GF'] >= 0.05)]
AI4KF = df_total[(df_total['binom_4GF'] < 0.05)&(df_total['binom_4KF'] >= 0.05)]

AI5GF = df_total[(df_total['binom_5KF'] < 0.05)&(df_total['binom_5GF'] >= 0.05)]
AI5KF = df_total[(df_total['binom_5GF'] < 0.05)&(df_total['binom_5KF'] >= 0.05)]

AI6GF = df_total[(df_total['binom_6KF'] < 0.05)&(df_total['binom_6GF'] >= 0.05)]
AI6KF = df_total[(df_total['binom_6GF'] < 0.05)&(df_total['binom_6KF'] >= 0.05)]


# In[144]:


#AI1GS = df_total[(df_total['binom_1KS'] < 0.05)&(df_total['binom_1GS'] >= 0.05)]
#AI1KS = df_total[(df_total['binom_1GS'] < 0.05)&(df_total['binom_1KS'] >= 0.05)]

AI2GS = df_total[(df_total['binom_2KS'] < 0.05)&(df_total['binom_2GS'] >= 0.05)]
AI2KS = df_total[(df_total['binom_2GS'] < 0.05)&(df_total['binom_2KS'] >= 0.05)]

AI3GS = df_total[(df_total['binom_3KS'] < 0.05)&(df_total['binom_3GS'] >= 0.05)]
AI3KS = df_total[(df_total['binom_3GS'] < 0.05)&(df_total['binom_3KS'] >= 0.05)]

AI4GS = df_total[(df_total['binom_4KS'] < 0.05)&(df_total['binom_4GS'] >= 0.05)]
AI4KS = df_total[(df_total['binom_4GS'] < 0.05)&(df_total['binom_4KS'] >= 0.05)]

AI5GS = df_total[(df_total['binom_5KS'] < 0.05)&(df_total['binom_5GS'] >= 0.05)]
AI5KS = df_total[(df_total['binom_5GS'] < 0.05)&(df_total['binom_5KS'] >= 0.05)]

AI6GS = df_total[(df_total['binom_6KS'] < 0.05)&(df_total['binom_6GS'] >= 0.05)]
AI6KS = df_total[(df_total['binom_6GS'] < 0.05)&(df_total['binom_6KS'] >= 0.05)]


# Collect the SNP_ID AI independently of ASE of each individual in arrays to make a table that can be printed and translated.

# In[145]:


AI1GF=np.array(AI1GF.SNP_ID)
AI1KF=np.array(AI1KF.SNP_ID)

AI2GF=np.array(AI2GF.SNP_ID)
AI2KF=np.array(AI2KF.SNP_ID)

AI3GF=np.array(AI3GF.SNP_ID)
AI3KF=np.array(AI3KF.SNP_ID)

AI4GF=np.array(AI4GF.SNP_ID)
AI4KF=np.array(AI4KF.SNP_ID)

AI5GF=np.array(AI5GF.SNP_ID)
AI5KF=np.array(AI5KF.SNP_ID)

AI6GF=np.array(AI6GF.SNP_ID)
AI6KF=np.array(AI6KF.SNP_ID)


# In[146]:


#AI1GS=np.array(AI1GS.SNP_ID)
#AI1KS=np.array(AI1KS.SNP_ID)

AI2GS=np.array(AI2GS.SNP_ID)
AI2KS=np.array(AI2KS.SNP_ID)

AI3GS=np.array(AI3GS.SNP_ID)
AI3KS=np.array(AI3KS.SNP_ID)

AI4GS=np.array(AI4GS.SNP_ID)
AI4KS=np.array(AI4KS.SNP_ID)

AI5GS=np.array(AI5GS.SNP_ID)
AI5KS=np.array(AI5KS.SNP_ID)

AI6GS=np.array(AI6GS.SNP_ID)
AI6KS=np.array(AI6KS.SNP_ID)


# In[147]:


AI_ID = pd.DataFrame(data = [AI1GF,
                             AI1KF,
                             AI2GF,
                             AI2KF,
                             AI3GF,
                             AI3KF,
                             AI4GF,
                             AI4KF,
                             AI5GF,
                             AI5KF,
                             AI6GF,
                             AI6KF,
                             AI2GS,
                             AI2KS,
                             AI3GS,
                             AI3KS,
                             AI4GS,
                             AI4KS,
                             AI5GS,
                             AI5KS,
                             AI6GS,
                             AI6KS]).T
AI_ID = AI_ID.rename(columns={0:'AI1GF',
                             1:'AI1KF',
                             2:'AI2GF',
                             3:'AI2KF',
                             4:'AI3GF',
                             5:'AI3KF',
                             6:'AI4GF',
                             7:'AI4KF',
                             8:'AI5GF',
                             9:'AI5KF',
                             10:'AI6GF',
                             11:'AI6KF',
                             12:'AI2GS',
                             13:'AI2KS',
                             14:'AI3GS',
                             15:'AI3KS',
                             16:'AI4GS',
                             17:'AI4KS',
                             18:'AI5GS',
                             19:'AI5KS',
                             20:'AI6GS',
                             21:'AI6KS'})
AI_ID.head()


# In[148]:


AI_ID.to_csv('AI_tissue_ID_FW_GF3KF6.csv')


# In[149]:


np.savetxt("AI1GF.csv", AI1GF, delimiter=",", header = 'SNP_ID')


# In[150]:


np.savetxt("AI1KF.csv", AI1KF, delimiter=",", header = 'SNP_ID')


# In[151]:


np.savetxt("AI2GF.csv", AI2GF, delimiter=",", header = 'SNP_ID')


# In[152]:


np.savetxt("AI2KF.csv", AI2KF, delimiter=",", header = 'SNP_ID')


# In[153]:


np.savetxt("AI3GF.csv", AI3GF, delimiter=",", header = 'SNP_ID')


# In[154]:


np.savetxt("AI3KF.csv", AI3KF, delimiter=",", header = 'SNP_ID')


# In[155]:


np.savetxt("AI4GF.csv", AI4GF, delimiter=",", header = 'SNP_ID')


# In[156]:


np.savetxt("AI4KF.csv", AI4KF, delimiter=",", header = 'SNP_ID')


# In[157]:


np.savetxt("AI5GF.csv", AI5GF, delimiter=",", header = 'SNP_ID')


# In[158]:


np.savetxt("AI5KF.csv", AI5KF, delimiter=",", header = 'SNP_ID')


# In[159]:


np.savetxt("AI6GF.csv", AI6GF, delimiter=",", header = 'SNP_ID')


# In[160]:


np.savetxt("AI6KF.csv", AI6KF, delimiter=",", header = 'SNP_ID')


# These tables don't work for a Venn diagram but to chek the SNPs in AI by individual and tissue with the help of the dictionaries (vlookup).

# ### Intersection of treatment

# We group the common SNPs in AI for gills and expressed in liver for each treatment.

# In[161]:


binom = pd.read_csv('SNPsID_AI_binomial_test_tissues_GF3KF6.csv')
binom.columns


# In[162]:


binom_1GF = binom['binom_1GF']
binom_2GF = binom['binom_2GF']
binom_3GF = binom['binom_3GF']
binom_4GF = binom['binom_4GF']
binom_5GF = binom['binom_5GF']
binom_6GF = binom['binom_6GF']

binom_1GS = binom['binom_1GS']
binom_2GS = binom['binom_2GS']
binom_3GS = binom['binom_3GS']
binom_4GS = binom['binom_4GS']
binom_5GS = binom['binom_5GS']
binom_6GS = binom['binom_6GS']

binom_1KF = binom['binom_1KF']
binom_2KF = binom['binom_2KF']
binom_3KF = binom['binom_3KF']
binom_4KF = binom['binom_4KF']
binom_5KF = binom['binom_5KF']
binom_6KF = binom['binom_6KF']

binom_2KS = binom['binom_2KS']
binom_3KS = binom['binom_3KS']
binom_4KS = binom['binom_4KS']
binom_5KS = binom['binom_5KS']
binom_6KS = binom['binom_6KS']


# In[163]:


from functools import reduce
AI_gillsF = reduce(np.intersect1d, (binom_1GF, binom_2GF, binom_3GF, binom_4GF, binom_5GF, binom_6GF))
AI_gillsS = reduce(np.intersect1d, (binom_1GS, binom_2GS, binom_3GS, binom_4GS, binom_5GS, binom_6GS))

AI_kidneyF = reduce(np.intersect1d, (binom_1KF, binom_2KF, binom_3KF, binom_4KF, binom_5KF, binom_6KF))
AI_kidneyS = reduce(np.intersect1d, (binom_2KS, binom_3KS, binom_4KS, binom_5KS, binom_6KS))

AI_gillsFS = np.intersect1d(AI_gillsF, AI_gillsS)
AI_kidneyFS = np.intersect1d(AI_kidneyF, AI_kidneyS)

AI_gk_F = np.intersect1d(AI_gillsF, AI_kidneyF)
AI_gk_S = np.intersect1d(AI_gillsS, AI_kidneyS)
AI_gF_kS = np.intersect1d(AI_gillsF, AI_kidneyS)
AI_gS_kF = np.intersect1d(AI_gillsS, AI_kidneyF)



# In[164]:


AI_treatment = pd.DataFrame(data = [AI_gillsF,
                                    AI_gillsS,
                                    AI_kidneyF,
                                    AI_kidneyS,
                                    AI_gillsFS,
                                    AI_kidneyFS,
                                    AI_gk_F,
                                    AI_gk_S,
                                    AI_gF_kS,
                                    AI_gS_kF]).T
AI_treatment = AI_treatment.rename(columns={0:'Gills FW',
                                            1:'Gills SW',
                                            2:'Kidney FW',
                                            3:'Kidney SW',
                                            4:'Gills',
                                            5:'Kidney',
                                            6:'Fresh water g&k',
                                            7:'Sea water g&k',
                                            8:'Gills FW & kidney SW',
                                            9:'Gills SW & kidney FW'})


# In[165]:


AI_treatment.head()


# In[166]:


AI_treatment.to_csv('AI_treatment_GF3KF6.csv')


# In[167]:


df_binom.to_csv('SNPs_AF_binom_GF3KF6.csv')


# Now you have all the input needed for the script vlookup_geneID

#!/usr/bin/env python
# coding: utf-8

# # ASE with data from two pseudogenomes

# This script contains the indications for data wrangling of two tables containing biallelic SNPs. These tables have been created from mapping to two different pseudogenomes. Each pseudogenome has been constructed from the haplotype of a different sample in different experimental groups. The objective is to merge the SNPs collected from each pseudogenome and to do the average of the counts in order to obtain an unbiased calculation of the reads in each gene isoform thus avoiding the mapping bias. 
# 
# In this script, the reference pseudogenome is the one created with the sample named GF3 (gills tissue exposed to fresh water corresponding to the individual 3) and the alternative pseudogenome is the one created with the sample named KF6 (kidney tissue, fresh water from individual 6).
# 
# In order to process the data, we will use the python libraries pandas and numpy.

# In[1]:


import pandas as pd
import numpy as np


# Read the data files

# In[4]:


GF3 = pd.read_csv('AD_GT_counts_bi_GF3.csv', low_memory=False)
KF6 = pd.read_csv('AD_GT_counts_bi_KF6.csv', low_memory=False)


# Merge the dataframes in order to achieve a reference SNP ID for a particular chromosome and position.
# 
# First, convert them into dataframes

# In[3]:


GF3 = pd.DataFrame(GF3)
KF6 = pd.DataFrame(KF6)


# Now collect the SNPs by chromosome and position into dataframes that will provide indexes

# In[4]:


GF3_SNPs = pd.DataFrame(GF3, columns = ['CHROM','POS','Gene.refGene','Func.refGene','ExonicFunc.refGene'])
KF6_SNPs = pd.DataFrame(KF6, columns = ['CHROM','POS','Gene.refGene','Func.refGene','ExonicFunc.refGene'])


# Check the common SNP between all four SNP tables

# In[5]:


SNPs_common = pd.merge(GF3_SNPs, KF6_SNPs, how = 'inner')
SNPs_common.head()
len(SNPs_common)


# Check all the SNPs that were called by the different mapping methods and not in common

# In[6]:


SNPs_out = pd.merge(GF3_SNPs, KF6_SNPs, how = 'outer')
SNPs_out.head()


# Explore the common SNPs in ref and pseudo to do the average: 

# In[7]:


len(SNPs_out)


# In[8]:


len(SNPs_out.CHROM.unique())


# In[9]:


len(SNPs_common.CHROM.unique())


# In[10]:


len(SNPs_common)


# Merge the two datafiles on the cordinates CHROM (chromosome) and POS (position) for retrieving the SNPs in common from the two mapping methods

# In[11]:


df = pd.merge(GF3,KF6, on = ('CHROM','POS'))
len(df)


# In[12]:


df.head()


# ## Compare the genotypes for assessing the counts on each allele

# Set the data for the mapping against the pseudogenome of sample GF3 because it hast the variants of the reference genotype. The SNPs called by this sample are set as reference and alternative alleles for all the other samples. 
# 
# The tables containing the SNPs from mapping to GF3 and KF6 contained only biallelic sites. However multiallelic possibilities can appear if one SNP was called for a different genotpe in each of the pseudogenomes. The multiallelic sites will be deleted from this dataframe to continue the analysis with bialleleic sites only.
# 

# ### Delete multiallelic sites

# Let's start droping the multiallelic. The criteria compares if the reference allele in the SNPs resulting from mapping against GF3 pseudogenome is different form the reference and the alternative alleles from the mapping against KF6 pseudogenome. This is later repeated for the alternative allele. Note that the SNPs mapped against GF3 pseudogenome have been pre-filtered for biallelic sites for all the samples.
# 
# Collect the indexes where the refrence allele in GF3 is different from both alleles in the other mappings or the alternative allele in GF3 is different from both alleles in the KF6 mapping.

# In[13]:


df_indexmultiGF3 = df[
              (((df['1GF_R_GF3.GT'] != df['1GF_R_KF6.GT']) & (df['1GF_R_GF3.GT'] != df['1GF_A_KF6.GT'])) |
               ((df['1GF_A_GF3.GT'] != df['1GF_R_KF6.GT']) & (df['1GF_A_GF3.GT'] != df['1GF_A_KF6.GT'])) | 
               ((df['2GF_R_GF3.GT'] != df['2GF_R_KF6.GT']) & (df['2GF_R_GF3.GT'] != df['2GF_A_KF6.GT'])) |
               ((df['2GF_A_GF3.GT'] != df['2GF_R_KF6.GT']) & (df['2GF_A_GF3.GT'] != df['2GF_A_KF6.GT'])) | 
               ((df['3GF_R_GF3.GT'] != df['3GF_R_KF6.GT']) & (df['3GF_R_GF3.GT'] != df['3GF_A_KF6.GT'])) |
               ((df['3GF_A_GF3.GT'] != df['3GF_R_KF6.GT']) & (df['3GF_A_GF3.GT'] != df['3GF_A_KF6.GT'])) | 
               ((df['4GF_R_GF3.GT'] != df['4GF_R_KF6.GT']) & (df['4GF_R_GF3.GT'] != df['4GF_A_KF6.GT'])) |
               ((df['4GF_A_GF3.GT'] != df['4GF_R_KF6.GT']) & (df['4GF_A_GF3.GT'] != df['4GF_A_KF6.GT'])) | 
               ((df['5GF_R_GF3.GT'] != df['5GF_R_KF6.GT']) & (df['5GF_R_GF3.GT'] != df['5GF_A_KF6.GT'])) |
               ((df['5GF_A_GF3.GT'] != df['5GF_R_KF6.GT']) & (df['5GF_A_GF3.GT'] != df['5GF_A_KF6.GT'])) |
               ((df['6GF_R_GF3.GT'] != df['6GF_R_KF6.GT']) & (df['6GF_R_GF3.GT'] != df['6GF_A_KF6.GT'])) |
               ((df['6GF_A_GF3.GT'] != df['6GF_R_KF6.GT']) & (df['6GF_A_GF3.GT'] != df['6GF_A_KF6.GT'])) | 
               ((df['1GS_R_GF3.GT'] != df['1GS_R_KF6.GT']) & (df['1GS_R_GF3.GT'] != df['1GS_A_KF6.GT'])) |
               ((df['1GS_A_GF3.GT'] != df['1GS_A_KF6.GT']) & (df['1GS_A_GF3.GT'] != df['1GS_A_KF6.GT'])) | 
               ((df['2GS_R_GF3.GT'] != df['2GS_R_KF6.GT']) & (df['2GS_R_GF3.GT'] != df['2GS_A_KF6.GT'])) |
               ((df['2GS_A_GF3.GT'] != df['2GS_R_KF6.GT']) & (df['2GS_A_GF3.GT'] != df['2GS_A_KF6.GT'])) |
               ((df['3GS_R_GF3.GT'] != df['3GS_R_KF6.GT']) & (df['3GS_R_GF3.GT'] != df['3GS_A_KF6.GT'])) |
               ((df['3GS_A_GF3.GT'] != df['3GS_R_KF6.GT']) & (df['3GS_A_GF3.GT'] != df['3GS_A_KF6.GT'])) | 
               ((df['4GS_R_GF3.GT'] != df['4GS_R_KF6.GT']) & (df['4GS_R_GF3.GT'] != df['4GS_A_KF6.GT'])) |
               ((df['4GS_A_GF3.GT'] != df['4GS_R_KF6.GT']) & (df['4GS_A_GF3.GT'] != df['4GS_A_KF6.GT'])) | 
               ((df['5GS_R_GF3.GT'] != df['5GS_R_KF6.GT']) & (df['5GS_R_GF3.GT'] != df['5GS_A_KF6.GT'])) |
               ((df['5GS_A_GF3.GT'] != df['5GS_R_KF6.GT']) & (df['5GS_A_GF3.GT'] != df['5GS_A_KF6.GT'])) | 
               ((df['6GS_R_GF3.GT'] != df['6GS_R_KF6.GT']) & (df['6GS_R_GF3.GT'] != df['6GS_A_KF6.GT'])) |
               ((df['6GS_A_GF3.GT'] != df['6GS_R_KF6.GT']) & (df['6GS_A_GF3.GT'] != df['6GS_A_KF6.GT'])) | 
               ((df['1KF_R_GF3.GT'] != df['1KF_R_KF6.GT']) & (df['1KF_R_GF3.GT'] != df['1KF_A_KF6.GT'])) |
               ((df['1KF_A_GF3.GT'] != df['1KF_R_KF6.GT']) & (df['1KF_A_GF3.GT'] != df['1KF_A_KF6.GT'])) | 
               ((df['2KF_R_GF3.GT'] != df['2KF_R_KF6.GT']) & (df['2KF_R_GF3.GT'] != df['2KF_A_KF6.GT'])) |
               ((df['2KF_A_GF3.GT'] != df['2KF_R_KF6.GT']) & (df['2KF_A_GF3.GT'] != df['2KF_A_KF6.GT'])) | 
               ((df['3KF_R_GF3.GT'] != df['3KF_R_KF6.GT']) & (df['3KF_R_GF3.GT'] != df['3KF_A_KF6.GT'])) |
               ((df['3KF_A_GF3.GT'] != df['3KF_R_KF6.GT']) & (df['3KF_A_GF3.GT'] != df['3KF_A_KF6.GT'])) | 
               ((df['4KF_R_GF3.GT'] != df['4KF_R_KF6.GT']) & (df['4KF_R_GF3.GT'] != df['4KF_A_KF6.GT'])) |
               ((df['4KF_A_GF3.GT'] != df['4KF_R_KF6.GT']) & (df['4KF_A_GF3.GT'] != df['4KF_A_KF6.GT'])) | 
               ((df['5KF_R_GF3.GT'] != df['5KF_R_KF6.GT']) & (df['5KF_R_GF3.GT'] != df['5KF_A_KF6.GT'])) |
               ((df['5KF_A_GF3.GT'] != df['5KF_R_KF6.GT']) & (df['5KF_A_GF3.GT'] != df['5KF_A_KF6.GT'])) | 
               ((df['6KF_R_GF3.GT'] != df['6KF_R_KF6.GT']) & (df['6KF_R_GF3.GT'] != df['6KF_A_KF6.GT'])) |
               ((df['6KF_A_GF3.GT'] != df['6KF_R_KF6.GT']) & (df['6KF_A_GF3.GT'] != df['6KF_A_KF6.GT'])) | 
               ((df['2KS_R_GF3.GT'] != df['2KS_R_KF6.GT']) & (df['2KS_R_GF3.GT'] != df['2KS_A_KF6.GT'])) |
               ((df['2KS_A_GF3.GT'] != df['2KS_R_KF6.GT']) & (df['2KS_A_GF3.GT'] != df['2KS_A_KF6.GT'])) | 
               ((df['3KS_R_GF3.GT'] != df['3KS_R_KF6.GT']) & (df['3KS_R_GF3.GT'] != df['3KS_A_KF6.GT'])) |
               ((df['3KS_A_GF3.GT'] != df['3KS_R_KF6.GT']) & (df['3KS_A_GF3.GT'] != df['3KS_A_KF6.GT'])) | 
               ((df['4KS_R_GF3.GT'] != df['4KS_R_KF6.GT']) & (df['4KS_R_GF3.GT'] != df['4KS_A_KF6.GT'])) |
               ((df['4KS_A_GF3.GT'] != df['4KS_R_KF6.GT']) & (df['4KS_A_GF3.GT'] != df['4KS_A_KF6.GT'])) | 
               ((df['5KS_R_GF3.GT'] != df['5KS_R_KF6.GT']) & (df['5KS_R_GF3.GT'] != df['5KS_A_KF6.GT'])) |
               ((df['5KS_A_GF3.GT'] != df['5KS_R_KF6.GT']) & (df['5KS_A_GF3.GT'] != df['5KS_A_KF6.GT'])) | 
               ((df['6KS_R_GF3.GT'] != df['6KS_R_KF6.GT']) & (df['6KS_R_GF3.GT'] != df['6KS_A_KF6.GT'])) |
               ((df['6KS_A_GF3.GT'] != df['6KS_R_KF6.GT']) & (df['6KS_A_GF3.GT'] != df['6KS_A_KF6.GT'])))].index


# In[14]:


df_bi = df.drop(index = df_indexmultiGF3)
len(df_bi)


# In[15]:


df_bi.head(20)


# To decrease the load of the big file, select the columns that are to be analysed

# In[16]:


for column in df_bi.columns :
    print (column)


# In[17]:


df_backup = df_bi


# In[18]:


df_bi = df_bi[['CHROM',
'POS',
'Gene.refGene_x',
'Func.refGene_x',
'ExonicFunc.refGene_x',
'AF_x',
'1GF_R_GF3.AD',
'1GF_A_GF3.AD',
'1GF_R_GF3.GT',
'1GF_A_GF3.GT',
'1GS_R_GF3.AD',
'1GS_A_GF3.AD',
'1GS_R_GF3.GT',
'1GS_A_GF3.GT',
'1KF_R_GF3.AD',
'1KF_A_GF3.AD',
'1KF_R_GF3.GT',
'1KF_A_GF3.GT',
'1KS_R_GF3.AD',
'1KS_A_GF3.AD',
'1KS_R_GF3.GT',
'1KS_A_GF3.GT',
'2GF_R_GF3.AD',
'2GF_A_GF3.AD',
'2GF_R_GF3.GT',
'2GF_A_GF3.GT',
'2GS_R_GF3.AD',
'2GS_A_GF3.AD',
'2GS_R_GF3.GT',
'2GS_A_GF3.GT',
'2KF_R_GF3.AD',
'2KF_A_GF3.AD',
'2KF_R_GF3.GT',
'2KF_A_GF3.GT',
'2KS_R_GF3.AD',
'2KS_A_GF3.AD',
'2KS_R_GF3.GT',
'2KS_A_GF3.GT',
'3GF_R_GF3.AD',
'3GF_A_GF3.AD',
'3GF_R_GF3.GT',
'3GF_A_GF3.GT',
'3GS_R_GF3.AD',
'3GS_A_GF3.AD',
'3GS_R_GF3.GT',
'3GS_A_GF3.GT',
'3KF_R_GF3.AD',
'3KF_A_GF3.AD',
'3KF_R_GF3.GT',
'3KF_A_GF3.GT',
'3KS_R_GF3.AD',
'3KS_A_GF3.AD',
'3KS_R_GF3.GT',
'3KS_A_GF3.GT',
'4GF_R_GF3.AD',
'4GF_A_GF3.AD',
'4GF_R_GF3.GT',
'4GF_A_GF3.GT',
'4GS_R_GF3.AD',
'4GS_A_GF3.AD',
'4GS_R_GF3.GT',
'4GS_A_GF3.GT',
'4KF_R_GF3.AD',
'4KF_A_GF3.AD',
'4KF_R_GF3.GT',
'4KF_A_GF3.GT',
'4KS_R_GF3.AD',
'4KS_A_GF3.AD',
'4KS_R_GF3.GT',
'4KS_A_GF3.GT',
'5GF_R_GF3.AD',
'5GF_A_GF3.AD',
'5GF_R_GF3.GT',
'5GF_A_GF3.GT',
'5GS_R_GF3.AD',
'5GS_A_GF3.AD',
'5GS_R_GF3.GT',
'5GS_A_GF3.GT',
'5KF_R_GF3.AD',
'5KF_A_GF3.AD',
'5KF_R_GF3.GT',
'5KF_A_GF3.GT',
'5KS_R_GF3.AD',
'5KS_A_GF3.AD',
'5KS_R_GF3.GT',
'5KS_A_GF3.GT',
'6GF_R_GF3.AD',
'6GF_A_GF3.AD',
'6GF_R_GF3.GT',
'6GF_A_GF3.GT',
'6GS_R_GF3.AD',
'6GS_A_GF3.AD',
'6GS_R_GF3.GT',
'6GS_A_GF3.GT',
'6KF_R_GF3.AD',
'6KF_A_GF3.AD',
'6KF_R_GF3.GT',
'6KF_A_GF3.GT',
'6KS_R_GF3.AD',
'6KS_A_GF3.AD',
'6KS_R_GF3.GT',
'6KS_A_GF3.GT',
'Gene.refGene_y',
'Func.refGene_y',
'ExonicFunc.refGene_y',
'AF_y',
'1GF_R_KF6.AD',
'1GF_A_KF6.AD',
'1GF_R_KF6.GT',
'1GF_A_KF6.GT',
'1GS_R_KF6.AD',
'1GS_A_KF6.AD',
'1GS_R_KF6.GT',
'1GS_A_KF6.GT',
'1KF_R_KF6.AD',
'1KF_A_KF6.AD',
'1KF_R_KF6.GT',
'1KF_A_KF6.GT',
'1KS_R_KF6.AD',
'1KS_A_KF6.AD',
'1KS_R_KF6.GT',
'1KS_A_KF6.GT',
'2GF_R_KF6.AD',
'2GF_A_KF6.AD',
'2GF_R_KF6.GT',
'2GF_A_KF6.GT',
'2GS_R_KF6.AD',
'2GS_A_KF6.AD',
'2GS_R_KF6.GT',
'2GS_A_KF6.GT',
'2KF_R_KF6.AD',
'2KF_A_KF6.AD',
'2KF_R_KF6.GT',
'2KF_A_KF6.GT',
'2KS_R_KF6.AD',
'2KS_A_KF6.AD',
'2KS_R_KF6.GT',
'2KS_A_KF6.GT',
'3GF_R_KF6.AD',
'3GF_A_KF6.AD',
'3GF_R_KF6.GT',
'3GF_A_KF6.GT',
'3GS_R_KF6.AD',
'3GS_A_KF6.AD',
'3GS_R_KF6.GT',
'3GS_A_KF6.GT',
'3KF_R_KF6.AD',
'3KF_A_KF6.AD',
'3KF_R_KF6.GT',
'3KF_A_KF6.GT',
'3KS_R_KF6.AD',
'3KS_A_KF6.AD',
'3KS_R_KF6.GT',
'3KS_A_KF6.GT',
'4GF_R_KF6.AD',
'4GF_A_KF6.AD',
'4GF_R_KF6.GT',
'4GF_A_KF6.GT',
'4GS_R_KF6.AD',
'4GS_A_KF6.AD',
'4GS_R_KF6.GT',
'4GS_A_KF6.GT',
'4KF_R_KF6.AD',
'4KF_A_KF6.AD',
'4KF_R_KF6.GT',
'4KF_A_KF6.GT',
'4KS_R_KF6.AD',
'4KS_A_KF6.AD',
'4KS_R_KF6.GT',
'4KS_A_KF6.GT',
'5GF_R_KF6.AD',
'5GF_A_KF6.AD',
'5GF_R_KF6.GT',
'5GF_A_KF6.GT',
'5GS_R_KF6.AD',
'5GS_A_KF6.AD',
'5GS_R_KF6.GT',
'5GS_A_KF6.GT',
'5KF_R_KF6.AD',
'5KF_A_KF6.AD',
'5KF_R_KF6.GT',
'5KF_A_KF6.GT',
'5KS_R_KF6.AD',
'5KS_A_KF6.AD',
'5KS_R_KF6.GT',
'5KS_A_KF6.GT',
'6GF_R_KF6.AD',
'6GF_A_KF6.AD',
'6GF_R_KF6.GT',
'6GF_A_KF6.GT',
'6GS_R_KF6.AD',
'6GS_A_KF6.AD',
'6GS_R_KF6.GT',
'6GS_A_KF6.GT',
'6KF_R_KF6.AD',
'6KF_A_KF6.AD',
'6KF_R_KF6.GT',
'6KF_A_KF6.GT',
'6KS_R_KF6.AD',
'6KS_A_KF6.AD',
'6KS_R_KF6.GT',
'6KS_A_KF6.GT']]


# Write a backup copy of the datafile with the biallelic SNPs.

# In[19]:


df_bi.to_csv('Big_df_2map.csv')


# ## Average for alleles

# When using a big dataset sometimes the kernel may block. Therefore, for a fresh start charge the file containing the biallelic sites

# In[20]:


import pandas as pd
import numpy as np


# In[21]:


df_bi = pd.read_csv('Big_df_2map.csv', low_memory = False)


# Now let's make the average on the counts for the reference and the alternative alleles. For each sample, select the genotype of the reference allele established by the reads mapped against GF3 and compare it to both alternative and reference genotypes of the values in the mapping with KF6 genome. Then proceed to make the average of the counts. The possible cases are the next:
# 
# 1.- Both homozigots for the reference allele: In that case set the same genotype to both alleles and sum all the counts from both sites, making the average for two individuals.
# 
# 2.- Reference and alternative allele in the SNP from the GF3 pseudogenome are placed in the same position for the SNP from the KF6 pseudogenome. Set the genotype to the reference and alternative alleles to the ones from the SNP in GF3 and make the average of the counts for reference and alternative genome
# 
# 3.- The reference allele in the reference pseudogenome is called as alternative allele in the alternative pseudogenome. In that case, set the genotype to the reference allele and make the average of the counts for the corresponding variant.
# 
# 4.- While the mapping against the reference genome resulted in an homozygot for the reference allele, the mapping against the alternative genome resulted in an heterozygot with the alternative allele in the alternative position. The genotype of the reference allele will be set from genotype of the reference allele in the reference mapping and the alternative allele will be set from the genotype of the alternative allele in the alternative mapping. The counts of each allele will be set accordingly to this distribution, thus dividing the counts by two individuals in the reference allele. 
# 
# 5.- The mapping against the reference genome resulted in an homozygot for the reference allele and the mapping against the alternative genome is heterozygot. In that case, the genotype of the alternative allele is set as the reference allele on the mapping against the alternative genome. The genotypes must be set accordingly, being the reference allele the one found in the reference genome and the alternative allele the one found in the alternative genome. The average of the counts will follow this pattern, being the average in the erference allele.
# 
# 6.- The mapping against the reference genome produced an heterozygot with reference and alternative alleles. The mapping against the alternative pseudogenome produced an homozygot for the reference allele. In that case, reference and alternative allele must be set according to the reference pseudogenome and the counts must be averaged for both individuals in the alternative allele.
# 
# 7.- The mapping against the reference genome produced an heterozygot with reference and alternative alleles. The mapping against the alternative pseudogenome produced an homozygot for the alternative allele. In that case, reference and alternative allele must be set according to the reference pseudogenome and the counts must be averaged for both individuals in the alternative allele.
# 
# This analysis is repeated for each sample
# 

# In[22]:


GF1_GF3KF6_R_GT = []
GF1_GF3KF6_R_AD = []
GF1_GF3KF6_A_GT = []
GF1_GF3KF6_A_AD = []
Homozygous = []


# In[23]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])'''
    '''Case heterozygots'''
    if ((row['1GF_R_GF3.GT'] == row['1GF_A_GF3.GT']) and (row['1GF_R_KF6.GT'] == row['1GF_A_KF6.GT']) and (row['1GF_R_GF3.GT'] == row['1GF_R_KF6.GT']) and (row['1GF_A_GF3.GT'] == row['1GF_A_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = ((row['1GF_R_GF3.AD'] + row['1GF_R_KF6.AD'] + row['1GF_A_GF3.AD'] + row['1GF_A_KF6.AD'])/2)
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_R_KF6.GT'])
            y = 0
            GF1_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Case 1 Both homozgots',i)'''
    elif ((row['1GF_R_GF3.GT'] == row['1GF_R_KF6.GT']) and (row['1GF_A_GF3.GT'] == row['1GF_A_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = (row['1GF_R_GF3.AD'] + row['1GF_R_KF6.AD'])/2
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_A_GF3.GT'])
            y = (row['1GF_A_GF3.AD'] + row['1GF_A_KF6.AD'])/2
            GF1_GF3KF6_A_AD.append(y)
            '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])
            print('Case 2 Same situation Ref and Alt in both for row',i)'''
    elif ((row['1GF_R_GF3.GT'] == row['1GF_A_KF6.GT']) and (row['1GF_A_GF3.GT'] == row['1GF_R_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = (row['1GF_R_GF3.AD'] + row['1GF_A_KF6.AD'])/2
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_A_GF3.GT'])
            y = (row['1GF_A_GF3.AD'] + row['1GF_R_KF6.AD'])/2
            GF1_GF3KF6_A_AD.append(y)
            '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])
            print('Case 3 Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['1GF_R_GF3.GT'] == row['1GF_A_GF3.GT']) and (row['1GF_R_KF6.GT'] != row['1GF_A_KF6.GT']) and (row['1GF_R_GF3.GT'] == row['1GF_R_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = (row['1GF_R_GF3.AD'] + row['1GF_A_GF3.AD'] + row['1GF_R_KF6.AD'])/2 # counts from 3 cells from which one is 0, from 2 individuals
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_A_KF6.GT'])
            y = row['1GF_A_KF6.AD']
            GF1_GF3KF6_A_AD.append(y)
            '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])
            print('Case 4 First homozygot, second with alternative allele in second position',i)'''
    elif ((row['1GF_R_GF3.GT'] == row['1GF_A_GF3.GT']) and (row['1GF_R_KF6.GT'] != row['1GF_A_KF6.GT']) and (row['1GF_R_GF3.GT'] == row['1GF_A_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = (row['1GF_R_GF3.AD'] + row['1GF_A_GF3.AD'] + row['1GF_A_KF6.AD'])/2 # counts from 3 cells from which one is 0, from 2 individuals
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_R_KF6.GT'])
            y = row['1GF_R_KF6.AD']
            GF1_GF3KF6_A_AD.append(y)
            '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])
            print('Case 5 First homozygot, second with alternative allele in first position of the second genome',i)'''
    elif ((row['1GF_R_GF3.GT'] != row['1GF_A_GF3.GT']) and (row['1GF_R_KF6.GT'] == row['1GF_A_KF6.GT']) and (row['1GF_R_GF3.GT'] == row['1GF_R_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = row['1GF_R_GF3.AD'] 
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_A_GF3.GT'])
            y = (row['1GF_A_GF3.AD'] + row['1GF_R_KF6.AD'] + row['1GF_A_KF6.AD'])/2  # counts from 3 cells from which one is 0, from 2 individuals
            GF1_GF3KF6_A_AD.append(y)
            '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])
            print('Case 6 Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['1GF_R_GF3.GT'] != row['1GF_A_GF3.GT']) and (row['1GF_R_KF6.GT'] == row['1GF_A_KF6.GT']) and (row['1GF_R_GF3.GT'] == row['1GF_A_KF6.GT'])):
            GF1_GF3KF6_R_GT.append(row['1GF_R_GF3.GT'])
            x = row['1GF_R_GF3.AD']
            GF1_GF3KF6_R_AD.append(x)
            GF1_GF3KF6_A_GT.append(row['1GF_R_KF6.GT'])
            y = (row['1GF_A_GF3.AD'] + row['1GF_R_KF6.AD'] + row['1GF_A_KF6.AD'])/2  # counts from 3 cells from which one is 0, from 2 individuals
            GF1_GF3KF6_A_AD.append(y)
            '''print ('index',i,'length',len(GF1_GF3KF6_R_GT),'Ref ',row['1GF_R_GF3.GT'],row['1GF_R_KF6.GT'],'Alt ',row['1GF_A_GF3.GT'],row['1GF_A_KF6.GT'])
            print('Case 7 Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GF1_GF3KF6_R_GT.append('error')
        GF1_GF3KF6_R_AD.append('error')
        GF1_GF3KF6_A_GT.append('error')
        GF1_GF3KF6_A_AD.append('error')


# In[24]:


# print(GF1_GF3KF6_R_AD)


# In[25]:


df_backup['1GF_R_GF3KF6.GT'] = GF1_GF3KF6_R_GT
df_backup['1GF_A_GF3KF6.GT'] = GF1_GF3KF6_A_GT
df_backup['1GF_R_GF3KF6.AD'] = GF1_GF3KF6_R_AD
df_backup['1GF_A_GF3KF6.AD'] = GF1_GF3KF6_A_AD


# In[26]:


df_backup.head()


# Now repeat for all other samples
# GF2

# In[27]:


GF2_GF3KF6_R_GT = []
GF2_GF3KF6_R_AD = []
GF2_GF3KF6_A_GT = []
GF2_GF3KF6_A_AD = []


# In[28]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GF2_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['2GF_R_GF3.GT'] == row['2GF_A_GF3.GT']) and (row['2GF_R_KF6.GT'] == row['2GF_A_KF6.GT']) and (row['2GF_R_GF3.GT'] == row['2GF_R_KF6.GT']) and (row['2GF_A_GF3.GT'] == row['2GF_A_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = (row['2GF_R_GF3.AD'] + row['2GF_R_KF6.AD'] + row['2GF_A_GF3.AD'] + row['2GF_A_KF6.AD'])/2
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_R_KF6.GT'])
            y = 0
            GF2_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['2GF_R_GF3.GT'] == row['2GF_A_KF6.GT']) and (row['2GF_A_GF3.GT'] == row['2GF_R_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = (row['2GF_R_GF3.AD'] + row['2GF_A_KF6.AD'])/2
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_A_GF3.GT'])
            y = (row['2GF_A_GF3.AD'] + row['2GF_R_KF6.AD'])/2
            GF2_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['2GF_R_GF3.GT'] == row['2GF_R_KF6.GT']) and (row['2GF_A_GF3.GT'] == row['2GF_A_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = (row['2GF_R_GF3.AD'] + row['2GF_R_KF6.AD'])/2
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_A_GF3.GT'])
            y = (row['2GF_A_GF3.AD'] + row['2GF_A_KF6.AD'])/2
            GF2_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['2GF_R_GF3.GT'] == row['2GF_A_GF3.GT']) and (row['2GF_R_KF6.GT'] != row['2GF_A_KF6.GT']) and (row['2GF_R_GF3.GT'] == row['2GF_R_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = (row['2GF_R_GF3.AD'] + row['2GF_A_GF3.AD'] + row['2GF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_A_KF6.GT'])
            y = row['2GF_A_KF6.AD']
            GF2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['2GF_R_GF3.GT'] == row['2GF_A_GF3.GT']) and (row['2GF_R_KF6.GT'] != row['2GF_A_KF6.GT']) and (row['2GF_R_GF3.GT'] == row['2GF_A_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = (row['2GF_R_GF3.AD'] + row['2GF_A_GF3.AD'] + row['2GF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_R_KF6.GT'])
            y = row['2GF_R_KF6.AD']
            GF2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['2GF_R_GF3.GT'] != row['2GF_A_GF3.GT']) and (row['2GF_R_KF6.GT'] == row['2GF_A_KF6.GT']) and (row['2GF_R_GF3.GT'] == row['2GF_R_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = row['2GF_R_GF3.AD'] 
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_A_GF3.GT'])
            y = (row['2GF_A_GF3.AD'] + row['2GF_R_KF6.AD'] + row['2GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['2GF_R_GF3.GT'] != row['2GF_A_GF3.GT']) and (row['2GF_R_KF6.GT'] == row['2GF_A_KF6.GT']) and (row['2GF_R_GF3.GT'] == row['2GF_A_KF6.GT'])):
            GF2_GF3KF6_R_GT.append(row['2GF_R_GF3.GT'])
            x = row['2GF_R_GF3.AD']
            GF2_GF3KF6_R_AD.append(x)
            GF2_GF3KF6_A_GT.append(row['2GF_R_KF6.GT'])
            y = (row['2GF_A_GF3.AD'] + row['2GF_R_KF6.AD'] + row['2GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GF2_GF3KF6_R_GT.append('error')
        GF2_GF3KF6_R_AD.append('error')
        GF2_GF3KF6_A_GT.append('error')
        GF2_GF3KF6_A_AD.append('error')


# In[29]:


df_backup['2GF_R_GF3KF6.GT'] = GF2_GF3KF6_R_GT
df_backup['2GF_A_GF3KF6.GT'] = GF2_GF3KF6_A_GT
df_backup['2GF_R_GF3KF6.AD'] = GF2_GF3KF6_R_AD
df_backup['2GF_A_GF3KF6.AD'] = GF2_GF3KF6_A_AD


# GF3

# In[30]:


GF3_GF3KF6_R_GT = []
GF3_GF3KF6_R_AD = []
GF3_GF3KF6_A_GT = []
GF3_GF3KF6_A_AD = []


# In[31]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GF3_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_R_KF6.GT']) and (row['3GF_A_KF6.GT'] == row['3GF_A_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = (row['3GF_R_KF6.AD'] + row['3GF_R_KF6.AD'] + row['3GF_A_KF6.AD'] + row['3GF_A_KF6.AD'])/2
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_R_KF6.GT'])
            y = 0
            GF3_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_A_KF6.GT'] == row['3GF_R_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = (row['3GF_R_KF6.AD'] + row['3GF_A_KF6.AD'])/2
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_A_KF6.GT'])
            y = (row['3GF_A_KF6.AD'] + row['3GF_R_KF6.AD'])/2
            GF3_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['3GF_R_KF6.GT'] == row['3GF_R_KF6.GT']) and (row['3GF_A_KF6.GT'] == row['3GF_A_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = (row['3GF_R_KF6.AD'] + row['3GF_R_KF6.AD'])/2
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_A_KF6.GT'])
            y = (row['3GF_A_KF6.AD'] + row['3GF_A_KF6.AD'])/2
            GF3_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] != row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_R_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = (row['3GF_R_KF6.AD'] + row['3GF_A_KF6.AD'] + row['3GF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_A_KF6.GT'])
            y = row['3GF_A_KF6.AD']
            GF3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] != row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = (row['3GF_R_KF6.AD'] + row['3GF_A_KF6.AD'] + row['3GF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_R_KF6.GT'])
            y = row['3GF_R_KF6.AD']
            GF3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['3GF_R_KF6.GT'] != row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_R_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = row['3GF_R_KF6.AD'] 
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_A_KF6.GT'])
            y = (row['3GF_A_KF6.AD'] + row['3GF_R_KF6.AD'] + row['3GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['3GF_R_KF6.GT'] != row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT']) and (row['3GF_R_KF6.GT'] == row['3GF_A_KF6.GT'])):
            GF3_GF3KF6_R_GT.append(row['3GF_R_KF6.GT'])
            x = row['3GF_R_KF6.AD']
            GF3_GF3KF6_R_AD.append(x)
            GF3_GF3KF6_A_GT.append(row['3GF_R_KF6.GT'])
            y = (row['3GF_A_KF6.AD'] + row['3GF_R_KF6.AD'] + row['3GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GF3_GF3KF6_R_GT.append('error')
        GF3_GF3KF6_R_AD.append('error')
        GF3_GF3KF6_A_GT.append('error')
        GF3_GF3KF6_A_AD.append('error')


# In[32]:


df_backup['3GF_R_GF3KF6.GT'] = GF3_GF3KF6_R_GT
df_backup['3GF_A_GF3KF6.GT'] = GF3_GF3KF6_A_GT
df_backup['3GF_R_GF3KF6.AD'] = GF3_GF3KF6_R_AD
df_backup['3GF_A_GF3KF6.AD'] = GF3_GF3KF6_A_AD


# GF4

# In[33]:


GF4_GF3KF6_R_GT = []
GF4_GF3KF6_R_AD = []
GF4_GF3KF6_A_GT = []
GF4_GF3KF6_A_AD = []


# In[34]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GF4_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['4GF_R_GF3.GT'] == row['4GF_A_GF3.GT']) and (row['4GF_R_KF6.GT'] == row['4GF_A_KF6.GT']) and (row['4GF_R_GF3.GT'] == row['4GF_R_KF6.GT']) and (row['4GF_A_GF3.GT'] == row['4GF_A_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = (row['4GF_R_GF3.AD'] + row['4GF_R_KF6.AD'] + row['4GF_A_GF3.AD'] + row['4GF_A_KF6.AD'])/2
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_R_KF6.GT'])
            y = 0
            GF4_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['4GF_R_GF3.GT'] == row['4GF_A_KF6.GT']) and (row['4GF_A_GF3.GT'] == row['4GF_R_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = (row['4GF_R_GF3.AD'] + row['4GF_A_KF6.AD'])/2
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_A_GF3.GT'])
            y = (row['4GF_A_GF3.AD'] + row['4GF_R_KF6.AD'])/2
            GF4_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['4GF_R_GF3.GT'] == row['4GF_R_KF6.GT']) and (row['4GF_A_GF3.GT'] == row['4GF_A_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = (row['4GF_R_GF3.AD'] + row['4GF_R_KF6.AD'])/2
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_A_GF3.GT'])
            y = (row['4GF_A_GF3.AD'] + row['4GF_A_KF6.AD'])/2
            GF4_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['4GF_R_GF3.GT'] == row['4GF_A_GF3.GT']) and (row['4GF_R_KF6.GT'] != row['4GF_A_KF6.GT']) and (row['4GF_R_GF3.GT'] == row['4GF_R_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = (row['4GF_R_GF3.AD'] + row['4GF_A_GF3.AD'] + row['4GF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_A_KF6.GT'])
            y = row['4GF_A_KF6.AD']
            GF4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['4GF_R_GF3.GT'] == row['4GF_A_GF3.GT']) and (row['4GF_R_KF6.GT'] != row['4GF_A_KF6.GT']) and (row['4GF_R_GF3.GT'] == row['4GF_A_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = (row['4GF_R_GF3.AD'] + row['4GF_A_GF3.AD'] + row['4GF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_R_KF6.GT'])
            y = row['4GF_R_KF6.AD']
            GF4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['4GF_R_GF3.GT'] != row['4GF_A_GF3.GT']) and (row['4GF_R_KF6.GT'] == row['4GF_A_KF6.GT']) and (row['4GF_R_GF3.GT'] == row['4GF_R_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = row['4GF_R_GF3.AD'] 
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_A_GF3.GT'])
            y = (row['4GF_A_GF3.AD'] + row['4GF_R_KF6.AD'] + row['4GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['4GF_R_GF3.GT'] != row['4GF_A_GF3.GT']) and (row['4GF_R_KF6.GT'] == row['4GF_A_KF6.GT']) and (row['4GF_R_GF3.GT'] == row['4GF_A_KF6.GT'])):
            GF4_GF3KF6_R_GT.append(row['4GF_R_GF3.GT'])
            x = row['4GF_R_GF3.AD']
            GF4_GF3KF6_R_AD.append(x)
            GF4_GF3KF6_A_GT.append(row['4GF_R_KF6.GT'])
            y = (row['4GF_A_GF3.AD'] + row['4GF_R_KF6.AD'] + row['4GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GF4_GF3KF6_R_GT.append('error')
        GF4_GF3KF6_R_AD.append('error')
        GF4_GF3KF6_A_GT.append('error')
        GF4_GF3KF6_A_AD.append('error')


# In[35]:


df_backup['4GF_R_GF3KF6.GT'] = GF4_GF3KF6_R_GT
df_backup['4GF_A_GF3KF6.GT'] = GF4_GF3KF6_A_GT
df_backup['4GF_R_GF3KF6.AD'] = GF4_GF3KF6_R_AD
df_backup['4GF_A_GF3KF6.AD'] = GF4_GF3KF6_A_AD


# GF5

# In[36]:


GF5_GF3KF6_R_GT = []
GF5_GF3KF6_R_AD = []
GF5_GF3KF6_A_GT = []
GF5_GF3KF6_A_AD = []


# In[37]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GF5_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['5GF_R_GF3.GT'] == row['5GF_A_GF3.GT']) and (row['5GF_R_KF6.GT'] == row['5GF_A_KF6.GT']) and (row['5GF_R_GF3.GT'] == row['5GF_R_KF6.GT']) and (row['5GF_A_GF3.GT'] == row['5GF_A_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = (row['5GF_R_GF3.AD'] + row['5GF_R_KF6.AD'] + row['5GF_A_GF3.AD'] + row['5GF_A_KF6.AD'])/2
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_R_KF6.GT'])
            y = 0
            GF5_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['5GF_R_GF3.GT'] == row['5GF_A_KF6.GT']) and (row['5GF_A_GF3.GT'] == row['5GF_R_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = (row['5GF_R_GF3.AD'] + row['5GF_A_KF6.AD'])/2
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_A_GF3.GT'])
            y = (row['5GF_A_GF3.AD'] + row['5GF_R_KF6.AD'])/2
            GF5_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['5GF_R_GF3.GT'] == row['5GF_R_KF6.GT']) and (row['5GF_A_GF3.GT'] == row['5GF_A_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = (row['5GF_R_GF3.AD'] + row['5GF_R_KF6.AD'])/2
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_A_GF3.GT'])
            y = (row['5GF_A_GF3.AD'] + row['5GF_A_KF6.AD'])/2
            GF5_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['5GF_R_GF3.GT'] == row['5GF_A_GF3.GT']) and (row['5GF_R_KF6.GT'] != row['5GF_A_KF6.GT']) and (row['5GF_R_GF3.GT'] == row['5GF_R_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = (row['5GF_R_GF3.AD'] + row['5GF_A_GF3.AD'] + row['5GF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_A_KF6.GT'])
            y = row['5GF_A_KF6.AD']
            GF5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['5GF_R_GF3.GT'] == row['5GF_A_GF3.GT']) and (row['5GF_R_KF6.GT'] != row['5GF_A_KF6.GT']) and (row['5GF_R_GF3.GT'] == row['5GF_A_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = (row['5GF_R_GF3.AD'] + row['5GF_A_GF3.AD'] + row['5GF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_R_KF6.GT'])
            y = row['5GF_R_KF6.AD']
            GF5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['5GF_R_GF3.GT'] != row['5GF_A_GF3.GT']) and (row['5GF_R_KF6.GT'] == row['5GF_A_KF6.GT']) and (row['5GF_R_GF3.GT'] == row['5GF_R_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = row['5GF_R_GF3.AD'] 
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_A_GF3.GT'])
            y = (row['5GF_A_GF3.AD'] + row['5GF_R_KF6.AD'] + row['5GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['5GF_R_GF3.GT'] != row['5GF_A_GF3.GT']) and (row['5GF_R_KF6.GT'] == row['5GF_A_KF6.GT']) and (row['5GF_R_GF3.GT'] == row['5GF_A_KF6.GT'])):
            GF5_GF3KF6_R_GT.append(row['5GF_R_GF3.GT'])
            x = row['5GF_R_GF3.AD']
            GF5_GF3KF6_R_AD.append(x)
            GF5_GF3KF6_A_GT.append(row['5GF_R_KF6.GT'])
            y = (row['5GF_A_GF3.AD'] + row['5GF_R_KF6.AD'] + row['5GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GF5_GF3KF6_R_GT.append('error')
        GF5_GF3KF6_R_AD.append('error')
        GF5_GF3KF6_A_GT.append('error')
        GF5_GF3KF6_A_AD.append('error')


# In[38]:


df_backup['5GF_R_GF3KF6.GT'] = GF5_GF3KF6_R_GT
df_backup['5GF_A_GF3KF6.GT'] = GF5_GF3KF6_A_GT
df_backup['5GF_R_GF3KF6.AD'] = GF5_GF3KF6_R_AD
df_backup['5GF_A_GF3KF6.AD'] = GF5_GF3KF6_A_AD


# GF6

# In[39]:


GF6_GF3KF6_R_GT = []
GF6_GF3KF6_R_AD = []
GF6_GF3KF6_A_GT = []
GF6_GF3KF6_A_AD = []


# In[40]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GF6_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['6GF_R_GF3.GT'] == row['6GF_A_GF3.GT']) and (row['6GF_R_KF6.GT'] == row['6GF_A_KF6.GT']) and (row['6GF_R_GF3.GT'] == row['6GF_R_KF6.GT']) and (row['6GF_A_GF3.GT'] == row['6GF_A_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = (row['6GF_R_GF3.AD'] + row['6GF_R_KF6.AD'] + row['6GF_A_GF3.AD'] + row['6GF_A_KF6.AD'])/2
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_R_KF6.GT'])
            y = 0
            GF6_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['6GF_R_GF3.GT'] == row['6GF_A_KF6.GT']) and (row['6GF_A_GF3.GT'] == row['6GF_R_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = (row['6GF_R_GF3.AD'] + row['6GF_A_KF6.AD'])/2
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_A_GF3.GT'])
            y = (row['6GF_A_GF3.AD'] + row['6GF_R_KF6.AD'])/2
            GF6_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['6GF_R_GF3.GT'] == row['6GF_R_KF6.GT']) and (row['6GF_A_GF3.GT'] == row['6GF_A_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = (row['6GF_R_GF3.AD'] + row['6GF_R_KF6.AD'])/2
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_A_GF3.GT'])
            y = (row['6GF_A_GF3.AD'] + row['6GF_A_KF6.AD'])/2
            GF6_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['6GF_R_GF3.GT'] == row['6GF_A_GF3.GT']) and (row['6GF_R_KF6.GT'] != row['6GF_A_KF6.GT']) and (row['6GF_R_GF3.GT'] == row['6GF_R_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = (row['6GF_R_GF3.AD'] + row['6GF_A_GF3.AD'] + row['6GF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_A_KF6.GT'])
            y = row['6GF_A_KF6.AD']
            GF6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['6GF_R_GF3.GT'] == row['6GF_A_GF3.GT']) and (row['6GF_R_KF6.GT'] != row['6GF_A_KF6.GT']) and (row['6GF_R_GF3.GT'] == row['6GF_A_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = (row['6GF_R_GF3.AD'] + row['6GF_A_GF3.AD'] + row['6GF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_R_KF6.GT'])
            y = row['6GF_R_KF6.AD']
            GF6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['6GF_R_GF3.GT'] != row['6GF_A_GF3.GT']) and (row['6GF_R_KF6.GT'] == row['6GF_A_KF6.GT']) and (row['6GF_R_GF3.GT'] == row['6GF_R_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = row['6GF_R_GF3.AD'] 
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_A_GF3.GT'])
            y = (row['6GF_A_GF3.AD'] + row['6GF_R_KF6.AD'] + row['6GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['6GF_R_GF3.GT'] != row['6GF_A_GF3.GT']) and (row['6GF_R_KF6.GT'] == row['6GF_A_KF6.GT']) and (row['6GF_R_GF3.GT'] == row['6GF_A_KF6.GT'])):
            GF6_GF3KF6_R_GT.append(row['6GF_R_GF3.GT'])
            x = row['6GF_R_GF3.AD']
            GF6_GF3KF6_R_AD.append(x)
            GF6_GF3KF6_A_GT.append(row['6GF_R_KF6.GT'])
            y = (row['6GF_A_GF3.AD'] + row['6GF_R_KF6.AD'] + row['6GF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GF6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GF6_GF3KF6_R_GT.append('error')
        GF6_GF3KF6_R_AD.append('error')
        GF6_GF3KF6_A_GT.append('error')
        GF6_GF3KF6_A_AD.append('error')


# In[41]:


df_backup['6GF_R_GF3KF6.GT'] = GF6_GF3KF6_R_GT
df_backup['6GF_A_GF3KF6.GT'] = GF6_GF3KF6_A_GT
df_backup['6GF_R_GF3KF6.AD'] = GF6_GF3KF6_R_AD
df_backup['6GF_A_GF3KF6.AD'] = GF6_GF3KF6_A_AD


# GS1

# In[42]:


GS1_GF3KF6_R_GT = []
GS1_GF3KF6_R_AD = []
GS1_GF3KF6_A_GT = []
GS1_GF3KF6_A_AD = []


# In[43]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GS1_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['1GS_R_GF3.GT'] == row['1GS_A_GF3.GT']) and (row['1GS_R_KF6.GT'] == row['1GS_A_KF6.GT']) and (row['1GS_R_GF3.GT'] == row['1GS_R_KF6.GT']) and (row['1GS_A_GF3.GT'] == row['1GS_A_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = (row['1GS_R_GF3.AD'] + row['1GS_R_KF6.AD'] + row['1GS_A_GF3.AD'] + row['1GS_A_KF6.AD'])/2
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_R_KF6.GT'])
            y = 0
            GS1_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['1GS_R_GF3.GT'] == row['1GS_A_KF6.GT']) and (row['1GS_A_GF3.GT'] == row['1GS_R_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = (row['1GS_R_GF3.AD'] + row['1GS_A_KF6.AD'])/2
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_A_GF3.GT'])
            y = (row['1GS_A_GF3.AD'] + row['1GS_R_KF6.AD'])/2
            GS1_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['1GS_R_GF3.GT'] == row['1GS_R_KF6.GT']) and (row['1GS_A_GF3.GT'] == row['1GS_A_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = (row['1GS_R_GF3.AD'] + row['1GS_R_KF6.AD'])/2
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_A_GF3.GT'])
            y = (row['1GS_A_GF3.AD'] + row['1GS_A_KF6.AD'])/2
            GS1_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['1GS_R_GF3.GT'] == row['1GS_A_GF3.GT']) and (row['1GS_R_KF6.GT'] != row['1GS_A_KF6.GT']) and (row['1GS_R_GF3.GT'] == row['1GS_R_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = (row['1GS_R_GF3.AD'] + row['1GS_A_GF3.AD'] + row['1GS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_A_KF6.GT'])
            y = row['1GS_A_KF6.AD']
            GS1_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['1GS_R_GF3.GT'] == row['1GS_A_GF3.GT']) and (row['1GS_R_KF6.GT'] != row['1GS_A_KF6.GT']) and (row['1GS_R_GF3.GT'] == row['1GS_A_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = (row['1GS_R_GF3.AD'] + row['1GS_A_GF3.AD'] + row['1GS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_R_KF6.GT'])
            y = row['1GS_R_KF6.AD']
            GS1_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['1GS_R_GF3.GT'] != row['1GS_A_GF3.GT']) and (row['1GS_R_KF6.GT'] == row['1GS_A_KF6.GT']) and (row['1GS_R_GF3.GT'] == row['1GS_R_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = row['1GS_R_GF3.AD'] 
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_A_GF3.GT'])
            y = (row['1GS_A_GF3.AD'] + row['1GS_R_KF6.AD'] + row['1GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS1_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['1GS_R_GF3.GT'] != row['1GS_A_GF3.GT']) and (row['1GS_R_KF6.GT'] == row['1GS_A_KF6.GT']) and (row['1GS_R_GF3.GT'] == row['1GS_A_KF6.GT'])):
            GS1_GF3KF6_R_GT.append(row['1GS_R_GF3.GT'])
            x = row['1GS_R_GF3.AD']
            GS1_GF3KF6_R_AD.append(x)
            GS1_GF3KF6_A_GT.append(row['1GS_R_KF6.GT'])
            y = (row['1GS_A_GF3.AD'] + row['1GS_R_KF6.AD'] + row['1GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS1_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GS1_GF3KF6_R_GT.append('error')
        GS1_GF3KF6_R_AD.append('error')
        GS1_GF3KF6_A_GT.append('error')
        GS1_GF3KF6_A_AD.append('error')


# In[44]:


df_backup['1GS_R_GF3KF6.GT'] = GS1_GF3KF6_R_GT
df_backup['1GS_A_GF3KF6.GT'] = GS1_GF3KF6_A_GT
df_backup['1GS_R_GF3KF6.AD'] = GS1_GF3KF6_R_AD
df_backup['1GS_A_GF3KF6.AD'] = GS1_GF3KF6_A_AD


# GS2

# In[45]:


GS2_GF3KF6_R_GT = []
GS2_GF3KF6_R_AD = []
GS2_GF3KF6_A_GT = []
GS2_GF3KF6_A_AD = []


# In[46]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GS2_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['2GS_R_GF3.GT'] == row['2GS_A_GF3.GT']) and (row['2GS_R_KF6.GT'] == row['2GS_A_KF6.GT']) and (row['2GS_R_GF3.GT'] == row['2GS_R_KF6.GT']) and (row['2GS_A_GF3.GT'] == row['2GS_A_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = (row['2GS_R_GF3.AD'] + row['2GS_R_KF6.AD'] + row['2GS_A_GF3.AD'] + row['2GS_A_KF6.AD'])/2
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_R_KF6.GT'])
            y = 0
            GS2_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['2GS_R_GF3.GT'] == row['2GS_A_KF6.GT']) and (row['2GS_A_GF3.GT'] == row['2GS_R_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = (row['2GS_R_GF3.AD'] + row['2GS_A_KF6.AD'])/2
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_A_GF3.GT'])
            y = (row['2GS_A_GF3.AD'] + row['2GS_R_KF6.AD'])/2
            GS2_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['2GS_R_GF3.GT'] == row['2GS_R_KF6.GT']) and (row['2GS_A_GF3.GT'] == row['2GS_A_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = (row['2GS_R_GF3.AD'] + row['2GS_R_KF6.AD'])/2
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_A_GF3.GT'])
            y = (row['2GS_A_GF3.AD'] + row['2GS_A_KF6.AD'])/2
            GS2_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['2GS_R_GF3.GT'] == row['2GS_A_GF3.GT']) and (row['2GS_R_KF6.GT'] != row['2GS_A_KF6.GT']) and (row['2GS_R_GF3.GT'] == row['2GS_R_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = (row['2GS_R_GF3.AD'] + row['2GS_A_GF3.AD'] + row['2GS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_A_KF6.GT'])
            y = row['2GS_A_KF6.AD']
            GS2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['2GS_R_GF3.GT'] == row['2GS_A_GF3.GT']) and (row['2GS_R_KF6.GT'] != row['2GS_A_KF6.GT']) and (row['2GS_R_GF3.GT'] == row['2GS_A_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = (row['2GS_R_GF3.AD'] + row['2GS_A_GF3.AD'] + row['2GS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_R_KF6.GT'])
            y = row['2GS_R_KF6.AD']
            GS2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['2GS_R_GF3.GT'] != row['2GS_A_GF3.GT']) and (row['2GS_R_KF6.GT'] == row['2GS_A_KF6.GT']) and (row['2GS_R_GF3.GT'] == row['2GS_R_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = row['2GS_R_GF3.AD'] 
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_A_GF3.GT'])
            y = (row['2GS_A_GF3.AD'] + row['2GS_R_KF6.AD'] + row['2GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['2GS_R_GF3.GT'] != row['2GS_A_GF3.GT']) and (row['2GS_R_KF6.GT'] == row['2GS_A_KF6.GT']) and (row['2GS_R_GF3.GT'] == row['2GS_A_KF6.GT'])):
            GS2_GF3KF6_R_GT.append(row['2GS_R_GF3.GT'])
            x = row['2GS_R_GF3.AD']
            GS2_GF3KF6_R_AD.append(x)
            GS2_GF3KF6_A_GT.append(row['2GS_R_KF6.GT'])
            y = (row['2GS_A_GF3.AD'] + row['2GS_R_KF6.AD'] + row['2GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GS2_GF3KF6_R_GT.append('error')
        GS2_GF3KF6_R_AD.append('error')
        GS2_GF3KF6_A_GT.append('error')
        GS2_GF3KF6_A_AD.append('error')


# In[47]:


df_backup['2GS_R_GF3KF6.GT'] = GS2_GF3KF6_R_GT
df_backup['2GS_A_GF3KF6.GT'] = GS2_GF3KF6_A_GT
df_backup['2GS_R_GF3KF6.AD'] = GS2_GF3KF6_R_AD
df_backup['2GS_A_GF3KF6.AD'] = GS2_GF3KF6_A_AD


# GS3

# In[48]:


GS3_GF3KF6_R_GT = []
GS3_GF3KF6_R_AD = []
GS3_GF3KF6_A_GT = []
GS3_GF3KF6_A_AD = []


# In[49]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GS3_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['3GS_R_GF3.GT'] == row['3GS_A_GF3.GT']) and (row['3GS_R_KF6.GT'] == row['3GS_A_KF6.GT']) and (row['3GS_R_GF3.GT'] == row['3GS_R_KF6.GT']) and (row['3GS_A_GF3.GT'] == row['3GS_A_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = (row['3GS_R_GF3.AD'] + row['3GS_R_KF6.AD'] + row['3GS_A_GF3.AD'] + row['3GS_A_KF6.AD'])/2
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_R_KF6.GT'])
            y = 0
            GS3_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['3GS_R_GF3.GT'] == row['3GS_A_KF6.GT']) and (row['3GS_A_GF3.GT'] == row['3GS_R_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = (row['3GS_R_GF3.AD'] + row['3GS_A_KF6.AD'])/2
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_A_GF3.GT'])
            y = (row['3GS_A_GF3.AD'] + row['3GS_R_KF6.AD'])/2
            GS3_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['3GS_R_GF3.GT'] == row['3GS_R_KF6.GT']) and (row['3GS_A_GF3.GT'] == row['3GS_A_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = (row['3GS_R_GF3.AD'] + row['3GS_R_KF6.AD'])/2
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_A_GF3.GT'])
            y = (row['3GS_A_GF3.AD'] + row['3GS_A_KF6.AD'])/2
            GS3_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['3GS_R_GF3.GT'] == row['3GS_A_GF3.GT']) and (row['3GS_R_KF6.GT'] != row['3GS_A_KF6.GT']) and (row['3GS_R_GF3.GT'] == row['3GS_R_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = (row['3GS_R_GF3.AD'] + row['3GS_A_GF3.AD'] + row['3GS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_A_KF6.GT'])
            y = row['3GS_A_KF6.AD']
            GS3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['3GS_R_GF3.GT'] == row['3GS_A_GF3.GT']) and (row['3GS_R_KF6.GT'] != row['3GS_A_KF6.GT']) and (row['3GS_R_GF3.GT'] == row['3GS_A_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = (row['3GS_R_GF3.AD'] + row['3GS_A_GF3.AD'] + row['3GS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_R_KF6.GT'])
            y = row['3GS_R_KF6.AD']
            GS3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['3GS_R_GF3.GT'] != row['3GS_A_GF3.GT']) and (row['3GS_R_KF6.GT'] == row['3GS_A_KF6.GT']) and (row['3GS_R_GF3.GT'] == row['3GS_R_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = row['3GS_R_GF3.AD'] 
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_A_GF3.GT'])
            y = (row['3GS_A_GF3.AD'] + row['3GS_R_KF6.AD'] + row['3GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['3GS_R_GF3.GT'] != row['3GS_A_GF3.GT']) and (row['3GS_R_KF6.GT'] == row['3GS_A_KF6.GT']) and (row['3GS_R_GF3.GT'] == row['3GS_A_KF6.GT'])):
            GS3_GF3KF6_R_GT.append(row['3GS_R_GF3.GT'])
            x = row['3GS_R_GF3.AD']
            GS3_GF3KF6_R_AD.append(x)
            GS3_GF3KF6_A_GT.append(row['3GS_R_KF6.GT'])
            y = (row['3GS_A_GF3.AD'] + row['3GS_R_KF6.AD'] + row['3GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GS3_GF3KF6_R_GT.append('error')
        GS3_GF3KF6_R_AD.append('error')
        GS3_GF3KF6_A_GT.append('error')
        GS3_GF3KF6_A_AD.append('error')


# In[50]:


df_backup['3GS_R_GF3KF6.GT'] = GS3_GF3KF6_R_GT
df_backup['3GS_A_GF3KF6.GT'] = GS3_GF3KF6_A_GT
df_backup['3GS_R_GF3KF6.AD'] = GS3_GF3KF6_R_AD
df_backup['3GS_A_GF3KF6.AD'] = GS3_GF3KF6_A_AD


# GS4

# In[51]:


GS4_GF3KF6_R_GT = []
GS4_GF3KF6_R_AD = []
GS4_GF3KF6_A_GT = []
GS4_GF3KF6_A_AD = []


# In[52]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GS4_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['4GS_R_GF3.GT'] == row['4GS_A_GF3.GT']) and (row['4GS_R_KF6.GT'] == row['4GS_A_KF6.GT']) and (row['4GS_R_GF3.GT'] == row['4GS_R_KF6.GT']) and (row['4GS_A_GF3.GT'] == row['4GS_A_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = (row['4GS_R_GF3.AD'] + row['4GS_R_KF6.AD'] + row['4GS_A_GF3.AD'] + row['4GS_A_KF6.AD'])/2
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_R_KF6.GT'])
            y = 0
            GS4_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['4GS_R_GF3.GT'] == row['4GS_A_KF6.GT']) and (row['4GS_A_GF3.GT'] == row['4GS_R_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = (row['4GS_R_GF3.AD'] + row['4GS_A_KF6.AD'])/2
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_A_GF3.GT'])
            y = (row['4GS_A_GF3.AD'] + row['4GS_R_KF6.AD'])/2
            GS4_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['4GS_R_GF3.GT'] == row['4GS_R_KF6.GT']) and (row['4GS_A_GF3.GT'] == row['4GS_A_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = (row['4GS_R_GF3.AD'] + row['4GS_R_KF6.AD'])/2
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_A_GF3.GT'])
            y = (row['4GS_A_GF3.AD'] + row['4GS_A_KF6.AD'])/2
            GS4_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['4GS_R_GF3.GT'] == row['4GS_A_GF3.GT']) and (row['4GS_R_KF6.GT'] != row['4GS_A_KF6.GT']) and (row['4GS_R_GF3.GT'] == row['4GS_R_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = (row['4GS_R_GF3.AD'] + row['4GS_A_GF3.AD'] + row['4GS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_A_KF6.GT'])
            y = row['4GS_A_KF6.AD']
            GS4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['4GS_R_GF3.GT'] == row['4GS_A_GF3.GT']) and (row['4GS_R_KF6.GT'] != row['4GS_A_KF6.GT']) and (row['4GS_R_GF3.GT'] == row['4GS_A_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = (row['4GS_R_GF3.AD'] + row['4GS_A_GF3.AD'] + row['4GS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_R_KF6.GT'])
            y = row['4GS_R_KF6.AD']
            GS4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['4GS_R_GF3.GT'] != row['4GS_A_GF3.GT']) and (row['4GS_R_KF6.GT'] == row['4GS_A_KF6.GT']) and (row['4GS_R_GF3.GT'] == row['4GS_R_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = row['4GS_R_GF3.AD'] 
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_A_GF3.GT'])
            y = (row['4GS_A_GF3.AD'] + row['4GS_R_KF6.AD'] + row['4GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['4GS_R_GF3.GT'] != row['4GS_A_GF3.GT']) and (row['4GS_R_KF6.GT'] == row['4GS_A_KF6.GT']) and (row['4GS_R_GF3.GT'] == row['4GS_A_KF6.GT'])):
            GS4_GF3KF6_R_GT.append(row['4GS_R_GF3.GT'])
            x = row['4GS_R_GF3.AD']
            GS4_GF3KF6_R_AD.append(x)
            GS4_GF3KF6_A_GT.append(row['4GS_R_KF6.GT'])
            y = (row['4GS_A_GF3.AD'] + row['4GS_R_KF6.AD'] + row['4GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GS4_GF3KF6_R_GT.append('error')
        GS4_GF3KF6_R_AD.append('error')
        GS4_GF3KF6_A_GT.append('error')
        GS4_GF3KF6_A_AD.append('error')


# In[53]:


df_backup['4GS_R_GF3KF6.GT'] = GS4_GF3KF6_R_GT
df_backup['4GS_A_GF3KF6.GT'] = GS4_GF3KF6_A_GT
df_backup['4GS_R_GF3KF6.AD'] = GS4_GF3KF6_R_AD
df_backup['4GS_A_GF3KF6.AD'] = GS4_GF3KF6_A_AD


# GS5

# In[54]:


GS5_GF3KF6_R_GT = []
GS5_GF3KF6_R_AD = []
GS5_GF3KF6_A_GT = []
GS5_GF3KF6_A_AD = []


# In[55]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GS5_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['5GS_R_GF3.GT'] == row['5GS_A_GF3.GT']) and (row['5GS_R_KF6.GT'] == row['5GS_A_KF6.GT']) and (row['5GS_R_GF3.GT'] == row['5GS_R_KF6.GT']) and (row['5GS_A_GF3.GT'] == row['5GS_A_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = (row['5GS_R_GF3.AD'] + row['5GS_R_KF6.AD'] + row['5GS_A_GF3.AD'] + row['5GS_A_KF6.AD'])/2
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_R_KF6.GT'])
            y = 0
            GS5_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['5GS_R_GF3.GT'] == row['5GS_A_KF6.GT']) and (row['5GS_A_GF3.GT'] == row['5GS_R_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = (row['5GS_R_GF3.AD'] + row['5GS_A_KF6.AD'])/2
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_A_GF3.GT'])
            y = (row['5GS_A_GF3.AD'] + row['5GS_R_KF6.AD'])/2
            GS5_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['5GS_R_GF3.GT'] == row['5GS_R_KF6.GT']) and (row['5GS_A_GF3.GT'] == row['5GS_A_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = (row['5GS_R_GF3.AD'] + row['5GS_R_KF6.AD'])/2
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_A_GF3.GT'])
            y = (row['5GS_A_GF3.AD'] + row['5GS_A_KF6.AD'])/2
            GS5_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['5GS_R_GF3.GT'] == row['5GS_A_GF3.GT']) and (row['5GS_R_KF6.GT'] != row['5GS_A_KF6.GT']) and (row['5GS_R_GF3.GT'] == row['5GS_R_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = (row['5GS_R_GF3.AD'] + row['5GS_A_GF3.AD'] + row['5GS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_A_KF6.GT'])
            y = row['5GS_A_KF6.AD']
            GS5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['5GS_R_GF3.GT'] == row['5GS_A_GF3.GT']) and (row['5GS_R_KF6.GT'] != row['5GS_A_KF6.GT']) and (row['5GS_R_GF3.GT'] == row['5GS_A_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = (row['5GS_R_GF3.AD'] + row['5GS_A_GF3.AD'] + row['5GS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_R_KF6.GT'])
            y = row['5GS_R_KF6.AD']
            GS5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['5GS_R_GF3.GT'] != row['5GS_A_GF3.GT']) and (row['5GS_R_KF6.GT'] == row['5GS_A_KF6.GT']) and (row['5GS_R_GF3.GT'] == row['5GS_R_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = row['5GS_R_GF3.AD'] 
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_A_GF3.GT'])
            y = (row['5GS_A_GF3.AD'] + row['5GS_R_KF6.AD'] + row['5GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['5GS_R_GF3.GT'] != row['5GS_A_GF3.GT']) and (row['5GS_R_KF6.GT'] == row['5GS_A_KF6.GT']) and (row['5GS_R_GF3.GT'] == row['5GS_A_KF6.GT'])):
            GS5_GF3KF6_R_GT.append(row['5GS_R_GF3.GT'])
            x = row['5GS_R_GF3.AD']
            GS5_GF3KF6_R_AD.append(x)
            GS5_GF3KF6_A_GT.append(row['5GS_R_KF6.GT'])
            y = (row['5GS_A_GF3.AD'] + row['5GS_R_KF6.AD'] + row['5GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GS5_GF3KF6_R_GT.append('error')
        GS5_GF3KF6_R_AD.append('error')
        GS5_GF3KF6_A_GT.append('error')
        GS5_GF3KF6_A_AD.append('error')


# In[56]:


df_backup['5GS_R_GF3KF6.GT'] = GS5_GF3KF6_R_GT
df_backup['5GS_A_GF3KF6.GT'] = GS5_GF3KF6_A_GT
df_backup['5GS_R_GF3KF6.AD'] = GS5_GF3KF6_R_AD
df_backup['5GS_A_GF3KF6.AD'] = GS5_GF3KF6_A_AD


# GS6

# In[57]:


GS6_GF3KF6_R_GT = []
GS6_GF3KF6_R_AD = []
GS6_GF3KF6_A_GT = []
GS6_GF3KF6_A_AD = []


# In[58]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(GS6_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['6GS_R_GF3.GT'] == row['6GS_A_GF3.GT']) and (row['6GS_R_KF6.GT'] == row['6GS_A_KF6.GT']) and (row['6GS_R_GF3.GT'] == row['6GS_R_KF6.GT']) and (row['6GS_A_GF3.GT'] == row['6GS_A_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = (row['6GS_R_GF3.AD'] + row['6GS_R_KF6.AD'] + row['6GS_A_GF3.AD'] + row['6GS_A_KF6.AD'])/2
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_R_KF6.GT'])
            y = 0
            GS6_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['6GS_R_GF3.GT'] == row['6GS_A_KF6.GT']) and (row['6GS_A_GF3.GT'] == row['6GS_R_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = (row['6GS_R_GF3.AD'] + row['6GS_A_KF6.AD'])/2
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_A_GF3.GT'])
            y = (row['6GS_A_GF3.AD'] + row['6GS_R_KF6.AD'])/2
            GS6_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['6GS_R_GF3.GT'] == row['6GS_R_KF6.GT']) and (row['6GS_A_GF3.GT'] == row['6GS_A_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = (row['6GS_R_GF3.AD'] + row['6GS_R_KF6.AD'])/2
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_A_GF3.GT'])
            y = (row['6GS_A_GF3.AD'] + row['6GS_A_KF6.AD'])/2
            GS6_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['6GS_R_GF3.GT'] == row['6GS_A_GF3.GT']) and (row['6GS_R_KF6.GT'] != row['6GS_A_KF6.GT']) and (row['6GS_R_GF3.GT'] == row['6GS_R_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = (row['6GS_R_GF3.AD'] + row['6GS_A_GF3.AD'] + row['6GS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_A_KF6.GT'])
            y = row['6GS_A_KF6.AD']
            GS6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['6GS_R_GF3.GT'] == row['6GS_A_GF3.GT']) and (row['6GS_R_KF6.GT'] != row['6GS_A_KF6.GT']) and (row['6GS_R_GF3.GT'] == row['6GS_A_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = (row['6GS_R_GF3.AD'] + row['6GS_A_GF3.AD'] + row['6GS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_R_KF6.GT'])
            y = row['6GS_R_KF6.AD']
            GS6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['6GS_R_GF3.GT'] != row['6GS_A_GF3.GT']) and (row['6GS_R_KF6.GT'] == row['6GS_A_KF6.GT']) and (row['6GS_R_GF3.GT'] == row['6GS_R_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = row['6GS_R_GF3.AD'] 
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_A_GF3.GT'])
            y = (row['6GS_A_GF3.AD'] + row['6GS_R_KF6.AD'] + row['6GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['6GS_R_GF3.GT'] != row['6GS_A_GF3.GT']) and (row['6GS_R_KF6.GT'] == row['6GS_A_KF6.GT']) and (row['6GS_R_GF3.GT'] == row['6GS_A_KF6.GT'])):
            GS6_GF3KF6_R_GT.append(row['6GS_R_GF3.GT'])
            x = row['6GS_R_GF3.AD']
            GS6_GF3KF6_R_AD.append(x)
            GS6_GF3KF6_A_GT.append(row['6GS_R_KF6.GT'])
            y = (row['6GS_A_GF3.AD'] + row['6GS_R_KF6.AD'] + row['6GS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            GS6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        GS6_GF3KF6_R_GT.append('error')
        GS6_GF3KF6_R_AD.append('error')
        GS6_GF3KF6_A_GT.append('error')
        GS6_GF3KF6_A_AD.append('error')


# In[59]:


df_backup['6GS_R_GF3KF6.GT'] = GS6_GF3KF6_R_GT
df_backup['6GS_A_GF3KF6.GT'] = GS6_GF3KF6_A_GT
df_backup['6GS_R_GF3KF6.AD'] = GS6_GF3KF6_R_AD
df_backup['6GS_A_GF3KF6.AD'] = GS6_GF3KF6_A_AD


# KF1

# In[60]:


KF1_GF3KF6_R_GT = []
KF1_GF3KF6_R_AD = []
KF1_GF3KF6_A_GT = []
KF1_GF3KF6_A_AD = []


# In[61]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KF1_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['1KF_R_GF3.GT'] == row['1KF_A_GF3.GT']) and (row['1KF_R_KF6.GT'] == row['1KF_A_KF6.GT']) and (row['1KF_R_GF3.GT'] == row['1KF_R_KF6.GT']) and (row['1KF_A_GF3.GT'] == row['1KF_A_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = (row['1KF_R_GF3.AD'] + row['1KF_R_KF6.AD'] + row['1KF_A_GF3.AD'] + row['1KF_A_KF6.AD'])/2
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_R_KF6.GT'])
            y = 0
            KF1_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['1KF_R_GF3.GT'] == row['1KF_A_KF6.GT']) and (row['1KF_A_GF3.GT'] == row['1KF_R_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = (row['1KF_R_GF3.AD'] + row['1KF_A_KF6.AD'])/2
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_A_GF3.GT'])
            y = (row['1KF_A_GF3.AD'] + row['1KF_R_KF6.AD'])/2
            KF1_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['1KF_R_GF3.GT'] == row['1KF_R_KF6.GT']) and (row['1KF_A_GF3.GT'] == row['1KF_A_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = (row['1KF_R_GF3.AD'] + row['1KF_R_KF6.AD'])/2
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_A_GF3.GT'])
            y = (row['1KF_A_GF3.AD'] + row['1KF_A_KF6.AD'])/2
            KF1_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['1KF_R_GF3.GT'] == row['1KF_A_GF3.GT']) and (row['1KF_R_KF6.GT'] != row['1KF_A_KF6.GT']) and (row['1KF_R_GF3.GT'] == row['1KF_R_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = (row['1KF_R_GF3.AD'] + row['1KF_A_GF3.AD'] + row['1KF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_A_KF6.GT'])
            y = row['1KF_A_KF6.AD']
            KF1_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['1KF_R_GF3.GT'] == row['1KF_A_GF3.GT']) and (row['1KF_R_KF6.GT'] != row['1KF_A_KF6.GT']) and (row['1KF_R_GF3.GT'] == row['1KF_A_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = (row['1KF_R_GF3.AD'] + row['1KF_A_GF3.AD'] + row['1KF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_R_KF6.GT'])
            y = row['1KF_R_KF6.AD']
            KF1_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['1KF_R_GF3.GT'] != row['1KF_A_GF3.GT']) and (row['1KF_R_KF6.GT'] == row['1KF_A_KF6.GT']) and (row['1KF_R_GF3.GT'] == row['1KF_R_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = row['1KF_R_GF3.AD'] 
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_A_GF3.GT'])
            y = (row['1KF_A_GF3.AD'] + row['1KF_R_KF6.AD'] + row['1KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF1_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['1KF_R_GF3.GT'] != row['1KF_A_GF3.GT']) and (row['1KF_R_KF6.GT'] == row['1KF_A_KF6.GT']) and (row['1KF_R_GF3.GT'] == row['1KF_A_KF6.GT'])):
            KF1_GF3KF6_R_GT.append(row['1KF_R_GF3.GT'])
            x = row['1KF_R_GF3.AD']
            KF1_GF3KF6_R_AD.append(x)
            KF1_GF3KF6_A_GT.append(row['1KF_R_KF6.GT'])
            y = (row['1KF_A_GF3.AD'] + row['1KF_R_KF6.AD'] + row['1KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF1_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KF1_GF3KF6_R_GT.append('error')
        KF1_GF3KF6_R_AD.append('error')
        KF1_GF3KF6_A_GT.append('error')
        KF1_GF3KF6_A_AD.append('error')


# In[62]:


df_backup['1KF_R_GF3KF6.GT'] = KF1_GF3KF6_R_GT
df_backup['1KF_A_GF3KF6.GT'] = KF1_GF3KF6_A_GT
df_backup['1KF_R_GF3KF6.AD'] = KF1_GF3KF6_R_AD
df_backup['1KF_A_GF3KF6.AD'] = KF1_GF3KF6_A_AD


# KF2

# In[63]:


KF2_GF3KF6_R_GT = []
KF2_GF3KF6_R_AD = []
KF2_GF3KF6_A_GT = []
KF2_GF3KF6_A_AD = []


# In[64]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KF2_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['2KF_R_GF3.GT'] == row['2KF_A_GF3.GT']) and (row['2KF_R_KF6.GT'] == row['2KF_A_KF6.GT']) and (row['2KF_R_GF3.GT'] == row['2KF_R_KF6.GT']) and (row['2KF_A_GF3.GT'] == row['2KF_A_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = (row['2KF_R_GF3.AD'] + row['2KF_R_KF6.AD'] + row['2KF_A_GF3.AD'] + row['2KF_A_KF6.AD'])/2
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_R_KF6.GT'])
            y = 0
            KF2_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['2KF_R_GF3.GT'] == row['2KF_A_KF6.GT']) and (row['2KF_A_GF3.GT'] == row['2KF_R_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = (row['2KF_R_GF3.AD'] + row['2KF_A_KF6.AD'])/2
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_A_GF3.GT'])
            y = (row['2KF_A_GF3.AD'] + row['2KF_R_KF6.AD'])/2
            KF2_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['2KF_R_GF3.GT'] == row['2KF_R_KF6.GT']) and (row['2KF_A_GF3.GT'] == row['2KF_A_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = (row['2KF_R_GF3.AD'] + row['2KF_R_KF6.AD'])/2
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_A_GF3.GT'])
            y = (row['2KF_A_GF3.AD'] + row['2KF_A_KF6.AD'])/2
            KF2_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['2KF_R_GF3.GT'] == row['2KF_A_GF3.GT']) and (row['2KF_R_KF6.GT'] != row['2KF_A_KF6.GT']) and (row['2KF_R_GF3.GT'] == row['2KF_R_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = (row['2KF_R_GF3.AD'] + row['2KF_A_GF3.AD'] + row['2KF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_A_KF6.GT'])
            y = row['2KF_A_KF6.AD']
            KF2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['2KF_R_GF3.GT'] == row['2KF_A_GF3.GT']) and (row['2KF_R_KF6.GT'] != row['2KF_A_KF6.GT']) and (row['2KF_R_GF3.GT'] == row['2KF_A_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = (row['2KF_R_GF3.AD'] + row['2KF_A_GF3.AD'] + row['2KF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_R_KF6.GT'])
            y = row['2KF_R_KF6.AD']
            KF2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['2KF_R_GF3.GT'] != row['2KF_A_GF3.GT']) and (row['2KF_R_KF6.GT'] == row['2KF_A_KF6.GT']) and (row['2KF_R_GF3.GT'] == row['2KF_R_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = row['2KF_R_GF3.AD'] 
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_A_GF3.GT'])
            y = (row['2KF_A_GF3.AD'] + row['2KF_R_KF6.AD'] + row['2KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['2KF_R_GF3.GT'] != row['2KF_A_GF3.GT']) and (row['2KF_R_KF6.GT'] == row['2KF_A_KF6.GT']) and (row['2KF_R_GF3.GT'] == row['2KF_A_KF6.GT'])):
            KF2_GF3KF6_R_GT.append(row['2KF_R_GF3.GT'])
            x = row['2KF_R_GF3.AD']
            KF2_GF3KF6_R_AD.append(x)
            KF2_GF3KF6_A_GT.append(row['2KF_R_KF6.GT'])
            y = (row['2KF_A_GF3.AD'] + row['2KF_R_KF6.AD'] + row['2KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KF2_GF3KF6_R_GT.append('error')
        KF2_GF3KF6_R_AD.append('error')
        KF2_GF3KF6_A_GT.append('error')
        KF2_GF3KF6_A_AD.append('error')


# In[65]:


df_backup['2KF_R_GF3KF6.GT'] = KF2_GF3KF6_R_GT
df_backup['2KF_A_GF3KF6.GT'] = KF2_GF3KF6_A_GT
df_backup['2KF_R_GF3KF6.AD'] = KF2_GF3KF6_R_AD
df_backup['2KF_A_GF3KF6.AD'] = KF2_GF3KF6_A_AD


# KF3

# In[66]:


KF3_GF3KF6_R_GT = []
KF3_GF3KF6_R_AD = []
KF3_GF3KF6_A_GT = []
KF3_GF3KF6_A_AD = []


# In[67]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KF3_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['3KF_R_GF3.GT'] == row['3KF_A_GF3.GT']) and (row['3KF_R_KF6.GT'] == row['3KF_A_KF6.GT']) and (row['3KF_R_GF3.GT'] == row['3KF_R_KF6.GT']) and (row['3KF_A_GF3.GT'] == row['3KF_A_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = (row['3KF_R_GF3.AD'] + row['3KF_R_KF6.AD'] + row['3KF_A_GF3.AD'] + row['3KF_A_KF6.AD'])/2
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_R_KF6.GT'])
            y = 0
            KF3_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['3KF_R_GF3.GT'] == row['3KF_A_KF6.GT']) and (row['3KF_A_GF3.GT'] == row['3KF_R_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = (row['3KF_R_GF3.AD'] + row['3KF_A_KF6.AD'])/2
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_A_GF3.GT'])
            y = (row['3KF_A_GF3.AD'] + row['3KF_R_KF6.AD'])/2
            KF3_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['3KF_R_GF3.GT'] == row['3KF_R_KF6.GT']) and (row['3KF_A_GF3.GT'] == row['3KF_A_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = (row['3KF_R_GF3.AD'] + row['3KF_R_KF6.AD'])/2
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_A_GF3.GT'])
            y = (row['3KF_A_GF3.AD'] + row['3KF_A_KF6.AD'])/2
            KF3_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['3KF_R_GF3.GT'] == row['3KF_A_GF3.GT']) and (row['3KF_R_KF6.GT'] != row['3KF_A_KF6.GT']) and (row['3KF_R_GF3.GT'] == row['3KF_R_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = (row['3KF_R_GF3.AD'] + row['3KF_A_GF3.AD'] + row['3KF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_A_KF6.GT'])
            y = row['3KF_A_KF6.AD']
            KF3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['3KF_R_GF3.GT'] == row['3KF_A_GF3.GT']) and (row['3KF_R_KF6.GT'] != row['3KF_A_KF6.GT']) and (row['3KF_R_GF3.GT'] == row['3KF_A_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = (row['3KF_R_GF3.AD'] + row['3KF_A_GF3.AD'] + row['3KF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_R_KF6.GT'])
            y = row['3KF_R_KF6.AD']
            KF3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['3KF_R_GF3.GT'] != row['3KF_A_GF3.GT']) and (row['3KF_R_KF6.GT'] == row['3KF_A_KF6.GT']) and (row['3KF_R_GF3.GT'] == row['3KF_R_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = row['3KF_R_GF3.AD'] 
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_A_GF3.GT'])
            y = (row['3KF_A_GF3.AD'] + row['3KF_R_KF6.AD'] + row['3KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['3KF_R_GF3.GT'] != row['3KF_A_GF3.GT']) and (row['3KF_R_KF6.GT'] == row['3KF_A_KF6.GT']) and (row['3KF_R_GF3.GT'] == row['3KF_A_KF6.GT'])):
            KF3_GF3KF6_R_GT.append(row['3KF_R_GF3.GT'])
            x = row['3KF_R_GF3.AD']
            KF3_GF3KF6_R_AD.append(x)
            KF3_GF3KF6_A_GT.append(row['3KF_R_KF6.GT'])
            y = (row['3KF_A_GF3.AD'] + row['3KF_R_KF6.AD'] + row['3KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KF3_GF3KF6_R_GT.append('error')
        KF3_GF3KF6_R_AD.append('error')
        KF3_GF3KF6_A_GT.append('error')
        KF3_GF3KF6_A_AD.append('error')


# In[68]:


df_backup['3KF_R_GF3KF6.GT'] = KF3_GF3KF6_R_GT
df_backup['3KF_A_GF3KF6.GT'] = KF3_GF3KF6_A_GT
df_backup['3KF_R_GF3KF6.AD'] = KF3_GF3KF6_R_AD
df_backup['3KF_A_GF3KF6.AD'] = KF3_GF3KF6_A_AD


# KF4

# In[69]:


KF4_GF3KF6_R_GT = []
KF4_GF3KF6_R_AD = []
KF4_GF3KF6_A_GT = []
KF4_GF3KF6_A_AD = []


# In[70]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KF4_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['4KF_R_GF3.GT'] == row['4KF_A_GF3.GT']) and (row['4KF_R_KF6.GT'] == row['4KF_A_KF6.GT']) and (row['4KF_R_GF3.GT'] == row['4KF_R_KF6.GT']) and (row['4KF_A_GF3.GT'] == row['4KF_A_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = (row['4KF_R_GF3.AD'] + row['4KF_R_KF6.AD'] + row['4KF_A_GF3.AD'] + row['4KF_A_KF6.AD'])/2
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_R_KF6.GT'])
            y = 0
            KF4_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['4KF_R_GF3.GT'] == row['4KF_A_KF6.GT']) and (row['4KF_A_GF3.GT'] == row['4KF_R_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = (row['4KF_R_GF3.AD'] + row['4KF_A_KF6.AD'])/2
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_A_GF3.GT'])
            y = (row['4KF_A_GF3.AD'] + row['4KF_R_KF6.AD'])/2
            KF4_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['4KF_R_GF3.GT'] == row['4KF_R_KF6.GT']) and (row['4KF_A_GF3.GT'] == row['4KF_A_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = (row['4KF_R_GF3.AD'] + row['4KF_R_KF6.AD'])/2
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_A_GF3.GT'])
            y = (row['4KF_A_GF3.AD'] + row['4KF_A_KF6.AD'])/2
            KF4_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['4KF_R_GF3.GT'] == row['4KF_A_GF3.GT']) and (row['4KF_R_KF6.GT'] != row['4KF_A_KF6.GT']) and (row['4KF_R_GF3.GT'] == row['4KF_R_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = (row['4KF_R_GF3.AD'] + row['4KF_A_GF3.AD'] + row['4KF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_A_KF6.GT'])
            y = row['4KF_A_KF6.AD']
            KF4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['4KF_R_GF3.GT'] == row['4KF_A_GF3.GT']) and (row['4KF_R_KF6.GT'] != row['4KF_A_KF6.GT']) and (row['4KF_R_GF3.GT'] == row['4KF_A_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = (row['4KF_R_GF3.AD'] + row['4KF_A_GF3.AD'] + row['4KF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_R_KF6.GT'])
            y = row['4KF_R_KF6.AD']
            KF4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['4KF_R_GF3.GT'] != row['4KF_A_GF3.GT']) and (row['4KF_R_KF6.GT'] == row['4KF_A_KF6.GT']) and (row['4KF_R_GF3.GT'] == row['4KF_R_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = row['4KF_R_GF3.AD'] 
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_A_GF3.GT'])
            y = (row['4KF_A_GF3.AD'] + row['4KF_R_KF6.AD'] + row['4KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['4KF_R_GF3.GT'] != row['4KF_A_GF3.GT']) and (row['4KF_R_KF6.GT'] == row['4KF_A_KF6.GT']) and (row['4KF_R_GF3.GT'] == row['4KF_A_KF6.GT'])):
            KF4_GF3KF6_R_GT.append(row['4KF_R_GF3.GT'])
            x = row['4KF_R_GF3.AD']
            KF4_GF3KF6_R_AD.append(x)
            KF4_GF3KF6_A_GT.append(row['4KF_R_KF6.GT'])
            y = (row['4KF_A_GF3.AD'] + row['4KF_R_KF6.AD'] + row['4KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KF4_GF3KF6_R_GT.append('error')
        KF4_GF3KF6_R_AD.append('error')
        KF4_GF3KF6_A_GT.append('error')
        KF4_GF3KF6_A_AD.append('error')


# In[71]:


df_backup['4KF_R_GF3KF6.GT'] = KF4_GF3KF6_R_GT
df_backup['4KF_A_GF3KF6.GT'] = KF4_GF3KF6_A_GT
df_backup['4KF_R_GF3KF6.AD'] = KF4_GF3KF6_R_AD
df_backup['4KF_A_GF3KF6.AD'] = KF4_GF3KF6_A_AD


# KF5

# In[72]:


KF5_GF3KF6_R_GT = []
KF5_GF3KF6_R_AD = []
KF5_GF3KF6_A_GT = []
KF5_GF3KF6_A_AD = []


# In[73]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KF5_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['5KF_R_GF3.GT'] == row['5KF_A_GF3.GT']) and (row['5KF_R_KF6.GT'] == row['5KF_A_KF6.GT']) and (row['5KF_R_GF3.GT'] == row['5KF_R_KF6.GT']) and (row['5KF_A_GF3.GT'] == row['5KF_A_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = (row['5KF_R_GF3.AD'] + row['5KF_R_KF6.AD'] + row['5KF_A_GF3.AD'] + row['5KF_A_KF6.AD'])/2
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_R_KF6.GT'])
            y = 0
            KF5_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['5KF_R_GF3.GT'] == row['5KF_A_KF6.GT']) and (row['5KF_A_GF3.GT'] == row['5KF_R_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = (row['5KF_R_GF3.AD'] + row['5KF_A_KF6.AD'])/2
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_A_GF3.GT'])
            y = (row['5KF_A_GF3.AD'] + row['5KF_R_KF6.AD'])/2
            KF5_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['5KF_R_GF3.GT'] == row['5KF_R_KF6.GT']) and (row['5KF_A_GF3.GT'] == row['5KF_A_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = (row['5KF_R_GF3.AD'] + row['5KF_R_KF6.AD'])/2
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_A_GF3.GT'])
            y = (row['5KF_A_GF3.AD'] + row['5KF_A_KF6.AD'])/2
            KF5_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['5KF_R_GF3.GT'] == row['5KF_A_GF3.GT']) and (row['5KF_R_KF6.GT'] != row['5KF_A_KF6.GT']) and (row['5KF_R_GF3.GT'] == row['5KF_R_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = (row['5KF_R_GF3.AD'] + row['5KF_A_GF3.AD'] + row['5KF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_A_KF6.GT'])
            y = row['5KF_A_KF6.AD']
            KF5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['5KF_R_GF3.GT'] == row['5KF_A_GF3.GT']) and (row['5KF_R_KF6.GT'] != row['5KF_A_KF6.GT']) and (row['5KF_R_GF3.GT'] == row['5KF_A_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = (row['5KF_R_GF3.AD'] + row['5KF_A_GF3.AD'] + row['5KF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_R_KF6.GT'])
            y = row['5KF_R_KF6.AD']
            KF5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['5KF_R_GF3.GT'] != row['5KF_A_GF3.GT']) and (row['5KF_R_KF6.GT'] == row['5KF_A_KF6.GT']) and (row['5KF_R_GF3.GT'] == row['5KF_R_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = row['5KF_R_GF3.AD'] 
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_A_GF3.GT'])
            y = (row['5KF_A_GF3.AD'] + row['5KF_R_KF6.AD'] + row['5KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['5KF_R_GF3.GT'] != row['5KF_A_GF3.GT']) and (row['5KF_R_KF6.GT'] == row['5KF_A_KF6.GT']) and (row['5KF_R_GF3.GT'] == row['5KF_A_KF6.GT'])):
            KF5_GF3KF6_R_GT.append(row['5KF_R_GF3.GT'])
            x = row['5KF_R_GF3.AD']
            KF5_GF3KF6_R_AD.append(x)
            KF5_GF3KF6_A_GT.append(row['5KF_R_KF6.GT'])
            y = (row['5KF_A_GF3.AD'] + row['5KF_R_KF6.AD'] + row['5KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KF5_GF3KF6_R_GT.append('error')
        KF5_GF3KF6_R_AD.append('error')
        KF5_GF3KF6_A_GT.append('error')
        KF5_GF3KF6_A_AD.append('error')


# In[74]:


df_backup['5KF_R_GF3KF6.GT'] = KF5_GF3KF6_R_GT
df_backup['5KF_A_GF3KF6.GT'] = KF5_GF3KF6_A_GT
df_backup['5KF_R_GF3KF6.AD'] = KF5_GF3KF6_R_AD
df_backup['5KF_A_GF3KF6.AD'] = KF5_GF3KF6_A_AD


# KF6

# In[75]:


KF6_GF3KF6_R_GT = []
KF6_GF3KF6_R_AD = []
KF6_GF3KF6_A_GT = []
KF6_GF3KF6_A_AD = []


# In[76]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KF6_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['6KF_R_GF3.GT'] == row['6KF_A_GF3.GT']) and (row['6KF_R_KF6.GT'] == row['6KF_A_KF6.GT']) and (row['6KF_R_GF3.GT'] == row['6KF_R_KF6.GT']) and (row['6KF_A_GF3.GT'] == row['6KF_A_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = (row['6KF_R_GF3.AD'] + row['6KF_R_KF6.AD'] + row['6KF_A_GF3.AD'] + row['6KF_A_KF6.AD'])/2
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_R_KF6.GT'])
            y = 0
            KF6_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['6KF_R_GF3.GT'] == row['6KF_A_KF6.GT']) and (row['6KF_A_GF3.GT'] == row['6KF_R_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = (row['6KF_R_GF3.AD'] + row['6KF_A_KF6.AD'])/2
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_A_GF3.GT'])
            y = (row['6KF_A_GF3.AD'] + row['6KF_R_KF6.AD'])/2
            KF6_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['6KF_R_GF3.GT'] == row['6KF_R_KF6.GT']) and (row['6KF_A_GF3.GT'] == row['6KF_A_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = (row['6KF_R_GF3.AD'] + row['6KF_R_KF6.AD'])/2
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_A_GF3.GT'])
            y = (row['6KF_A_GF3.AD'] + row['6KF_A_KF6.AD'])/2
            KF6_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['6KF_R_GF3.GT'] == row['6KF_A_GF3.GT']) and (row['6KF_R_KF6.GT'] != row['6KF_A_KF6.GT']) and (row['6KF_R_GF3.GT'] == row['6KF_R_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = (row['6KF_R_GF3.AD'] + row['6KF_A_GF3.AD'] + row['6KF_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_A_KF6.GT'])
            y = row['6KF_A_KF6.AD']
            KF6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['6KF_R_GF3.GT'] == row['6KF_A_GF3.GT']) and (row['6KF_R_KF6.GT'] != row['6KF_A_KF6.GT']) and (row['6KF_R_GF3.GT'] == row['6KF_A_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = (row['6KF_R_GF3.AD'] + row['6KF_A_GF3.AD'] + row['6KF_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_R_KF6.GT'])
            y = row['6KF_R_KF6.AD']
            KF6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['6KF_R_GF3.GT'] != row['6KF_A_GF3.GT']) and (row['6KF_R_KF6.GT'] == row['6KF_A_KF6.GT']) and (row['6KF_R_GF3.GT'] == row['6KF_R_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = row['6KF_R_GF3.AD'] 
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_A_GF3.GT'])
            y = (row['6KF_A_GF3.AD'] + row['6KF_R_KF6.AD'] + row['6KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['6KF_R_GF3.GT'] != row['6KF_A_GF3.GT']) and (row['6KF_R_KF6.GT'] == row['6KF_A_KF6.GT']) and (row['6KF_R_GF3.GT'] == row['6KF_A_KF6.GT'])):
            KF6_GF3KF6_R_GT.append(row['6KF_R_GF3.GT'])
            x = row['6KF_R_GF3.AD']
            KF6_GF3KF6_R_AD.append(x)
            KF6_GF3KF6_A_GT.append(row['6KF_R_KF6.GT'])
            y = (row['6KF_A_GF3.AD'] + row['6KF_R_KF6.AD'] + row['6KF_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KF6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KF6_GF3KF6_R_GT.append('error')
        KF6_GF3KF6_R_AD.append('error')
        KF6_GF3KF6_A_GT.append('error')
        KF6_GF3KF6_A_AD.append('error')


# In[77]:


df_backup['6KF_R_GF3KF6.GT'] = KF6_GF3KF6_R_GT
df_backup['6KF_A_GF3KF6.GT'] = KF6_GF3KF6_A_GT
df_backup['6KF_R_GF3KF6.AD'] = KF6_GF3KF6_R_AD
df_backup['6KF_A_GF3KF6.AD'] = KF6_GF3KF6_A_AD


# KS1, bad quality sample, we skip it

# KS2

# In[78]:


KS2_GF3KF6_R_GT = []
KS2_GF3KF6_R_AD = []
KS2_GF3KF6_A_GT = []
KS2_GF3KF6_A_AD = []


# In[79]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KS2_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['2KS_R_GF3.GT'] == row['2KS_A_GF3.GT']) and (row['2KS_R_KF6.GT'] == row['2KS_A_KF6.GT']) and (row['2KS_R_GF3.GT'] == row['2KS_R_KF6.GT']) and (row['2KS_A_GF3.GT'] == row['2KS_A_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = (row['2KS_R_GF3.AD'] + row['2KS_R_KF6.AD'] + row['2KS_A_GF3.AD'] + row['2KS_A_KF6.AD'])/2
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_R_KF6.GT'])
            y = 0
            KS2_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['2KS_R_GF3.GT'] == row['2KS_A_KF6.GT']) and (row['2KS_A_GF3.GT'] == row['2KS_R_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = (row['2KS_R_GF3.AD'] + row['2KS_A_KF6.AD'])/2
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_A_GF3.GT'])
            y = (row['2KS_A_GF3.AD'] + row['2KS_R_KF6.AD'])/2
            KS2_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['2KS_R_GF3.GT'] == row['2KS_R_KF6.GT']) and (row['2KS_A_GF3.GT'] == row['2KS_A_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = (row['2KS_R_GF3.AD'] + row['2KS_R_KF6.AD'])/2
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_A_GF3.GT'])
            y = (row['2KS_A_GF3.AD'] + row['2KS_A_KF6.AD'])/2
            KS2_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['2KS_R_GF3.GT'] == row['2KS_A_GF3.GT']) and (row['2KS_R_KF6.GT'] != row['2KS_A_KF6.GT']) and (row['2KS_R_GF3.GT'] == row['2KS_R_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = (row['2KS_R_GF3.AD'] + row['2KS_A_GF3.AD'] + row['2KS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_A_KF6.GT'])
            y = row['2KS_A_KF6.AD']
            KS2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['2KS_R_GF3.GT'] == row['2KS_A_GF3.GT']) and (row['2KS_R_KF6.GT'] != row['2KS_A_KF6.GT']) and (row['2KS_R_GF3.GT'] == row['2KS_A_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = (row['2KS_R_GF3.AD'] + row['2KS_A_GF3.AD'] + row['2KS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_R_KF6.GT'])
            y = row['2KS_R_KF6.AD']
            KS2_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['2KS_R_GF3.GT'] != row['2KS_A_GF3.GT']) and (row['2KS_R_KF6.GT'] == row['2KS_A_KF6.GT']) and (row['2KS_R_GF3.GT'] == row['2KS_R_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = row['2KS_R_GF3.AD'] 
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_A_GF3.GT'])
            y = (row['2KS_A_GF3.AD'] + row['2KS_R_KF6.AD'] + row['2KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['2KS_R_GF3.GT'] != row['2KS_A_GF3.GT']) and (row['2KS_R_KF6.GT'] == row['2KS_A_KF6.GT']) and (row['2KS_R_GF3.GT'] == row['2KS_A_KF6.GT'])):
            KS2_GF3KF6_R_GT.append(row['2KS_R_GF3.GT'])
            x = row['2KS_R_GF3.AD']
            KS2_GF3KF6_R_AD.append(x)
            KS2_GF3KF6_A_GT.append(row['2KS_R_KF6.GT'])
            y = (row['2KS_A_GF3.AD'] + row['2KS_R_KF6.AD'] + row['2KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS2_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KS2_GF3KF6_R_GT.append('error')
        KS2_GF3KF6_R_AD.append('error')
        KS2_GF3KF6_A_GT.append('error')
        KS2_GF3KF6_A_AD.append('error')


# In[80]:


df_backup['2KS_R_GF3KF6.GT'] = KS2_GF3KF6_R_GT
df_backup['2KS_A_GF3KF6.GT'] = KS2_GF3KF6_A_GT
df_backup['2KS_R_GF3KF6.AD'] = KS2_GF3KF6_R_AD
df_backup['2KS_A_GF3KF6.AD'] = KS2_GF3KF6_A_AD


# KS3

# In[81]:


KS3_GF3KF6_R_GT = []
KS3_GF3KF6_R_AD = []
KS3_GF3KF6_A_GT = []
KS3_GF3KF6_A_AD = []


# In[82]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KS3_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['3KS_R_GF3.GT'] == row['3KS_A_GF3.GT']) and (row['3KS_R_KF6.GT'] == row['3KS_A_KF6.GT']) and (row['3KS_R_GF3.GT'] == row['3KS_R_KF6.GT']) and (row['3KS_A_GF3.GT'] == row['3KS_A_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = (row['3KS_R_GF3.AD'] + row['3KS_R_KF6.AD'] + row['3KS_A_GF3.AD'] + row['3KS_A_KF6.AD'])/2
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_R_KF6.GT'])
            y = 0
            KS3_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['3KS_R_GF3.GT'] == row['3KS_A_KF6.GT']) and (row['3KS_A_GF3.GT'] == row['3KS_R_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = (row['3KS_R_GF3.AD'] + row['3KS_A_KF6.AD'])/2
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_A_GF3.GT'])
            y = (row['3KS_A_GF3.AD'] + row['3KS_R_KF6.AD'])/2
            KS3_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['3KS_R_GF3.GT'] == row['3KS_R_KF6.GT']) and (row['3KS_A_GF3.GT'] == row['3KS_A_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = (row['3KS_R_GF3.AD'] + row['3KS_R_KF6.AD'])/2
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_A_GF3.GT'])
            y = (row['3KS_A_GF3.AD'] + row['3KS_A_KF6.AD'])/2
            KS3_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['3KS_R_GF3.GT'] == row['3KS_A_GF3.GT']) and (row['3KS_R_KF6.GT'] != row['3KS_A_KF6.GT']) and (row['3KS_R_GF3.GT'] == row['3KS_R_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = (row['3KS_R_GF3.AD'] + row['3KS_A_GF3.AD'] + row['3KS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_A_KF6.GT'])
            y = row['3KS_A_KF6.AD']
            KS3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['3KS_R_GF3.GT'] == row['3KS_A_GF3.GT']) and (row['3KS_R_KF6.GT'] != row['3KS_A_KF6.GT']) and (row['3KS_R_GF3.GT'] == row['3KS_A_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = (row['3KS_R_GF3.AD'] + row['3KS_A_GF3.AD'] + row['3KS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_R_KF6.GT'])
            y = row['3KS_R_KF6.AD']
            KS3_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['3KS_R_GF3.GT'] != row['3KS_A_GF3.GT']) and (row['3KS_R_KF6.GT'] == row['3KS_A_KF6.GT']) and (row['3KS_R_GF3.GT'] == row['3KS_R_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = row['3KS_R_GF3.AD'] 
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_A_GF3.GT'])
            y = (row['3KS_A_GF3.AD'] + row['3KS_R_KF6.AD'] + row['3KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['3KS_R_GF3.GT'] != row['3KS_A_GF3.GT']) and (row['3KS_R_KF6.GT'] == row['3KS_A_KF6.GT']) and (row['3KS_R_GF3.GT'] == row['3KS_A_KF6.GT'])):
            KS3_GF3KF6_R_GT.append(row['3KS_R_GF3.GT'])
            x = row['3KS_R_GF3.AD']
            KS3_GF3KF6_R_AD.append(x)
            KS3_GF3KF6_A_GT.append(row['3KS_R_KF6.GT'])
            y = (row['3KS_A_GF3.AD'] + row['3KS_R_KF6.AD'] + row['3KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS3_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KS3_GF3KF6_R_GT.append('error')
        KS3_GF3KF6_R_AD.append('error')
        KS3_GF3KF6_A_GT.append('error')
        KS3_GF3KF6_A_AD.append('error')


# In[83]:


df_backup['3KS_R_GF3KF6.GT'] = KS3_GF3KF6_R_GT
df_backup['3KS_A_GF3KF6.GT'] = KS3_GF3KF6_A_GT
df_backup['3KS_R_GF3KF6.AD'] = KS3_GF3KF6_R_AD
df_backup['3KS_A_GF3KF6.AD'] = KS3_GF3KF6_A_AD


# KS4

# In[84]:


KS4_GF3KF6_R_GT = []
KS4_GF3KF6_R_AD = []
KS4_GF3KF6_A_GT = []
KS4_GF3KF6_A_AD = []


# In[85]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KS4_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['4KS_R_GF3.GT'] == row['4KS_A_GF3.GT']) and (row['4KS_R_KF6.GT'] == row['4KS_A_KF6.GT']) and (row['4KS_R_GF3.GT'] == row['4KS_R_KF6.GT']) and (row['4KS_A_GF3.GT'] == row['4KS_A_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = (row['4KS_R_GF3.AD'] + row['4KS_R_KF6.AD'] + row['4KS_A_GF3.AD'] + row['4KS_A_KF6.AD'])/2
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_R_KF6.GT'])
            y = 0
            KS4_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['4KS_R_GF3.GT'] == row['4KS_A_KF6.GT']) and (row['4KS_A_GF3.GT'] == row['4KS_R_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = (row['4KS_R_GF3.AD'] + row['4KS_A_KF6.AD'])/2
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_A_GF3.GT'])
            y = (row['4KS_A_GF3.AD'] + row['4KS_R_KF6.AD'])/2
            KS4_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['4KS_R_GF3.GT'] == row['4KS_R_KF6.GT']) and (row['4KS_A_GF3.GT'] == row['4KS_A_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = (row['4KS_R_GF3.AD'] + row['4KS_R_KF6.AD'])/2
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_A_GF3.GT'])
            y = (row['4KS_A_GF3.AD'] + row['4KS_A_KF6.AD'])/2
            KS4_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['4KS_R_GF3.GT'] == row['4KS_A_GF3.GT']) and (row['4KS_R_KF6.GT'] != row['4KS_A_KF6.GT']) and (row['4KS_R_GF3.GT'] == row['4KS_R_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = (row['4KS_R_GF3.AD'] + row['4KS_A_GF3.AD'] + row['4KS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_A_KF6.GT'])
            y = row['4KS_A_KF6.AD']
            KS4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['4KS_R_GF3.GT'] == row['4KS_A_GF3.GT']) and (row['4KS_R_KF6.GT'] != row['4KS_A_KF6.GT']) and (row['4KS_R_GF3.GT'] == row['4KS_A_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = (row['4KS_R_GF3.AD'] + row['4KS_A_GF3.AD'] + row['4KS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_R_KF6.GT'])
            y = row['4KS_R_KF6.AD']
            KS4_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['4KS_R_GF3.GT'] != row['4KS_A_GF3.GT']) and (row['4KS_R_KF6.GT'] == row['4KS_A_KF6.GT']) and (row['4KS_R_GF3.GT'] == row['4KS_R_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = row['4KS_R_GF3.AD'] 
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_A_GF3.GT'])
            y = (row['4KS_A_GF3.AD'] + row['4KS_R_KF6.AD'] + row['4KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['4KS_R_GF3.GT'] != row['4KS_A_GF3.GT']) and (row['4KS_R_KF6.GT'] == row['4KS_A_KF6.GT']) and (row['4KS_R_GF3.GT'] == row['4KS_A_KF6.GT'])):
            KS4_GF3KF6_R_GT.append(row['4KS_R_GF3.GT'])
            x = row['4KS_R_GF3.AD']
            KS4_GF3KF6_R_AD.append(x)
            KS4_GF3KF6_A_GT.append(row['4KS_R_KF6.GT'])
            y = (row['4KS_A_GF3.AD'] + row['4KS_R_KF6.AD'] + row['4KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS4_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KS4_GF3KF6_R_GT.append('error')
        KS4_GF3KF6_R_AD.append('error')
        KS4_GF3KF6_A_GT.append('error')
        KS4_GF3KF6_A_AD.append('error')


# In[86]:


df_backup['4KS_R_GF3KF6.GT'] = KS4_GF3KF6_R_GT
df_backup['4KS_A_GF3KF6.GT'] = KS4_GF3KF6_A_GT
df_backup['4KS_R_GF3KF6.AD'] = KS4_GF3KF6_R_AD
df_backup['4KS_A_GF3KF6.AD'] = KS4_GF3KF6_A_AD


# KS5

# In[87]:


KS5_GF3KF6_R_GT = []
KS5_GF3KF6_R_AD = []
KS5_GF3KF6_A_GT = []
KS5_GF3KF6_A_AD = []


# In[88]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KS5_GF3KF6_R_GT),'ref ',row['5KS_R_GF3.GT'],row['5KS_R_KF6.GT'],'alt',row['5KS_A_GF3.GT'],row['5KS_A_KF6.GT'])'''
    '''Case heterozygots'''
    if ((row['5KS_R_GF3.GT'] == row['5KS_A_GF3.GT']) and (row['5KS_R_KF6.GT'] == row['5KS_A_KF6.GT']) and (row['5KS_R_GF3.GT'] == row['5KS_R_KF6.GT']) and (row['5KS_A_GF3.GT'] == row['5KS_A_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = (row['5KS_R_GF3.AD'] + row['5KS_R_KF6.AD'] + row['5KS_A_GF3.AD'] + row['5KS_A_KF6.AD'])/2
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_R_KF6.GT'])
            y = 0
            KS5_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            ''''print('Both homozgots',i)'''
    elif ((row['5KS_R_GF3.GT'] == row['5KS_A_KF6.GT']) and (row['5KS_A_GF3.GT'] == row['5KS_R_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = (row['5KS_R_GF3.AD'] + row['5KS_A_KF6.AD'])/2
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_A_GF3.GT'])
            y = (row['5KS_A_GF3.AD'] + row['5KS_R_KF6.AD'])/2
            KS5_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['5KS_R_GF3.GT'] == row['5KS_R_KF6.GT']) and (row['5KS_A_GF3.GT'] == row['5KS_A_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = (row['5KS_R_GF3.AD'] + row['5KS_R_KF6.AD'])/2
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_A_GF3.GT'])
            y = (row['5KS_A_GF3.AD'] + row['5KS_A_KF6.AD'])/2
            KS5_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['5KS_R_GF3.GT'] == row['5KS_A_GF3.GT']) and (row['5KS_R_KF6.GT'] != row['5KS_A_KF6.GT']) and (row['5KS_R_GF3.GT'] == row['5KS_R_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = (row['5KS_R_GF3.AD'] + row['5KS_A_GF3.AD'] + row['5KS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_A_KF6.GT'])
            y = row['5KS_A_KF6.AD']
            KS5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['5KS_R_GF3.GT'] == row['5KS_A_GF3.GT']) and (row['5KS_R_KF6.GT'] != row['5KS_A_KF6.GT']) and (row['5KS_R_GF3.GT'] == row['5KS_A_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = (row['5KS_R_GF3.AD'] + row['5KS_A_GF3.AD'] + row['5KS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_R_KF6.GT'])
            y = row['5KS_R_KF6.AD']
            KS5_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['5KS_R_GF3.GT'] != row['5KS_A_GF3.GT']) and (row['5KS_R_KF6.GT'] == row['5KS_A_KF6.GT']) and (row['5KS_R_GF3.GT'] == row['5KS_R_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = row['5KS_R_GF3.AD'] 
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_A_GF3.GT'])
            y = (row['5KS_A_GF3.AD'] + row['5KS_R_KF6.AD'] + row['5KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['5KS_R_GF3.GT'] != row['5KS_A_GF3.GT']) and (row['5KS_R_KF6.GT'] == row['5KS_A_KF6.GT']) and (row['5KS_R_GF3.GT'] == row['5KS_A_KF6.GT'])):
            KS5_GF3KF6_R_GT.append(row['5KS_R_GF3.GT'])
            x = row['5KS_R_GF3.AD']
            KS5_GF3KF6_R_AD.append(x)
            KS5_GF3KF6_A_GT.append(row['5KS_R_KF6.GT'])
            y = (row['5KS_A_GF3.AD'] + row['5KS_R_KF6.AD'] + row['5KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS5_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KS5_GF3KF6_R_GT.append('error')
        KS5_GF3KF6_R_AD.append('error')
        KS5_GF3KF6_A_GT.append('error')
        KS5_GF3KF6_A_AD.append('error')


# In[89]:


df_backup['5KS_R_GF3KF6.GT'] = KS5_GF3KF6_R_GT
df_backup['5KS_A_GF3KF6.GT'] = KS5_GF3KF6_A_GT
df_backup['5KS_R_GF3KF6.AD'] = KS5_GF3KF6_R_AD
df_backup['5KS_A_GF3KF6.AD'] = KS5_GF3KF6_A_AD


# KS6

# In[90]:


KS6_GF3KF6_R_GT = []
KS6_GF3KF6_R_AD = []
KS6_GF3KF6_A_GT = []
KS6_GF3KF6_A_AD = []


# In[91]:


for index, row in df_bi.iterrows() : 
    i = index
    '''print ('index',i,'length',len(KS6_GF3KF6_R_GT))'''
    '''Case heterozygots'''
    if ((row['6KS_R_GF3.GT'] == row['6KS_A_GF3.GT']) and (row['6KS_R_KF6.GT'] == row['6KS_A_KF6.GT']) and (row['6KS_R_GF3.GT'] == row['6KS_R_KF6.GT']) and (row['6KS_A_GF3.GT'] == row['6KS_A_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = (row['6KS_R_GF3.AD'] + row['6KS_R_KF6.AD'] + row['6KS_A_GF3.AD'] + row['6KS_A_KF6.AD'])/2
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_R_KF6.GT'])
            y = 0
            KS6_GF3KF6_A_AD.append(y)
            Homozygous.append(i)
            '''print('Both homozgots',i)'''
    elif ((row['6KS_R_GF3.GT'] == row['6KS_A_KF6.GT']) and (row['6KS_A_GF3.GT'] == row['6KS_R_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = (row['6KS_R_GF3.AD'] + row['6KS_A_KF6.AD'])/2
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_A_GF3.GT'])
            y = (row['6KS_A_GF3.AD'] + row['6KS_R_KF6.AD'])/2
            KS6_GF3KF6_A_AD.append(y)
            '''print('Inverse situation Ref and Alt in both for row',i)'''
    elif ((row['6KS_R_GF3.GT'] == row['6KS_R_KF6.GT']) and (row['6KS_A_GF3.GT'] == row['6KS_A_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = (row['6KS_R_GF3.AD'] + row['6KS_R_KF6.AD'])/2
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_A_GF3.GT'])
            y = (row['6KS_A_GF3.AD'] + row['6KS_A_KF6.AD'])/2
            KS6_GF3KF6_A_AD.append(y)
            '''print('Same situation Ref and Alt in both for row',i)'''
    elif ((row['6KS_R_GF3.GT'] == row['6KS_A_GF3.GT']) and (row['6KS_R_KF6.GT'] != row['6KS_A_KF6.GT']) and (row['6KS_R_GF3.GT'] == row['6KS_R_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = (row['6KS_R_GF3.AD'] + row['6KS_A_GF3.AD'] + row['6KS_R_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_A_KF6.GT'])
            y = row['6KS_A_KF6.AD']
            KS6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in second position',i)'''
    elif ((row['6KS_R_GF3.GT'] == row['6KS_A_GF3.GT']) and (row['6KS_R_KF6.GT'] != row['6KS_A_KF6.GT']) and (row['6KS_R_GF3.GT'] == row['6KS_A_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = (row['6KS_R_GF3.AD'] + row['6KS_A_GF3.AD'] + row['6KS_A_KF6.AD'])/2 # counts from 3 cells, from 2 individuals
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_R_KF6.GT'])
            y = row['6KS_R_KF6.AD']
            KS6_GF3KF6_A_AD.append(y)
            '''print('First homozygot, second with alternative allele in first position vof the second genome',i)'''
    elif ((row['6KS_R_GF3.GT'] != row['6KS_A_GF3.GT']) and (row['6KS_R_KF6.GT'] == row['6KS_A_KF6.GT']) and (row['6KS_R_GF3.GT'] == row['6KS_R_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = row['6KS_R_GF3.AD'] 
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_A_GF3.GT'])
            y = (row['6KS_A_GF3.AD'] + row['6KS_R_KF6.AD'] + row['6KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in second position of the first genome',i)'''
    elif ((row['6KS_R_GF3.GT'] != row['6KS_A_GF3.GT']) and (row['6KS_R_KF6.GT'] == row['6KS_A_KF6.GT']) and (row['6KS_R_GF3.GT'] == row['6KS_A_KF6.GT'])):
            KS6_GF3KF6_R_GT.append(row['6KS_R_GF3.GT'])
            x = row['6KS_R_GF3.AD']
            KS6_GF3KF6_R_AD.append(x)
            KS6_GF3KF6_A_GT.append(row['6KS_R_KF6.GT'])
            y = (row['6KS_A_GF3.AD'] + row['6KS_R_KF6.AD'] + row['6KS_A_KF6.AD'])/2  # counts from 3 cells, from 2 individuals
            KS6_GF3KF6_A_AD.append(y)
            '''print('Second homozygot, first with alternative allele in first position of the first genome',i)'''
    else:
        print('Error in row ', i)
        KS6_GF3KF6_R_GT.append('error')
        KS6_GF3KF6_R_AD.append('error')
        KS6_GF3KF6_A_GT.append('error')
        KS6_GF3KF6_A_AD.append('error')


# In[92]:


df_backup['6KS_R_GF3KF6.GT'] = KS6_GF3KF6_R_GT
df_backup['6KS_A_GF3KF6.GT'] = KS6_GF3KF6_A_GT
df_backup['6KS_R_GF3KF6.AD'] = KS6_GF3KF6_R_AD
df_backup['6KS_A_GF3KF6.AD'] = KS6_GF3KF6_A_AD


# In[93]:


df_backup.head(20)


# In[94]:


df_backup.to_csv('GF3KF6.csv')
np.savetxt('Homozygous.csv', Homozygous, delimiter=",", fmt='%s')
#print(Homozygous)


# For a fresh start in case the Kernel is dead afer the analysis

# In[95]:


import pandas as pd
import numpy as np


# In[96]:


df_backup = pd.read_csv('GF3KF6.csv')


# ### Final dataframe after mixing both mappings

# Select the columns and clean up values you don't need

# In[97]:


for col in df_backup.columns: 
    print(col) 


# In[98]:


def select_columns(data_frame, column_names):
    new_frame = data_frame.loc[:, column_names]
    return new_frame


# In[99]:


selected_columns = ['CHROM',
                    'POS',
                    'Gene.refGene_x',
                    'Func.refGene_x',
                    'ExonicFunc.refGene_x',
                    'AF_x',
                    'Gene.refGene_y',
                    'Func.refGene_y',
                    'ExonicFunc.refGene_y',
                    'AF_y',
                    '1GF_R_GF3KF6.GT',
                    '1GF_A_GF3KF6.GT',
                    '1GF_R_GF3KF6.AD',
                    '1GF_A_GF3KF6.AD',
                    '2GF_R_GF3KF6.GT',
                    '2GF_A_GF3KF6.GT',
                    '2GF_R_GF3KF6.AD',
                    '2GF_A_GF3KF6.AD',
                    '3GF_R_GF3KF6.GT',
                    '3GF_A_GF3KF6.GT',
                    '3GF_R_GF3KF6.AD',
                    '3GF_A_GF3KF6.AD',
                    '4GF_R_GF3KF6.GT',
                    '4GF_A_GF3KF6.GT',
                    '4GF_R_GF3KF6.AD',
                    '4GF_A_GF3KF6.AD',
                    '5GF_R_GF3KF6.GT',
                    '5GF_A_GF3KF6.GT',
                    '5GF_R_GF3KF6.AD',
                    '5GF_A_GF3KF6.AD',
                    '6GF_R_GF3KF6.GT',
                    '6GF_A_GF3KF6.GT',
                    '6GF_R_GF3KF6.AD',
                    '6GF_A_GF3KF6.AD',
                    '1GS_R_GF3KF6.GT',
                    '1GS_A_GF3KF6.GT',
                    '1GS_R_GF3KF6.AD',
                    '1GS_A_GF3KF6.AD',
                    '2GS_R_GF3KF6.GT',
                    '2GS_A_GF3KF6.GT',
                    '2GS_R_GF3KF6.AD',
                    '2GS_A_GF3KF6.AD',
                    '3GS_R_GF3KF6.GT',
                    '3GS_A_GF3KF6.GT',
                    '3GS_R_GF3KF6.AD',
                    '3GS_A_GF3KF6.AD',
                    '4GS_R_GF3KF6.GT',
                    '4GS_A_GF3KF6.GT',
                    '4GS_R_GF3KF6.AD',
                    '4GS_A_GF3KF6.AD',
                    '5GS_R_GF3KF6.GT',
                    '5GS_A_GF3KF6.GT',
                    '5GS_R_GF3KF6.AD',
                    '5GS_A_GF3KF6.AD',
                    '6GS_R_GF3KF6.GT',
                    '6GS_A_GF3KF6.GT',
                    '6GS_R_GF3KF6.AD',
                    '6GS_A_GF3KF6.AD',
                    '1KF_R_GF3KF6.GT',
                    '1KF_A_GF3KF6.GT',
                    '1KF_R_GF3KF6.AD',
                    '1KF_A_GF3KF6.AD',
                    '2KF_R_GF3KF6.GT',
                    '2KF_A_GF3KF6.GT',
                    '2KF_R_GF3KF6.AD',
                    '2KF_A_GF3KF6.AD',
                    '3KF_R_GF3KF6.GT',
                    '3KF_A_GF3KF6.GT',
                    '3KF_R_GF3KF6.AD',
                    '3KF_A_GF3KF6.AD',
                    '4KF_R_GF3KF6.GT',
                    '4KF_A_GF3KF6.GT',
                    '4KF_R_GF3KF6.AD',
                    '4KF_A_GF3KF6.AD',
                    '5KF_R_GF3KF6.GT',
                    '5KF_A_GF3KF6.GT',
                    '5KF_R_GF3KF6.AD',
                    '5KF_A_GF3KF6.AD',
                    '6KF_R_GF3KF6.GT',
                    '6KF_A_GF3KF6.GT',
                    '6KF_R_GF3KF6.AD',
                    '6KF_A_GF3KF6.AD',
                    '2KS_R_GF3KF6.GT',
                    '2KS_A_GF3KF6.GT',
                    '2KS_R_GF3KF6.AD',
                    '2KS_A_GF3KF6.AD',
                    '3KS_R_GF3KF6.GT',
                    '3KS_A_GF3KF6.GT',
                    '3KS_R_GF3KF6.AD',
                    '3KS_A_GF3KF6.AD',
                    '4KS_R_GF3KF6.GT',
                    '4KS_A_GF3KF6.GT',
                    '4KS_R_GF3KF6.AD',
                    '4KS_A_GF3KF6.AD',
                    '5KS_R_GF3KF6.GT',
                    '5KS_A_GF3KF6.GT',
                    '5KS_R_GF3KF6.AD',
                    '5KS_A_GF3KF6.AD',
                    '6KS_R_GF3KF6.GT',
                    '6KS_A_GF3KF6.GT',
                    '6KS_R_GF3KF6.AD',
                    '6KS_A_GF3KF6.AD'] #Put the name of the columns you want to conserve
df_bi_final = select_columns(df_backup, selected_columns)


# Make a backup of the final file

# In[100]:


df_bi_final.to_csv('Allele_counts_2maps.csv')


# In[101]:


df = df_bi_final


# In[102]:


len(df_bi_final)


# ## Clean low AD

# Now that we have the averaged and unbiased SNPs for our sample population, we need to make a filter for allele counts that are superior to 10 in all the columns of ".AD" (allele depth). To do so, sum the counts on each SNP for the reference and alternative alleles ad drop the SNPs below 10 counts.

# In case the kernel is dead, loead the files for a fresh start:

# In[103]:


import pandas as pd
import numpy as np


# In[160]:


df = pd.read_csv('Allele_counts_2maps.csv')


# In[161]:


index_df1GF = df[(df['1GF_R_GF3KF6.AD']+df['1GF_A_GF3KF6.AD']) < 10].index
index_df2GF = df[(df['2GF_R_GF3KF6.AD']+df['2GF_A_GF3KF6.AD']) < 10].index
index_df3GF = df[(df['3GF_R_GF3KF6.AD']+df['3GF_A_GF3KF6.AD']) < 10].index
index_df4GF = df[(df['4GF_R_GF3KF6.AD']+df['4GF_A_GF3KF6.AD']) < 10].index
index_df5GF = df[(df['5GF_R_GF3KF6.AD']+df['5GF_A_GF3KF6.AD']) < 10].index
index_df6GF = df[(df['6GF_R_GF3KF6.AD']+df['6GF_A_GF3KF6.AD']) < 10].index
index_df1GS = df[(df['1GS_R_GF3KF6.AD']+df['1GS_A_GF3KF6.AD']) < 10].index
index_df2GS = df[(df['2GS_R_GF3KF6.AD']+df['2GS_A_GF3KF6.AD']) < 10].index
index_df3GS = df[(df['3GS_R_GF3KF6.AD']+df['3GS_A_GF3KF6.AD']) < 10].index
index_df4GS = df[(df['4GS_R_GF3KF6.AD']+df['4GS_A_GF3KF6.AD']) < 10].index
index_df5GS = df[(df['5GS_R_GF3KF6.AD']+df['5GS_A_GF3KF6.AD']) < 10].index
index_df6GS = df[(df['6GS_R_GF3KF6.AD']+df['6GS_A_GF3KF6.AD']) < 10].index
index_df1KF = df[(df['1KF_R_GF3KF6.AD']+df['1KF_A_GF3KF6.AD']) < 10].index
index_df2KF = df[(df['2KF_R_GF3KF6.AD']+df['2KF_A_GF3KF6.AD']) < 10].index
index_df3KF = df[(df['3KF_R_GF3KF6.AD']+df['3KF_A_GF3KF6.AD']) < 10].index
index_df4KF = df[(df['4KF_R_GF3KF6.AD']+df['4KF_A_GF3KF6.AD']) < 10].index
index_df5KF = df[(df['5KF_R_GF3KF6.AD']+df['5KF_A_GF3KF6.AD']) < 10].index
index_df6KF = df[(df['6KF_R_GF3KF6.AD']+df['6KF_A_GF3KF6.AD']) < 10].index
index_df2KS = df[(df['2KS_R_GF3KF6.AD']+df['2KS_A_GF3KF6.AD']) < 10].index
index_df3KS = df[(df['3KS_R_GF3KF6.AD']+df['3KS_A_GF3KF6.AD']) < 10].index
index_df4KS = df[(df['4KS_R_GF3KF6.AD']+df['4KS_A_GF3KF6.AD']) < 10].index
index_df5KS = df[(df['5KS_R_GF3KF6.AD']+df['5KS_A_GF3KF6.AD']) < 10].index
index_df6KS = df[(df['6KS_R_GF3KF6.AD']+df['6KS_A_GF3KF6.AD']) < 10].index


# Keep all the indexes that are >= 10 for the SNPs in the reference genome and make a new matrix 

# In[164]:


AD10_df = index_df1GF.intersection(index_df1GS)
AD10_df = AD10_df.intersection(index_df1KF)
AD10_df = AD10_df.intersection(index_df2GF)
AD10_df = AD10_df.intersection(index_df2GS)
AD10_df = AD10_df.intersection(index_df2KF)
AD10_df = AD10_df.intersection(index_df2KS)
AD10_df = AD10_df.intersection(index_df3GF)
AD10_df = AD10_df.intersection(index_df3GS)
AD10_df = AD10_df.intersection(index_df3KF)
AD10_df = AD10_df.intersection(index_df3KS)
AD10_df = AD10_df.intersection(index_df4GF)
AD10_df = AD10_df.intersection(index_df4GS)
AD10_df = AD10_df.intersection(index_df4KF)
AD10_df = AD10_df.intersection(index_df4KS)
AD10_df = AD10_df.intersection(index_df5GF)
AD10_df = AD10_df.intersection(index_df5GS)
AD10_df = AD10_df.intersection(index_df5KF)
AD10_df = AD10_df.intersection(index_df5KS)
AD10_df = AD10_df.intersection(index_df6GF)
AD10_df = AD10_df.intersection(index_df6GS)
AD10_df = AD10_df.intersection(index_df6KF)
AD10_df = AD10_df.intersection(index_df6KS)


# Make the new DFs with the >= 10 AD

# In[165]:


df10 = df.drop(index = AD10_df)  #df10 = df.loc[AD10_df,:]


# In[166]:


df10.columns


# In[167]:


len(df10)


# Make an SNP ID based on the index

# In[168]:


df10['SNP_ID']=df10.index


# In[169]:


df10.to_csv('SNPs_2map_clean.csv')


# # Assign reference and alternative alleles uniformly

# The reference and alternative allele in each of the different samples have been assigned according to a comparison between the different mappings of the same sample, but this process has been independent. Now the ref and alt alleles from each sample need to be assigned properly according to the other samples. The model will be the sample KF6 from which the alternative genome was constructed.

# In[170]:


import pandas as pd
import numpy as np


# In[171]:


df10 = pd.read_csv('SNPs_2map_clean.csv')


# In[172]:


df10.tail(10)


# In[173]:


len(df10)


# KF6

# In[174]:


KF6R_AD = []
KF6A_AD = []
KF6R_GT = []
KF6A_GT = []


# In[175]:


for index, row in df10.iterrows() : 
    i = index
    if (row['6KF_R_GF3KF6.GT'] != row['6KF_R_GF3KF6.GT']):
        KF6R_GT.append(row['6KF_A_GF3KF6.GT'])
        KF6A_GT.append(row['6KF_R_GF3KF6.GT'])
        KF6R_AD.append(row['6KF_A_GF3KF6.AD'])
        KF6A_AD.append(row['6KF_R_GF3KF6.AD'])   
    else :
        KF6R_GT.append(row['6KF_R_GF3KF6.GT'])
        KF6A_GT.append(row['6KF_A_GF3KF6.GT'])
        KF6R_AD.append(row['6KF_R_GF3KF6.AD'])
        KF6A_AD.append(row['6KF_A_GF3KF6.AD'])
        
print(i)


# In[176]:


GF3R_AD = []
GF3A_AD = []
GF3R_GT = []
GF3A_GT = []


# In[177]:


for index, row in df10.iterrows() : 
    i = index
    if (row['6KF_R_GF3KF6.GT'] != row['3GF_R_GF3KF6.GT']):
        GF3R_GT.append(row['3GF_A_GF3KF6.GT'])
        GF3A_GT.append(row['3GF_R_GF3KF6.GT'])
        GF3R_AD.append(row['3GF_A_GF3KF6.AD'])
        GF3A_AD.append(row['3GF_R_GF3KF6.AD'])   
    else :
        GF3R_GT.append(row['3GF_R_GF3KF6.GT'])
        GF3A_GT.append(row['3GF_A_GF3KF6.GT'])
        GF3R_AD.append(row['3GF_R_GF3KF6.AD'])
        GF3A_AD.append(row['3GF_A_GF3KF6.AD'])
        
print(i)


# In[178]:


len(GF3R_AD)


# Starting point

# GF1

# In[179]:


GF1R_AD = []
GF1A_AD = []
GF1R_GT = []
GF1A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['1GF_R_GF3KF6.GT']):
        GF1R_GT.append(row['1GF_A_GF3KF6.GT'])
        GF1A_GT.append(row['1GF_R_GF3KF6.GT'])
        GF1R_AD.append(row['1GF_A_GF3KF6.AD'])
        GF1A_AD.append(row['1GF_R_GF3KF6.AD'])
    else :
        GF1R_GT.append(row['1GF_R_GF3KF6.GT'])
        GF1A_GT.append(row['1GF_A_GF3KF6.GT'])
        GF1R_AD.append(row['1GF_R_GF3KF6.AD'])
        GF1A_AD.append(row['1GF_A_GF3KF6.AD'])
print(i)


# GF2

# In[180]:


GF2R_AD = []
GF2A_AD = []
GF2R_GT = []
GF2A_GT = []


# In[181]:


for index, row in df10.iterrows() : 
    i = index
    if (row['6KF_R_GF3KF6.GT'] != row['2GF_R_GF3KF6.GT']):
        GF2R_GT.append(row['2GF_A_GF3KF6.GT'])
        GF2A_GT.append(row['2GF_R_GF3KF6.GT'])
        GF2R_AD.append(row['2GF_A_GF3KF6.AD'])
        GF2A_AD.append(row['2GF_R_GF3KF6.AD'])   
    else :
        GF2R_GT.append(row['2GF_R_GF3KF6.GT'])
        GF2A_GT.append(row['2GF_A_GF3KF6.GT'])
        GF2R_AD.append(row['2GF_R_GF3KF6.AD'])
        GF2A_AD.append(row['2GF_A_GF3KF6.AD'])
print(i)


# GF4

# In[182]:


GF4R_AD = []
GF4A_AD = []
GF4R_GT = []
GF4A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['4GF_R_GF3KF6.GT']):
        GF4R_GT.append(row['4GF_A_GF3KF6.GT'])
        GF4A_GT.append(row['4GF_R_GF3KF6.GT'])
        GF4R_AD.append(row['4GF_A_GF3KF6.AD'])
        GF4A_AD.append(row['4GF_R_GF3KF6.AD'])
    else :
        GF4R_GT.append(row['4GF_R_GF3KF6.GT'])
        GF4A_GT.append(row['4GF_A_GF3KF6.GT'])
        GF4R_AD.append(row['4GF_R_GF3KF6.AD'])
        GF4A_AD.append(row['4GF_A_GF3KF6.AD'])
print(i)


# GF5

# In[183]:


GF5R_AD = []
GF5A_AD = []
GF5R_GT = []
GF5A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['5GF_R_GF3KF6.GT']):
        GF5R_GT.append(row['5GF_A_GF3KF6.GT'])
        GF5A_GT.append(row['5GF_R_GF3KF6.GT'])
        GF5R_AD.append(row['5GF_A_GF3KF6.AD'])
        GF5A_AD.append(row['5GF_R_GF3KF6.AD'])
   
    else :
        GF5R_GT.append(row['5GF_R_GF3KF6.GT'])
        GF5A_GT.append(row['5GF_A_GF3KF6.GT'])
        GF5R_AD.append(row['5GF_R_GF3KF6.AD'])
        GF5A_AD.append(row['5GF_A_GF3KF6.AD'])
print(i)


# GF6

# In[184]:


GF6R_AD = []
GF6A_AD = []
GF6R_GT = []
GF6A_GT = []
for index, row in df10.iterrows() : 
    i = index
    if (row['6KF_R_GF3KF6.GT'] != row['6GF_R_GF3KF6.GT']):
        GF6R_GT.append(row['6GF_A_GF3KF6.GT'])
        GF6A_GT.append(row['6GF_R_GF3KF6.GT'])
        GF6R_AD.append(row['6GF_A_GF3KF6.AD'])
        GF6A_AD.append(row['6GF_R_GF3KF6.AD'])
   
    else :
        GF6R_GT.append(row['6GF_R_GF3KF6.GT'])
        GF6A_GT.append(row['6GF_A_GF3KF6.GT'])
        GF6R_AD.append(row['6GF_R_GF3KF6.AD'])
        GF6A_AD.append(row['6GF_A_GF3KF6.AD'])
print(i)  


# GS1

# In[185]:


GS1R_AD = []
GS1A_AD = []
GS1R_GT = []
GS1A_GT = []
for index, row in df10.iterrows() : 
    i = index
    if (row['6KF_R_GF3KF6.GT'] != row['1GS_R_GF3KF6.GT']):
        GS1R_GT.append(row['1GS_A_GF3KF6.GT'])
        GS1A_GT.append(row['1GS_R_GF3KF6.GT'])
        GS1R_AD.append(row['1GS_A_GF3KF6.AD'])
        GS1A_AD.append(row['1GS_R_GF3KF6.AD'])
    else :
        GS1R_GT.append(row['1GS_R_GF3KF6.GT'])
        GS1A_GT.append(row['1GS_A_GF3KF6.GT'])
        GS1R_AD.append(row['1GS_R_GF3KF6.AD'])
        GS1A_AD.append(row['1GS_A_GF3KF6.AD'])
print(i)


# GS2

# In[186]:


GS2R_AD = []
GS2A_AD = []
GS2R_GT = []
GS2A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['2GS_R_GF3KF6.GT']):
        GS2R_GT.append(row['2GS_A_GF3KF6.GT'])
        GS2A_GT.append(row['2GS_R_GF3KF6.GT'])
        GS2R_AD.append(row['2GS_A_GF3KF6.AD'])
        GS2A_AD.append(row['2GS_R_GF3KF6.AD'])
   
    else :
        GS2R_GT.append(row['2GS_R_GF3KF6.GT'])
        GS2A_GT.append(row['2GS_A_GF3KF6.GT'])
        GS2R_AD.append(row['2GS_R_GF3KF6.AD'])
        GS2A_AD.append(row['2GS_A_GF3KF6.AD'])
print(i)


# GS3

# In[187]:


GS3R_AD = []
GS3A_AD = []
GS3R_GT = []
GS3A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['3GS_R_GF3KF6.GT']):
        GS3R_GT.append(row['3GS_A_GF3KF6.GT'])
        GS3A_GT.append(row['3GS_R_GF3KF6.GT'])
        GS3R_AD.append(row['3GS_A_GF3KF6.AD'])
        GS3A_AD.append(row['3GS_R_GF3KF6.AD'])
   
    else :
        GS3R_GT.append(row['3GS_R_GF3KF6.GT'])
        GS3A_GT.append(row['3GS_A_GF3KF6.GT'])
        GS3R_AD.append(row['3GS_R_GF3KF6.AD'])
        GS3A_AD.append(row['3GS_A_GF3KF6.AD'])
print(i)


# GS4

# In[188]:


GS4R_AD = []
GS4A_AD = []
GS4R_GT = []
GS4A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['4GS_R_GF3KF6.GT']):
        GS4R_GT.append(row['4GS_A_GF3KF6.GT'])
        GS4A_GT.append(row['4GS_R_GF3KF6.GT'])
        GS4R_AD.append(row['4GS_A_GF3KF6.AD'])
        GS4A_AD.append(row['4GS_R_GF3KF6.AD'])
   
    else :
        GS4R_GT.append(row['4GS_R_GF3KF6.GT'])
        GS4A_GT.append(row['4GS_A_GF3KF6.GT'])
        GS4R_AD.append(row['4GS_R_GF3KF6.AD'])
        GS4A_AD.append(row['4GS_A_GF3KF6.AD'])


# GS5

# In[189]:


GS5R_AD = []
GS5A_AD = []
GS5R_GT = []
GS5A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['5GS_R_GF3KF6.GT']):
        GS5R_GT.append(row['5GS_A_GF3KF6.GT'])
        GS5A_GT.append(row['5GS_R_GF3KF6.GT'])
        GS5R_AD.append(row['5GS_A_GF3KF6.AD'])
        GS5A_AD.append(row['5GS_R_GF3KF6.AD'])
   
    else :
        GS5R_GT.append(row['5GS_R_GF3KF6.GT'])
        GS5A_GT.append(row['5GS_A_GF3KF6.GT'])
        GS5R_AD.append(row['5GS_R_GF3KF6.AD'])
        GS5A_AD.append(row['5GS_A_GF3KF6.AD'])
print(i)


# GS6

# In[190]:


GS6R_AD = []
GS6A_AD = []
GS6R_GT = []
GS6A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['6GS_R_GF3KF6.GT']):
        GS6R_GT.append(row['6GS_A_GF3KF6.GT'])
        GS6A_GT.append(row['6GS_R_GF3KF6.GT'])
        GS6R_AD.append(row['6GS_A_GF3KF6.AD'])
        GS6A_AD.append(row['6GS_R_GF3KF6.AD'])
   
    else :
        GS6R_GT.append(row['6GS_R_GF3KF6.GT'])
        GS6A_GT.append(row['6GS_A_GF3KF6.GT'])
        GS6R_AD.append(row['6GS_R_GF3KF6.AD'])
        GS6A_AD.append(row['6GS_A_GF3KF6.AD'])
print(i)


# KF1

# In[191]:


KF1R_AD = []
KF1A_AD = []
KF1R_GT = []
KF1A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['1KF_R_GF3KF6.GT']):
        KF1R_GT.append(row['1KF_A_GF3KF6.GT'])
        KF1A_GT.append(row['1KF_R_GF3KF6.GT'])
        KF1R_AD.append(row['1KF_A_GF3KF6.AD'])
        KF1A_AD.append(row['1KF_R_GF3KF6.AD'])
   
    else :
        KF1R_GT.append(row['1KF_R_GF3KF6.GT'])
        KF1A_GT.append(row['1KF_A_GF3KF6.GT'])
        KF1R_AD.append(row['1KF_R_GF3KF6.AD'])
        KF1A_AD.append(row['1KF_A_GF3KF6.AD'])
print(i)


# KF2

# In[192]:


KF2R_AD = []
KF2A_AD = []
KF2R_GT = []
KF2A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['2KF_R_GF3KF6.GT']):
        KF2R_GT.append(row['2KF_A_GF3KF6.GT'])
        KF2A_GT.append(row['2KF_R_GF3KF6.GT'])
        KF2R_AD.append(row['2KF_A_GF3KF6.AD'])
        KF2A_AD.append(row['2KF_R_GF3KF6.AD'])
   
    else :
        KF2R_GT.append(row['2KF_R_GF3KF6.GT'])
        KF2A_GT.append(row['2KF_A_GF3KF6.GT'])
        KF2R_AD.append(row['2KF_R_GF3KF6.AD'])
        KF2A_AD.append(row['2KF_A_GF3KF6.AD'])
print(i)


# KF3

# In[193]:


KF3R_AD = []
KF3A_AD = []
KF3R_GT = []
KF3A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['3KF_R_GF3KF6.GT']):
        KF3R_GT.append(row['3KF_A_GF3KF6.GT'])
        KF3A_GT.append(row['3KF_R_GF3KF6.GT'])
        KF3R_AD.append(row['3KF_A_GF3KF6.AD'])
        KF3A_AD.append(row['3KF_R_GF3KF6.AD'])
   
    else :
        KF3R_GT.append(row['3KF_R_GF3KF6.GT'])
        KF3A_GT.append(row['3KF_A_GF3KF6.GT'])
        KF3R_AD.append(row['3KF_R_GF3KF6.AD'])
        KF3A_AD.append(row['3KF_A_GF3KF6.AD'])
print(i)


# KF4

# In[194]:


KF4R_AD = []
KF4A_AD = []
KF4R_GT = []
KF4A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['4KF_R_GF3KF6.GT']):
        KF4R_GT.append(row['4KF_A_GF3KF6.GT'])
        KF4A_GT.append(row['4KF_R_GF3KF6.GT'])
        KF4R_AD.append(row['4KF_A_GF3KF6.AD'])
        KF4A_AD.append(row['4KF_R_GF3KF6.AD'])
   
    else :
        KF4R_GT.append(row['4KF_R_GF3KF6.GT'])
        KF4A_GT.append(row['4KF_A_GF3KF6.GT'])
        KF4R_AD.append(row['4KF_R_GF3KF6.AD'])
        KF4A_AD.append(row['4KF_A_GF3KF6.AD'])
print(i)


# KF5

# In[195]:


KF5R_AD = []
KF5A_AD = []
KF5R_GT = []
KF5A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['5KF_R_GF3KF6.GT']):
        KF5R_GT.append(row['5KF_A_GF3KF6.GT'])
        KF5A_GT.append(row['5KF_R_GF3KF6.GT'])
        KF5R_AD.append(row['5KF_A_GF3KF6.AD'])
        KF5A_AD.append(row['5KF_R_GF3KF6.AD'])
   
    else :
        KF5R_GT.append(row['5KF_R_GF3KF6.GT'])
        KF5A_GT.append(row['5KF_A_GF3KF6.GT'])
        KF5R_AD.append(row['5KF_R_GF3KF6.AD'])
        KF5A_AD.append(row['5KF_A_GF3KF6.AD'])
print(i)


# KF6

# In[196]:


KF6R_AD = []
KF6A_AD = []
KF6R_GT = []
KF6A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['6KF_R_GF3KF6.GT']):
        KF6R_GT.append(row['6KF_A_GF3KF6.GT'])
        KF6A_GT.append(row['6KF_R_GF3KF6.GT'])
        KF6R_AD.append(row['6KF_A_GF3KF6.AD'])
        KF6A_AD.append(row['6KF_R_GF3KF6.AD'])
   
    else :
        KF6R_GT.append(row['6KF_R_GF3KF6.GT'])
        KF6A_GT.append(row['6KF_A_GF3KF6.GT'])
        KF6R_AD.append(row['6KF_R_GF3KF6.AD'])
        KF6A_AD.append(row['6KF_A_GF3KF6.AD'])
print(i)


# KS2

# In[197]:


KS2R_AD = []
KS2A_AD = []
KS2R_GT = []
KS2A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['2KS_R_GF3KF6.GT']):
        KS2R_GT.append(row['2KS_A_GF3KF6.GT'])
        KS2A_GT.append(row['2KS_R_GF3KF6.GT'])
        KS2R_AD.append(row['2KS_A_GF3KF6.AD'])
        KS2A_AD.append(row['2KS_R_GF3KF6.AD'])
   
    else :
        KS2R_GT.append(row['2KS_R_GF3KF6.GT'])
        KS2A_GT.append(row['2KS_A_GF3KF6.GT'])
        KS2R_AD.append(row['2KS_R_GF3KF6.AD'])
        KS2A_AD.append(row['2KS_A_GF3KF6.AD'])
print(i)


# KS3

# In[198]:


KS3R_AD = []
KS3A_AD = []
KS3R_GT = []
KS3A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['3KS_R_GF3KF6.GT']):
        KS3R_GT.append(row['3KS_A_GF3KF6.GT'])
        KS3A_GT.append(row['3KS_R_GF3KF6.GT'])
        KS3R_AD.append(row['3KS_A_GF3KF6.AD'])
        KS3A_AD.append(row['3KS_R_GF3KF6.AD'])
   
    else :
        KS3R_GT.append(row['3KS_R_GF3KF6.GT'])
        KS3A_GT.append(row['3KS_A_GF3KF6.GT'])
        KS3R_AD.append(row['3KS_R_GF3KF6.AD'])
        KS3A_AD.append(row['3KS_A_GF3KF6.AD'])


# KS4

# In[199]:


KS4R_AD = []
KS4A_AD = []
KS4R_GT = []
KS4A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['4KS_R_GF3KF6.GT']):
        KS4R_GT.append(row['4KS_A_GF3KF6.GT'])
        KS4A_GT.append(row['4KS_R_GF3KF6.GT'])
        KS4R_AD.append(row['4KS_A_GF3KF6.AD'])
        KS4A_AD.append(row['4KS_R_GF3KF6.AD'])
   
    else :
        KS4R_GT.append(row['4KS_R_GF3KF6.GT'])
        KS4A_GT.append(row['4KS_A_GF3KF6.GT'])
        KS4R_AD.append(row['4KS_R_GF3KF6.AD'])
        KS4A_AD.append(row['4KS_A_GF3KF6.AD'])


# KS5

# In[200]:


KS5R_AD = []
KS5A_AD = []
KS5R_GT = []
KS5A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['5KS_R_GF3KF6.GT']):
        KS5R_GT.append(row['5KS_A_GF3KF6.GT'])
        KS5A_GT.append(row['5KS_R_GF3KF6.GT'])
        KS5R_AD.append(row['5KS_A_GF3KF6.AD'])
        KS5A_AD.append(row['5KS_R_GF3KF6.AD'])
   
    else :
        KS5R_GT.append(row['5KS_R_GF3KF6.GT'])
        KS5A_GT.append(row['5KS_A_GF3KF6.GT'])
        KS5R_AD.append(row['5KS_R_GF3KF6.AD'])
        KS5A_AD.append(row['5KS_A_GF3KF6.AD'])


# KS6

# In[201]:


KS6R_AD = []
KS6A_AD = []
KS6R_GT = []
KS6A_GT = []
for index, row in df10.iterrows() : 
    if (row['6KF_R_GF3KF6.GT'] != row['6KS_R_GF3KF6.GT']):
        KS6R_GT.append(row['6KS_A_GF3KF6.GT'])
        KS6A_GT.append(row['6KS_R_GF3KF6.GT'])
        KS6R_AD.append(row['6KS_A_GF3KF6.AD'])
        KS6A_AD.append(row['6KS_R_GF3KF6.AD'])
   
    else :
        KS6R_GT.append(row['6KS_R_GF3KF6.GT'])
        KS6A_GT.append(row['6KS_A_GF3KF6.GT'])
        KS6R_AD.append(row['6KS_R_GF3KF6.AD'])
        KS6A_AD.append(row['6KS_A_GF3KF6.AD'])


# Set the different chains to the table column. Make a backup first.

# In[202]:


df_10_backup = df10 

# For quality control
df10['1GF_R.GT'] = df10['1GF_R_final.GT']
df10['2GF_R.GT'] = GF2R_GT 
df10['3GF_R.GT'] = GF3R_GT 
df10['4GF_R.GT'] = GF4R_GT 
df10['5GF_R.GT'] = GF5R_GT 
df10['6GF_R.GT'] = GF6R_GT 
df10['1GS_R.GT'] = GS1R_GT 
df10['2GS_R.GT'] = GS2R_GT 
df10['3GS_R.GT'] = GS3R_GT 
df10['4GS_R.GT'] = GS4R_GT 
df10['5GS_R.GT'] = GS5R_GT 
df10['6GS_R.GT'] = GS6R_GT 
df10['1KF_R.GT'] = KF1R_GT 
df10['2KF_R.GT'] = KF2R_GT 
df10['3KF_R.GT'] = KF3R_GT 
df10['4KF_R.GT'] = KF4R_GT 
df10['5KF_R.GT'] = KF5R_GT 
df10['6KF_R.GT'] = KF6R_GT 
df10['2KS_R.GT'] = KS2R_GT 
df10['3KS_R.GT'] = KS3R_GT 
df10['4KS_R.GT'] = KS4R_GT 
df10['5KS_R.GT'] = KS5R_GT 
df10['6KS_R.GT'] = KS6R_GT 

df10['1GF_A.GT'] = df10['1GF_A_final.GT']
df10['2GF_A.GT'] = GF2A_GT 
df10['3GF_A.GT'] = GF3A_GT 
df10['4GF_A.GT'] = GF4A_GT 
df10['5GF_A.GT'] = GF5A_GT 
df10['6GF_A.GT'] = GF6A_GT 
df10['1GS_A.GT'] = GS1A_GT 
df10['2GS_A.GT'] = GS2A_GT 
df10['3GS_A.GT'] = GS3A_GT 
df10['4GS_A.GT'] = GS4A_GT 
df10['5GS_A.GT'] = GS5A_GT 
df10['6GS_A.GT'] = GS6A_GT 
df10['1KF_A.GT'] = KF1A_GT 
df10['2KF_A.GT'] = KF2A_GT 
df10['3KF_A.GT'] = KF3A_GT 
df10['4KF_A.GT'] = KF4A_GT 
df10['5KF_A.GT'] = KF5A_GT 
df10['6KF_A.GT'] = KF6A_GT 
df10['2KS_A.GT'] = KS2A_GT 
df10['3KS_A.GT'] = KS3A_GT 
df10['4KS_A.GT'] = KS4A_GT 
df10['5KS_A.GT'] = KS5A_GT 
df10['6KS_A.GT'] = KS6A_GT 

df10['1GF_R.AD'] = df10['1GF_R_final.AD']
df10['2GF_R.AD'] = GF2R_AD 
df10['3GF_R.AD'] = GF3R_AD 
df10['4GF_R.AD'] = GF4R_AD 
df10['5GF_R.AD'] = GF5R_AD 
df10['6GF_R.AD'] = GF6R_AD 
df10['1GS_R.AD'] = GS1R_AD 
df10['2GS_R.AD'] = GS2R_AD 
df10['3GS_R.AD'] = GS3R_AD 
df10['4GS_R.AD'] = GS4R_AD 
df10['5GS_R.AD'] = GS5R_AD 
df10['6GS_R.AD'] = GS6R_AD 
df10['1KF_R.AD'] = KF1R_AD 
df10['2KF_R.AD'] = KF2R_AD 
df10['3KF_R.AD'] = KF3R_AD 
df10['4KF_R.AD'] = KF4R_AD 
df10['5KF_R.AD'] = KF5R_AD 
df10['6KF_R.AD'] = KF6R_AD 
df10['2KS_R.AD'] = KS2R_AD 
df10['3KS_R.AD'] = KS3R_AD 
df10['4KS_R.AD'] = KS4R_AD 
df10['5KS_R.AD'] = KS5R_AD 
df10['6KS_R.AD'] = KS6R_AD 

df10['1GF_A.AD'] = df10['1GF_A_final.AD']
df10['2GF_A.AD'] = GF2A_AD 
df10['3GF_A.AD'] = GF3A_AD 
df10['4GF_A.AD'] = GF4A_AD 
df10['5GF_A.AD'] = GF5A_AD 
df10['6GF_A.AD'] = GF6A_AD 
df10['1GS_A.AD'] = GS1A_AD 
df10['2GS_A.AD'] = GS2A_AD 
df10['3GS_A.AD'] = GS3A_AD 
df10['4GS_A.AD'] = GS4A_AD 
df10['5GS_A.AD'] = GS5A_AD 
df10['6GS_A.AD'] = GS6A_AD 
df10['1KF_A.AD'] = KF1A_AD 
df10['2KF_A.AD'] = KF2A_AD 
df10['3KF_A.AD'] = KF3A_AD 
df10['4KF_A.AD'] = KF4A_AD 
df10['5KF_A.AD'] = KF5A_AD 
df10['6KF_A.AD'] = KF6A_AD 
df10['2KS_A.AD'] = KS2A_AD 
df10['3KS_A.AD'] = KS3A_AD 
df10['4KS_A.AD'] = KS4A_AD 
df10['5KS_A.AD'] = KS5A_AD 
df10['6KS_A.AD'] = KS6A_AD 
# In[203]:


df10 = df_10_backup


# In[204]:


df10['1GF_R_GF3KF6.GT'] = GF1R_GT 
df10['2GF_R_GF3KF6.GT'] = GF2R_GT 
df10['3GF_R_GF3KF6.GT'] = GF3R_GT 
df10['4GF_R_GF3KF6.GT'] = GF4R_GT 
df10['5GF_R_GF3KF6.GT'] = GF5R_GT 
df10['6GF_R_GF3KF6.GT'] = GF6R_GT 
df10['1GS_R_GF3KF6.GT'] = GS1R_GT 
df10['2GS_R_GF3KF6.GT'] = GS2R_GT 
df10['3GS_R_GF3KF6.GT'] = GS3R_GT 
df10['4GS_R_GF3KF6.GT'] = GS4R_GT 
df10['5GS_R_GF3KF6.GT'] = GS5R_GT 
df10['6GS_R_GF3KF6.GT'] = GS6R_GT 
df10['1KF_R_GF3KF6.GT'] = KF1R_GT 
df10['2KF_R_GF3KF6.GT'] = KF2R_GT 
df10['3KF_R_GF3KF6.GT'] = KF3R_GT 
df10['4KF_R_GF3KF6.GT'] = KF4R_GT 
df10['5KF_R_GF3KF6.GT'] = KF5R_GT 
df10['6KF_R_GF3KF6.GT'] = KF6R_GT 
df10['2KS_R_GF3KF6.GT'] = KS2R_GT 
df10['3KS_R_GF3KF6.GT'] = KS3R_GT 
df10['4KS_R_GF3KF6.GT'] = KS4R_GT 
df10['5KS_R_GF3KF6.GT'] = KS5R_GT 
df10['6KS_R_GF3KF6.GT'] = KS6R_GT 

df10['1GF_A_GF3KF6.GT'] = GF1A_GT 
df10['2GF_A_GF3KF6.GT'] = GF2A_GT 
df10['3GF_A_GF3KF6.GT'] = GF3A_GT 
df10['4GF_A_GF3KF6.GT'] = GF4A_GT 
df10['5GF_A_GF3KF6.GT'] = GF5A_GT 
df10['6GF_A_GF3KF6.GT'] = GF6A_GT 
df10['1GS_A_GF3KF6.GT'] = GS1A_GT 
df10['2GS_A_GF3KF6.GT'] = GS2A_GT 
df10['3GS_A_GF3KF6.GT'] = GS3A_GT 
df10['4GS_A_GF3KF6.GT'] = GS4A_GT 
df10['5GS_A_GF3KF6.GT'] = GS5A_GT 
df10['6GS_A_GF3KF6.GT'] = GS6A_GT 
df10['1KF_A_GF3KF6.GT'] = KF1A_GT 
df10['2KF_A_GF3KF6.GT'] = KF2A_GT 
df10['3KF_A_GF3KF6.GT'] = KF3A_GT 
df10['4KF_A_GF3KF6.GT'] = KF4A_GT 
df10['5KF_A_GF3KF6.GT'] = KF5A_GT 
df10['6KF_A_GF3KF6.GT'] = KF6A_GT 
df10['2KS_A_GF3KF6.GT'] = KS2A_GT 
df10['3KS_A_GF3KF6.GT'] = KS3A_GT 
df10['4KS_A_GF3KF6.GT'] = KS4A_GT 
df10['5KS_A_GF3KF6.GT'] = KS5A_GT 
df10['6KS_A_GF3KF6.GT'] = KS6A_GT 

df10['1GF_R_GF3KF6.AD'] = GF1R_AD 
df10['2GF_R_GF3KF6.AD'] = GF2R_AD 
df10['3GF_R_GF3KF6.AD'] = GF3R_AD 
df10['4GF_R_GF3KF6.AD'] = GF4R_AD 
df10['5GF_R_GF3KF6.AD'] = GF5R_AD 
df10['6GF_R_GF3KF6.AD'] = GF6R_AD 
df10['1GS_R_GF3KF6.AD'] = GS1R_AD 
df10['2GS_R_GF3KF6.AD'] = GS2R_AD 
df10['3GS_R_GF3KF6.AD'] = GS3R_AD 
df10['4GS_R_GF3KF6.AD'] = GS4R_AD 
df10['5GS_R_GF3KF6.AD'] = GS5R_AD 
df10['6GS_R_GF3KF6.AD'] = GS6R_AD 
df10['1KF_R_GF3KF6.AD'] = KF1R_AD 
df10['2KF_R_GF3KF6.AD'] = KF2R_AD 
df10['3KF_R_GF3KF6.AD'] = KF3R_AD 
df10['4KF_R_GF3KF6.AD'] = KF4R_AD 
df10['5KF_R_GF3KF6.AD'] = KF5R_AD 
df10['6KF_R_GF3KF6.AD'] = KF6R_AD 
df10['2KS_R_GF3KF6.AD'] = KS2R_AD 
df10['3KS_R_GF3KF6.AD'] = KS3R_AD 
df10['4KS_R_GF3KF6.AD'] = KS4R_AD 
df10['5KS_R_GF3KF6.AD'] = KS5R_AD 
df10['6KS_R_GF3KF6.AD'] = KS6R_AD 

df10['1GF_A_GF3KF6.AD'] = GF1A_AD 
df10['2GF_A_GF3KF6.AD'] = GF2A_AD 
df10['3GF_A_GF3KF6.AD'] = GF3A_AD 
df10['4GF_A_GF3KF6.AD'] = GF4A_AD 
df10['5GF_A_GF3KF6.AD'] = GF5A_AD 
df10['6GF_A_GF3KF6.AD'] = GF6A_AD 
df10['1GS_A_GF3KF6.AD'] = GS1A_AD 
df10['2GS_A_GF3KF6.AD'] = GS2A_AD 
df10['3GS_A_GF3KF6.AD'] = GS3A_AD 
df10['4GS_A_GF3KF6.AD'] = GS4A_AD 
df10['5GS_A_GF3KF6.AD'] = GS5A_AD 
df10['6GS_A_GF3KF6.AD'] = GS6A_AD 
df10['1KF_A_GF3KF6.AD'] = KF1A_AD 
df10['2KF_A_GF3KF6.AD'] = KF2A_AD 
df10['3KF_A_GF3KF6.AD'] = KF3A_AD 
df10['4KF_A_GF3KF6.AD'] = KF4A_AD 
df10['5KF_A_GF3KF6.AD'] = KF5A_AD 
df10['6KF_A_GF3KF6.AD'] = KF6A_AD 
df10['2KS_A_GF3KF6.AD'] = KS2A_AD 
df10['3KS_A_GF3KF6.AD'] = KS3A_AD 
df10['4KS_A_GF3KF6.AD'] = KS4A_AD 
df10['5KS_A_GF3KF6.AD'] = KS5A_AD 
df10['6KS_A_GF3KF6.AD'] = KS6A_AD 


# In[205]:


df10.head(20)


# In[206]:


#df10.to_csv('Quality control refalt asignation.csv')


# In[207]:


len(df10)


# # Clean up the monoallelic expression data

# The monoallelic expression (MAE) can result from homozygots as well as imprinted genes. In order to distinguish each case, we need to know the genotype. However, that's not possible for the amount of SNPs that we are working with. Therefore we drop all the SNPs whose allele frequency is 0 or 1. 
# 
# Calculate the allele frequencies first.

# In[208]:


df10['AF_1GF'] = (df10['1GF_R_GF3KF6.AD']/(df10['1GF_R_GF3KF6.AD'] + df10['1GF_A_GF3KF6.AD']))
df10['AF_2GF'] = (df10['2GF_R_GF3KF6.AD']/(df10['2GF_R_GF3KF6.AD'] + df10['2GF_A_GF3KF6.AD']))
df10['AF_3GF'] = (df10['3GF_R_GF3KF6.AD']/(df10['3GF_R_GF3KF6.AD'] + df10['3GF_A_GF3KF6.AD']))
df10['AF_4GF'] = (df10['4GF_R_GF3KF6.AD']/(df10['4GF_R_GF3KF6.AD'] + df10['4GF_A_GF3KF6.AD']))
df10['AF_5GF'] = (df10['5GF_R_GF3KF6.AD']/(df10['5GF_R_GF3KF6.AD'] + df10['5GF_A_GF3KF6.AD']))
df10['AF_6GF'] = (df10['6GF_R_GF3KF6.AD']/(df10['6GF_R_GF3KF6.AD'] + df10['6GF_A_GF3KF6.AD']))

df10['AF_1GS'] = (df10['1GS_R_GF3KF6.AD']/(df10['1GS_R_GF3KF6.AD'] + df10['1GS_A_GF3KF6.AD']))
df10['AF_2GS'] = (df10['2GS_R_GF3KF6.AD']/(df10['2GS_R_GF3KF6.AD'] + df10['2GS_A_GF3KF6.AD']))
df10['AF_3GS'] = (df10['3GS_R_GF3KF6.AD']/(df10['3GS_R_GF3KF6.AD'] + df10['3GS_A_GF3KF6.AD']))
df10['AF_4GS'] = (df10['4GS_R_GF3KF6.AD']/(df10['4GS_R_GF3KF6.AD'] + df10['4GS_A_GF3KF6.AD']))
df10['AF_5GS'] = (df10['5GS_R_GF3KF6.AD']/(df10['5GS_R_GF3KF6.AD'] + df10['5GS_A_GF3KF6.AD']))
df10['AF_6GS'] = (df10['6GS_R_GF3KF6.AD']/(df10['6GS_R_GF3KF6.AD'] + df10['6GS_A_GF3KF6.AD']))

df10['AF_1KF'] = (df10['1KF_R_GF3KF6.AD']/(df10['1KF_R_GF3KF6.AD'] + df10['1KF_A_GF3KF6.AD']))
df10['AF_2KF'] = (df10['2KF_R_GF3KF6.AD']/(df10['2KF_R_GF3KF6.AD'] + df10['2KF_A_GF3KF6.AD']))
df10['AF_3KF'] = (df10['3KF_R_GF3KF6.AD']/(df10['3KF_R_GF3KF6.AD'] + df10['3KF_A_GF3KF6.AD']))
df10['AF_4KF'] = (df10['4KF_R_GF3KF6.AD']/(df10['4KF_R_GF3KF6.AD'] + df10['4KF_A_GF3KF6.AD']))
df10['AF_5KF'] = (df10['5KF_R_GF3KF6.AD']/(df10['5KF_R_GF3KF6.AD'] + df10['5KF_A_GF3KF6.AD']))
df10['AF_6KF'] = (df10['6KF_R_GF3KF6.AD']/(df10['6KF_R_GF3KF6.AD'] + df10['6KF_A_GF3KF6.AD']))

df10['AF_2KS'] = (df10['2KS_R_GF3KF6.AD']/(df10['2KS_R_GF3KF6.AD'] + df10['2KS_A_GF3KF6.AD']))
df10['AF_3KS'] = (df10['3KS_R_GF3KF6.AD']/(df10['3KS_R_GF3KF6.AD'] + df10['3KS_A_GF3KF6.AD']))
df10['AF_4KS'] = (df10['4KS_R_GF3KF6.AD']/(df10['4KS_R_GF3KF6.AD'] + df10['4KS_A_GF3KF6.AD']))
df10['AF_5KS'] = (df10['5KS_R_GF3KF6.AD']/(df10['5KS_R_GF3KF6.AD'] + df10['5KS_A_GF3KF6.AD']))
df10['AF_6KS'] = (df10['6KS_R_GF3KF6.AD']/(df10['6KS_R_GF3KF6.AD'] + df10['6KS_A_GF3KF6.AD']))


# In[209]:


df10 = df10.fillna(0)


# In[210]:


len(df10)


# Collect the indexes whose allele frequencies are 0 or 1 in both tissues of the same sample. These indexes correspond to the SNPs that are MAE for at least one sample.

# In[211]:


index_AF_1GF = df10[((df10['AF_1GF'] == 0) | (df10['AF_1GF'] == 1)) & ((df10['AF_1KF'] == 0) | (df10['AF_1KF'] == 1))].index
index_AF_2GF = df10[((df10['AF_2GF'] == 0) | (df10['AF_2GF'] == 1)) & ((df10['AF_2KF'] == 0) | (df10['AF_2KF'] == 1))].index
index_AF_3GF = df10[((df10['AF_3GF'] == 0) | (df10['AF_3GF'] == 1)) & ((df10['AF_3KF'] == 0) | (df10['AF_3KF'] == 1))].index
index_AF_4GF = df10[((df10['AF_4GF'] == 0) | (df10['AF_4GF'] == 1)) & ((df10['AF_4KF'] == 0) | (df10['AF_4KF'] == 1))].index
index_AF_5GF = df10[((df10['AF_5GF'] == 0) | (df10['AF_5GF'] == 1)) & ((df10['AF_5KF'] == 0) | (df10['AF_5KF'] == 1))].index
index_AF_6GF = df10[((df10['AF_6GF'] == 0) | (df10['AF_6GF'] == 1)) & ((df10['AF_6KF'] == 0) | (df10['AF_6KF'] == 1))].index

index_AF_2GS = df10[((df10['AF_2GS'] == 0) | (df10['AF_2GS'] == 1)) & ((df10['AF_2KS'] == 0) | (df10['AF_2KS'] == 1))].index
index_AF_3GS = df10[((df10['AF_3GS'] == 0) | (df10['AF_3GS'] == 1)) & ((df10['AF_3KS'] == 0) | (df10['AF_3KS'] == 1))].index
index_AF_4GS = df10[((df10['AF_4GS'] == 0) | (df10['AF_4GS'] == 1)) & ((df10['AF_4KS'] == 0) | (df10['AF_4KS'] == 1))].index
index_AF_5GS = df10[((df10['AF_5GS'] == 0) | (df10['AF_5GS'] == 1)) & ((df10['AF_5KS'] == 0) | (df10['AF_5KS'] == 1))].index
index_AF_6GS = df10[((df10['AF_6GS'] == 0) | (df10['AF_6GS'] == 1)) & ((df10['AF_6KS'] == 0) | (df10['AF_6KS'] == 1))].index

#index_AF_1KF = df10[((df10['AF_1KF'] == 0) | (df10['AF_1KF'] == 1) & (df10['AF_1GF'] != 0) | (df10['AF_1GF'] != 1))].index
#index_AF_2KF = df10[((df10['AF_2KF'] == 0) | (df10['AF_2KF'] == 1) & (df10['AF_2GF'] != 0) | (df10['AF_2GF'] != 1))].index
#index_AF_3KF = df10[((df10['AF_3KF'] == 0) | (df10['AF_3KF'] == 1) & (df10['AF_3GF'] != 0) | (df10['AF_3GF'] != 1))].index
#index_AF_4KF = df10[((df10['AF_4KF'] == 0) | (df10['AF_4KF'] == 1) & (df10['AF_4GF'] != 0) | (df10['AF_4GF'] != 1))].index
#index_AF_5KF = df10[((df10['AF_5KF'] == 0) | (df10['AF_5KF'] == 1) & (df10['AF_5GF'] != 0) | (df10['AF_5GF'] != 1))].index
#index_AF_6KF = df10[((df10['AF_6KF'] == 0) | (df10['AF_6KF'] == 1) & (df10['AF_6GF'] != 0) | (df10['AF_6GF'] != 1))].index

#index_AF_2KS = df10[((df10['AF_2KS'] == 0) | (df10['AF_2KS'] == 1) & (df10['AF_2GS'] != 0) | (df10['AF_2GS'] != 1))].index
#index_AF_3KS = df10[((df10['AF_3KS'] == 0) | (df10['AF_3KS'] == 1) & (df10['AF_3GS'] != 0) | (df10['AF_3GS'] != 1))].index
#index_AF_4KS = df10[((df10['AF_4KS'] == 0) | (df10['AF_4KS'] == 1) & (df10['AF_4GS'] != 0) | (df10['AF_4GS'] != 1))].index
#index_AF_5KS = df10[((df10['AF_5KS'] == 0) | (df10['AF_5KS'] == 1) & (df10['AF_5GS'] != 0) | (df10['AF_5GS'] != 1))].index
#index_AF_6KS = df10[((df10['AF_6KS'] == 0) | (df10['AF_6KS'] == 1) & (df10['AF_6GS'] != 0) | (df10['AF_6GS'] != 1))].index


# Keep all the indexes that are not MAE for the SNPs in the dference genome and make a new matrix 

# In[212]:


MAE_df = index_AF_1GF.union(index_AF_2GF)
MAE_df = MAE_df.union(index_AF_2GS)
#MAE_df = MAE_df.union(index_AF_2KF)
#MAE_df = MAE_df.union(index_AF_2KS)
MAE_df = MAE_df.union(index_AF_3GF)
MAE_df = MAE_df.union(index_AF_3GS)
#MAE_df = MAE_df.union(index_AF_3KF)
#MAE_df = MAE_df.union(index_AF_3KS)
MAE_df = MAE_df.union(index_AF_4GF)
MAE_df = MAE_df.union(index_AF_4GS)
#MAE_df = MAE_df.union(index_AF_4KF)
#MAE_df = MAE_df.union(index_AF_4KS)
MAE_df = MAE_df.union(index_AF_5GF)
MAE_df = MAE_df.union(index_AF_5GS)
#MAE_df = MAE_df.union(index_AF_5KF)
#MAE_df = MAE_df.union(index_AF_5KS)
MAE_df = MAE_df.union(index_AF_6GF)
MAE_df = MAE_df.union(index_AF_6GS)
#MAE_df = MAE_df.union(index_AF_6KF)
#MAE_df = MAE_df.union(index_AF_6KS)
#MAE_df = MAE_df.union(index_AF1)
#MAE_df = MAE_df.union(index_AF0)


# Make the new DFs with the non MAE SNPs

# In[213]:


MAE = df10.drop(index = MAE_df) #MAE = df10.loc[MAE_df,:]


# In[214]:


len(MAE)


# In[215]:


MAE.head(10)


# In[216]:


len(MAE)


# In[ ]:


MAE.drop(['Unnamed: 0', 'Unnamed: 0.1'], axis=1)


# In[217]:


MAE.to_csv('SNPs_2map_MAE_GF3KF6.csv')


# Finished! Now let's treat the data with the script "Manhattan allele distance plot" for controlling the mapping bias and doing the MAE analysis.

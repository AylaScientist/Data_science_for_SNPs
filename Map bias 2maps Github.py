#!/usr/bin/env python
# coding: utf-8

# # Analysis of the mapping bias

# This script will analyse the mapping bias on the file that was submitted to Data wrangling. It is necessary calculate the allele frequencies for this analysis

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


df=pd.read_csv('SNPs_2map_MAE_GF3KF6.csv')


# In[3]:


df.head()


# In[4]:


# get the column names
for col in df.columns: 
    print(col) 


# Calculate allele frequencies based on the reference allele over the total allele counts for the same SNP.

# In[5]:


df['AF']=((df['1GF_R_GF3KF6.AD']+
           df['2GF_R_GF3KF6.AD']+
           df['3GF_R_GF3KF6.AD']+
           df['4GF_R_GF3KF6.AD']+
           df['5GF_R_GF3KF6.AD']+
           df['6GF_R_GF3KF6.AD']+
           df['1GS_R_GF3KF6.AD']+
           df['2GS_R_GF3KF6.AD']+
           df['3GS_R_GF3KF6.AD']+
           df['4GS_R_GF3KF6.AD']+
           df['5GS_R_GF3KF6.AD']+
           df['6GS_R_GF3KF6.AD']+
           df['1KF_R_GF3KF6.AD']+
           df['2KF_R_GF3KF6.AD']+
           df['3KF_R_GF3KF6.AD']+
           df['4KF_R_GF3KF6.AD']+
           df['5KF_R_GF3KF6.AD']+
           df['6KF_R_GF3KF6.AD']+
           df['2KS_R_GF3KF6.AD']+
           df['3KS_R_GF3KF6.AD']+
           df['4KS_R_GF3KF6.AD']+
           df['5KS_R_GF3KF6.AD']+
           df['6KS_R_GF3KF6.AD'])/
          (df['1GF_R_GF3KF6.AD']+
           df['2GF_R_GF3KF6.AD']+
           df['3GF_R_GF3KF6.AD']+
           df['4GF_R_GF3KF6.AD']+
           df['5GF_R_GF3KF6.AD']+
           df['6GF_R_GF3KF6.AD']+
           df['1GS_R_GF3KF6.AD']+
           df['2GS_R_GF3KF6.AD']+
           df['3GS_R_GF3KF6.AD']+
           df['4GS_R_GF3KF6.AD']+
           df['5GS_R_GF3KF6.AD']+
           df['6GS_R_GF3KF6.AD']+
           df['1KF_R_GF3KF6.AD']+
           df['2KF_R_GF3KF6.AD']+
           df['3KF_R_GF3KF6.AD']+
           df['4KF_R_GF3KF6.AD']+
           df['5KF_R_GF3KF6.AD']+
           df['6KF_R_GF3KF6.AD']+
           df['2KS_R_GF3KF6.AD']+
           df['3KS_R_GF3KF6.AD']+
           df['4KS_R_GF3KF6.AD']+
           df['5KS_R_GF3KF6.AD']+
           df['6KS_R_GF3KF6.AD']+
           df['1GF_A_GF3KF6.AD']+
           df['2GF_A_GF3KF6.AD']+
           df['3GF_A_GF3KF6.AD']+
           df['4GF_A_GF3KF6.AD']+
           df['5GF_A_GF3KF6.AD']+
           df['6GF_A_GF3KF6.AD']+
           df['1GS_A_GF3KF6.AD']+
           df['2GS_A_GF3KF6.AD']+
           df['3GS_A_GF3KF6.AD']+
           df['4GS_A_GF3KF6.AD']+
           df['5GS_A_GF3KF6.AD']+
           df['6GS_A_GF3KF6.AD']+
           df['1KF_A_GF3KF6.AD']+
           df['2KF_A_GF3KF6.AD']+
           df['3KF_A_GF3KF6.AD']+
           df['4KF_A_GF3KF6.AD']+
           df['5KF_A_GF3KF6.AD']+
           df['6KF_A_GF3KF6.AD']+
           df['2KS_A_GF3KF6.AD']+
           df['3KS_A_GF3KF6.AD']+
           df['4KS_A_GF3KF6.AD']+
           df['5KS_A_GF3KF6.AD']+
           df['6KS_A_GF3KF6.AD']))


# In[6]:


df['AF_GF']=((df['1GF_R_GF3KF6.AD']+
           df['2GF_R_GF3KF6.AD']+
           df['3GF_R_GF3KF6.AD']+
           df['4GF_R_GF3KF6.AD']+
           df['5GF_R_GF3KF6.AD']+
           df['6GF_R_GF3KF6.AD'])/
          (df['1GF_R_GF3KF6.AD']+
           df['2GF_R_GF3KF6.AD']+
           df['3GF_R_GF3KF6.AD']+
           df['4GF_R_GF3KF6.AD']+
           df['5GF_R_GF3KF6.AD']+
           df['6GF_R_GF3KF6.AD']+
           df['1GF_A_GF3KF6.AD']+
           df['2GF_A_GF3KF6.AD']+
           df['3GF_A_GF3KF6.AD']+
           df['4GF_A_GF3KF6.AD']+
           df['5GF_A_GF3KF6.AD']+
           df['6GF_A_GF3KF6.AD']))


# In[7]:


df['AF_GS']=((df['1GS_R_GF3KF6.AD']+
           df['2GS_R_GF3KF6.AD']+
           df['3GS_R_GF3KF6.AD']+
           df['4GS_R_GF3KF6.AD']+
           df['5GS_R_GF3KF6.AD']+
           df['6GS_R_GF3KF6.AD'])/
          (df['1GS_R_GF3KF6.AD']+
           df['2GS_R_GF3KF6.AD']+
           df['3GS_R_GF3KF6.AD']+
           df['4GS_R_GF3KF6.AD']+
           df['5GS_R_GF3KF6.AD']+
           df['6GS_R_GF3KF6.AD']+
           df['1GS_A_GF3KF6.AD']+
           df['2GS_A_GF3KF6.AD']+
           df['3GS_A_GF3KF6.AD']+
           df['4GS_A_GF3KF6.AD']+
           df['5GS_A_GF3KF6.AD']+
           df['6GS_A_GF3KF6.AD']))


# In[8]:


df['AF_KF']=((df['1KF_R_GF3KF6.AD']+
           df['2KF_R_GF3KF6.AD']+
           df['3KF_R_GF3KF6.AD']+
           df['4KF_R_GF3KF6.AD']+
           df['5KF_R_GF3KF6.AD']+
           df['6KF_R_GF3KF6.AD'])/
          (df['1KF_R_GF3KF6.AD']+
           df['2KF_R_GF3KF6.AD']+
           df['3KF_R_GF3KF6.AD']+
           df['4KF_R_GF3KF6.AD']+
           df['5KF_R_GF3KF6.AD']+
           df['6KF_R_GF3KF6.AD']+
           df['1KF_A_GF3KF6.AD']+
           df['2KF_A_GF3KF6.AD']+
           df['3KF_A_GF3KF6.AD']+
           df['4KF_A_GF3KF6.AD']+
           df['5KF_A_GF3KF6.AD']+
           df['6KF_A_GF3KF6.AD']))


# In[9]:


df['AF_KS']=((df['2KS_R_GF3KF6.AD']+
           df['3KS_R_GF3KF6.AD']+
           df['4KS_R_GF3KF6.AD']+
           df['5KS_R_GF3KF6.AD']+
           df['6KS_R_GF3KF6.AD'])/
          (df['2KS_R_GF3KF6.AD']+
           df['3KS_R_GF3KF6.AD']+
           df['4KS_R_GF3KF6.AD']+
           df['5KS_R_GF3KF6.AD']+
           df['6KS_R_GF3KF6.AD']+
           df['2KS_A_GF3KF6.AD']+
           df['3KS_A_GF3KF6.AD']+
           df['4KS_A_GF3KF6.AD']+
           df['5KS_A_GF3KF6.AD']+
           df['6KS_A_GF3KF6.AD']))


# In[10]:


df.to_csv('SNPs2maps_GF3KF6.csv')


# This file can be submitted to the script "Statistical tests"

# ## Graphical analysis of the mapping bias

# This is a istogram of the allele frequencie values. It should follow a normal distribution

# In[74]:


hist_data = df['AF'].values
# the histogram of the data
plt.hist(hist_data, 50, density=False, facecolor='green')

plt.xlabel('Allele frequency')
plt.ylabel('Frequency')
plt.title('Histogram Example')

plt.grid(True)

plt.show()


# In[75]:


hist_data = df['AF_GF'].values
# the histogram of the data
plt.hist(hist_data, 50, density=False, facecolor='green')

plt.xlabel('Allele frequency') #plt.xlabel(df['AF'].iloc[0])
plt.ylabel('Frequency')
plt.title('Histogram Example')

plt.grid(True)

plt.show()


# In[76]:


hist_data = df['AF_GS'].values
# the histogram of the data
plt.hist(hist_data, 50, density=False, facecolor='green')

plt.xlabel('Allele frequency') #plt.xlabel(df['AF'].iloc[0])
plt.ylabel('Frequency')
plt.title('Histogram Example')

plt.grid(True)

plt.show()


# In[77]:


hist_data = df['AF_KF'].values
# the histogram of the data
plt.hist(hist_data, 50, density=False, facecolor='green')

plt.xlabel('Allele frequency') #plt.xlabel(df['AF'].iloc[0])
plt.ylabel('Frequency')
plt.title('Histogram Example')

plt.grid(True)

plt.show()


# In[78]:


hist_data = df['AF_KS'].values
# the histogram of the data
plt.hist(hist_data, 50, density=False, facecolor='green')

plt.xlabel('Allele frequency') #plt.xlabel(df['AF'].iloc[0])
plt.ylabel('Frequency')
plt.title('Histogram Example')

plt.grid(True)

plt.show()


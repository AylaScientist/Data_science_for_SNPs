# python 3.7
# QC.py

"""
python 3.7

    @autor : Ayla
    @version : 0.1
    @date : July 2020

This script is a quality control of the data filtered by the script "ASE_workflow.py".
These SNPs are unbiased and biallelic, and the monoallelic expression (MAE) has been discarded.
This script will plot the distribution of the data in general and for each experimental group.
It needs a file with the name of the experimental group as the acronyms used in the sample name. For example:
    The sample number 1 corresponds to the gills exposed to freshwater
    The two experimental factors are "gills" (G) and "freshwater" (F)
    Then the sample is called "1GF" and the group name will be called "GF"
This script also needs a csv file with the first column naming the experimental groups and the next columns with the
sample names of each sample in a group
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter


def allele_freqs(group_names, df):
    # Calculate the allele frequency for each experimental group in order to make a plot of each group:
    for group_name in group_names:
        af_group = "AF_" + group_name
        x_name = af_group + "_x"
        y_name = af_group + "_y"
        x = df.loc[:, df.columns.str.contains(group_name + "_R_.AD")]
        y = df.loc[:, df.columns.str.contains(group_name + "_A_.AD")]
        df[x_name] = x.sum(axis=1)
        df[y_name] = y.sum(axis=1)
        df[af_group] = df[x_name]/(df[x_name] + df[y_name])
        df.drop(columns=x_name, inplace=True)
        df.drop(columns=y_name, inplace=True)

    df = df.fillna(0)
    print("Allele frequencies by experimental group calculated")
    return df




def plot(df, group_names):
    # Plot all the histograms with the allele frequencies

    for group_name in group_names:
        af_group = "AF_" + group_name

        # Possibility 1
        # df.hist(column=af_group, bins=100, grid=False, figsize=(8,10), layout=(1,1), sharex=True, color='#86bf91', zorder=2, rwidth=0.9)

       # Possibility 3 previous in the script
        hist_data = df[af_group].values
        plt.hist(hist_data, 100, density=False, facecolor='#86bf91')
        plt.xlabel('Allele frequency')
        plt.ylabel('Frequency')
        plt.title(af_group)
        plt.grid(True)
        axes = plt.gca()
        axes.set_ylim([0, 35])
        plt.show()


def main():
    # Read the working path
    PATH = os.getcwd()

    # Import the files
    df = pd.read_csv("SNPs_ready.csv")
    groups_df = pd.read_csv("Experimental_groups.csv")

    # Create a numpy array of arrays with the samples of each group and the total samples of the experiment
    groups = groups_df.values

    # Create a numpy array with the name of each group
    group_names = list(groups_df["Group"])

    """ 
    ALLELE FREQUENCIES
    ------------------
    Compute allele frequencies for each experimental group.
    """
    df_qc = allele_freqs(group_names, df)

    # Copy the quality control file for further calculation of the stats
    df_qc.to_csv("QC.csv")

    # Plot the histograms
    plot(df_qc, group_names)

if __name__ == '__main__':
    main()
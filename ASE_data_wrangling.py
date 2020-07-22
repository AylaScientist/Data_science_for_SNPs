# Python 3.7
# ASE_workflow.py

"""
python 3.7

    @autor : Ayla
    @version : 0.1
    @date : July 2020

In this script we will process two tables containing SNPs that have been called against two sample-based pseudogenomes
from the same species. The input is produced after annotation of the vcf file with annovar and extraction of the table
using the "Variant to table" tool from the GATK package. The variants are biallelic and the structure of the csv file
includes four columns for each sample: allele depths from recessive and alternative alleles and genotypes from recessive
and alternative alleles.
"""

import pandas as pd
import numpy as np
import os
from os import path


# Make an iterator for the collection of multiallelic sites:
def multiallelic_sample(df, sample, PSGs):
    # For obtaining the column names:
    r = "_R_"
    a = "_A_"
    gt = ".GT"
    PSG1_code: str = PSGs[0]
    PSG2_code: str = PSGs[1]
    rPSG1_GT = int(df.columns.get_loc(str(sample + r + PSG1_code + gt)))
    aPSG1_GT = int(df.columns.get_loc(str(sample + a + PSG1_code + gt)))
    rPSG2_GT = int(df.columns.get_loc(str(sample + r + PSG2_code + gt)))
    aPSG2_GT = int(df.columns.get_loc(str(sample + a + PSG2_code + gt)))
    print("Multiallelic loop in sample: ", sample)

    sample_index = []
    for i in range(len(df)):
        if (((df.iloc[i, rPSG1_GT]) != df.iloc[i, rPSG2_GT]) & (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG2_GT])) | (
                (df.iloc[i, aPSG1_GT] != df.iloc[i, rPSG2_GT]) & (df.iloc[i, aPSG1_GT] != df.iloc[i, aPSG2_GT])):
            sample_index.append(i)
    return sample_index


def multiallelic(df, samples, PSGs, multi_index):
    for sample in samples:
        multiallelic_sample_index = multiallelic_sample(df, sample, PSGs)
        multi_index = multi_index + multiallelic_sample_index
    df_bi = df.drop(index=multi_index)
    return df_bi


def evaluation(df, sample, PSGs):
    # Evaluation of genotype and average. This is an iterator
    r = "_R_"
    a = "_A_"
    gt = ".GT"
    ad = ".AD"
    PSG1_code: str = PSGs[0]
    PSG2_code: str = PSGs[1]
    rPSG1_GT = int(df.columns.get_loc(str(sample + r + PSG1_code + gt)))
    aPSG1_GT = int(df.columns.get_loc(str(sample + a + PSG1_code + gt)))
    rPSG2_GT = int(df.columns.get_loc(str(sample + r + PSG2_code + gt)))
    aPSG2_GT = int(df.columns.get_loc(str(sample + a + PSG2_code + gt)))
    rPSG1_AD = int(df.columns.get_loc(str(sample + r + PSG1_code + ad)))
    aPSG1_AD = int(df.columns.get_loc(str(sample + a + PSG1_code + ad)))
    rPSG2_AD = int(df.columns.get_loc(str(sample + r + PSG2_code + ad)))
    aPSG2_AD = int(df.columns.get_loc(str(sample + a + PSG2_code + ad)))

    # Create the empty lists to collect the average variables and genotypes
    r_AD = []
    a_AD = []
    r_GT = []
    a_GT = []

    print("Evaluation loop in sample: ", sample)

    # Start iterator
    for i in range(len(df)):
        # Case heterozygots
        if (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]) and (df.iloc[i, aPSG1_GT] == df.iloc[i, aPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = ((df.iloc[i, rPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2)
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            y = 0
            a_AD.append(y)

        # Case 2 Same situation Ref and Alt in both for df
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]) and (df.iloc[i, aPSG1_GT] == df.iloc[i, aPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = (df.iloc[i, rPSG1_AD] + df.iloc[i, rPSG2_AD]) / 2
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG1_GT])
            y = (df.iloc[i, aPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2
            a_AD.append(y)
        # Case 3 Inverse situation Ref and Alt in both for df
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG2_GT]) and (df.iloc[i, aPSG1_GT] == df.iloc[i, rPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG1_GT])
            y = (df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD]) / 2
            a_AD.append(y)

        # Case 4 First homozygot, second with alternative allele in second position
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] != df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD] + df.iloc[
                i, rPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 individuals
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG2_GT])
            y = df.iloc[i, aPSG2_AD]
            a_AD.append(y)

        # Case 5 First homozygot, second with alternative allele in first position of the second genome
        elif ((df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] != df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG2_GT])):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD] + df.iloc[
                i, aPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 individuals
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            y = df.iloc[i, rPSG2_AD]
            a_AD.append(y)

        # Case 6 Second homozygot, first with alternative allele in second position of the first genome
        elif (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = df.iloc[i, rPSG1_AD]
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG1_GT])
            y = (df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[
                i, aPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 individuals
            a_AD.append(y)

        # Case 7 Second homozygot, first with alternative allele in first position of the first genome
        elif (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = df.iloc[i, rPSG1_AD]
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            y = (df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[
                i, aPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 individuals
            a_AD.append(y)

        else:
            print("Error in df row ", i, ", sample num ", sample)

    return (r_GT, r_AD, a_GT, a_AD)


def sample_average(df, samples, PSGs):
    cols = []  # Index of columns to be drop in the end
    # Create 4 columns for each averaged sample:
    for sample in samples:
        print("Sample average in sample ", sample,"\n")
        # Average each sample according to genotype:
        (r_GT, r_AD, a_GT, a_AD) = evaluation(df, sample, PSGs)

        # Add the new columns
        sample_r_GT = str(sample + "_R_.GT")
        sample_a_GT = str(sample + "_A_.GT")
        sample_r_AD = str(sample + "_R_.AD")
        sample_a_AD = str(sample + "_A_.AD")
        df[sample_r_GT] = r_GT
        df[sample_a_GT] = a_GT
        df[sample_r_AD] = r_AD
        df[sample_a_AD] = a_AD

        # Drop the columns that won't be used
        r = "_R_"
        a = "_A_"
        gt = ".GT"
        ad = ".AD"
        PSG1_code: str = PSGs[0]
        PSG2_code: str = PSGs[1]
        rPSG1_GT = int(df.columns.get_loc(str(sample + r + PSG1_code + gt)))
        aPSG1_GT = int(df.columns.get_loc(str(sample + a + PSG1_code + gt)))
        rPSG2_GT = int(df.columns.get_loc(str(sample + r + PSG2_code + gt)))
        aPSG2_GT = int(df.columns.get_loc(str(sample + a + PSG2_code + gt)))
        rPSG1_AD = int(df.columns.get_loc(str(sample + r + PSG1_code + ad)))
        aPSG1_AD = int(df.columns.get_loc(str(sample + a + PSG1_code + ad)))
        rPSG2_AD = int(df.columns.get_loc(str(sample + r + PSG2_code + ad)))
        aPSG2_AD = int(df.columns.get_loc(str(sample + a + PSG2_code + ad)))

        cols.append(rPSG1_GT)
        cols.append(aPSG1_GT)
        cols.append(rPSG2_GT)
        cols.append(aPSG2_GT)
        cols.append(rPSG1_AD)
        cols.append(aPSG1_AD)
        cols.append(rPSG2_AD)
        cols.append(aPSG2_AD)

    return df, cols


def AD10(df, samples):
    # The iterator collects the index were AD from both alleles is below 10
    # Keep all this analysis for the same individual including two tissues at a time
    i = 0
    indexes_ad10 = 0
    if str(path.exists("Samples_MAE.csv")):
        tissues = pd.read_csv("Samples_MAE.csv")
        tissues = pd.DataFrame(tissues)
        first_samples = tissues.iloc[:, 0]
        second_samples = tissues.iloc[:, 1]

        for first_sample in first_samples:
            print("AD10 filter in sample ", first_sample, " and sample ", second_samples[i])
            sample1_r_ad = str(first_sample + "_R_.AD")
            sample1_a_ad = str(first_sample + "_A_.AD")
            sample2_r_ad = str(second_samples[i] + "_R_.AD")
            sample2_a_ad = str(second_samples[i] + "_A_.AD")
            index_sample = df[
                (((df[sample1_a_ad] + df[sample1_r_ad]) < 10) | ((df[sample2_a_ad] + df[sample2_r_ad]) < 10))].index
            if i == 0:
                indexes_ad10 = np.array(index_sample)
                print("Number of rows to drop for sample ",first_sample," ", len(indexes_ad10))
                i = i + 1
            elif i != 0:
                index_sample = np.array(index_sample)
                indexes_ad10 = np.intersect1d(indexes_ad10, index_sample)
                print("Number accumulated of rows to drop for sample ",first_sample," ", len(indexes_ad10))
                i = i + 1
    else:
        for sample in samples:
            print("AD10 filter in sample ", sample)
            sample_r_ad = str(sample + "_R_.AD")
            sample_a_ad = str(sample + "_A_.AD")

            index_sample = df[(df[sample_a_ad] + df[sample_r_ad]) < 10].index
            if i == 0:
                indexes_ad10 = np.array(index_sample)
                print("Number of rows to drop", len(indexes_ad10))
                i = i + 1
            elif i != 0:
                index_sample = np.array(index_sample)
                indexes_ad10 = np.intersect1d(indexes_ad10, index_sample)
                print("Number of rows to drop", len(indexes_ad10))
                i = i + 1

    df.drop(index=indexes_ad10, inplace=True)
    return df


def compare(df, sample, sample1rgt):
    # This iterator will construct the new columns for the samples according to the reference and alternative alleles
    # assigned in sample1

    # Position of the sample columns
    sample_r_gt = int(df.columns.get_loc(str(sample + "_R_.GT")))
    sample_a_gt = int(df.columns.get_loc(str(sample + "_A_.GT")))
    sample_r_ad = int(df.columns.get_loc(str(sample + "_R_.AD")))
    sample_a_ad = int(df.columns.get_loc(str(sample + "_A_.AD")))

    # Create the empty lists to collect the average variables and genotypes
    r_AD = []
    a_AD = []
    r_GT = []
    a_GT = []

    print("Compare genotypes loop in sample: ", sample)

    # Start iterator
    for i in range(len(df)):
        if (df.iloc[i, sample1rgt] == df.iloc[i, sample_r_gt]):
            r_AD.append(df.iloc[i, sample_r_ad])
            a_AD.append(df.iloc[i, sample_a_ad])
            r_GT.append(df.iloc[i, sample_r_gt])
            a_GT.append(df.iloc[i, sample_a_gt])
        else:
            r_AD.append(df.iloc[i, sample_a_ad])
            a_AD.append(df.iloc[i, sample_r_ad])
            r_GT.append(df.iloc[i, sample_a_gt])
            a_GT.append(df.iloc[i, sample_r_gt])

    return (r_AD, a_AD, r_GT, a_GT)


def genotype(df, samples):
    # Determines the genotypes for sample1 and calls the function to compare genotypes
    i = 0
    sample1rgt = 0
    for sample in samples:
        i = i + 1
        if i == 1:
            sample1rgt = int(df.columns.get_loc(str(sample + "_R_.GT")))
        elif i != 1:
            (r_ad, a_ad, r_gt, a_gt) = compare(df, sample, sample1rgt)

            # Collection the name of the samples
            sample_r_gt = str(sample + "_R_.GT")
            sample_a_gt = str(sample + "_A_.GT")
            sample_r_ad = str(sample + "_R_.AD")
            sample_a_ad = str(sample + "_A_.AD")

            # Assign the list to the column name
            df[sample_r_gt] = r_gt
            df[sample_a_gt] = a_gt
            df[sample_r_ad] = r_ad
            df[sample_a_ad] = a_ad
    return df


def frequencies(df, samples):
    """

    :param df: dataframe with unbiased SNP counts without allele frequencies
    :param samples: list of strings including the sample names
    :return: dataframe with the allele frequencies by sample
    """
    # Calculate the allele frequencies
    for sample in samples:
        sample_r_ad = str(sample + "_R_.AD")
        sample_a_ad = str(sample + "_A_.AD")
        af_sample = str("AF_" + sample)
        df[af_sample] = df[sample_r_ad] / (df[sample_r_ad] + df[sample_a_ad])
    df = df.fillna(0)  # To fill cases where there are no counts and therefore AF is divided by 0
    print("Allele frequencies by sample calculated")

    return df


def MAE(df, samples):
    # Collect the indexes whose allele frequencies are 0 or 1 in both tissues of the same individual. These indexes
    # correspond to the SNPs that are MAE for at least one sample
    if str(path.exists("Samples_MAE.csv")):
        tissues = pd.read_csv("Samples_MAE.csv")
        tissues = pd.DataFrame(tissues)
        first_samples = tissues.iloc[:, 0]
        second_samples = tissues.iloc[:, 1]
        i = 0
        for first_sample in first_samples:
            print("MAE loop in sample ", first_sample, " and sample ", second_samples[i])
            af1 = str("AF_" + first_sample)
            af2 = str("AF_" + second_samples[i])
            index_af = df[((df[af1] == 0) & (df[af2] == 0)) | (df[af1] == 1) & (df[af2] == 1)].index
            df.drop(index=index_af, inplace=True)
            i = i + 1
    else:
        for sample in samples:
            af = str("AF_" + sample)
            index_af = df[(df[af] == 0) | (df[af] == 1)].index
            df.drop(index=index_af, inplace=True)

    return df


def main():
    """
    Input the files and convert them into data frames
    """

    # Read the dataframe and convert it into pandas dataframe. PSG stands for pseudogenome:

    PATH = os.getcwd()
    """
    print(PATH)
    PSG1_name: str = input("Please enter the file name of the AD_GT_counts_bi for the first pseudogenome:\n")
    PSG2_name: str = input("Please enter the file name of the AD_GT_counts_bi for the second pseudogenome:\n")
    
    PSG1 = pd.read_csv(PSG1_name, low_memory=False)
    PSG2 = pd.read_csv(PSG2_name, low_memory=False)
    """
    PSG1 = pd.read_csv("AD_GT_counts_bi_GF3.csv", low_memory=False)
    PSG2 = pd.read_csv("AD_GT_counts_bi_KF6.csv", low_memory=False)
    PSG1 = pd.DataFrame(PSG1)
    PSG2 = pd.DataFrame(PSG2)

    # Read the sample names and the codes for each pseudogenome:
    """sample_names = pd.read_csv(input("Please enter the file name that has the sample names:\n"))
    PSG_codes = pd.read_csv(input("Please enter the file name that has the pseudogenomes codes:\n"))"
    """
    sample_names = pd.read_csv("Sample_names.csv")
    PSG_codes = pd.read_csv("Pseudogenome_codes.csv")

    # Create arrays of the sample names and the pseudogenome codes
    samples = sample_names['Sample_name'].values
    PSGs = PSG_codes['PSGs'].values

    """
     MERGE THE DATAFRAMES 
    -----------------------
    
    The format of the table is the same for both dataframes. Columns from 1 to 6 are:
        CHROM: Chromosome,
        POS: position,
        Gene.refGene: Gene_ID
        Func.refGene: Annotation of the function for this SNP
        ExonicFunc.refGene: Annotation of the function if this SNP is exonic
        AF: Allele frequency estimated for all the samples
    
    From column 7th till the end there will be 4 columns for each sample. The structure of the column name is
    as follows:
        NAME (determined by the sample name)
        _R or _A (reference or alternative allele)
        _??? (The pseudogenome code)
        .AD (allele depth)
        .GT (Genotype)
    
    The next step is to merge the two dataframes on the common chromosome and position:
    """

    df = pd.merge(PSG1, PSG2, on=('CHROM', 'POS'))
    print('Files merged')

    """
     DELETE MULTIALLELIC SITES 
     -------------------------
    
    Compare the genotypes for assessing the counts on each allele
    Set the data for the mapping against the pseudogenome of sample PSG1 because it hast the variants of the reference genotype. The SNPs called by this sample are set as reference and alternative alleles for all the other samples.
    
    The tables containing the SNPs from mapping to PSG1 and PSG2 contained only biallelic sites. However multiallelic possibilities can appear if one SNP was called for a different genotpe in each of the pseudogenomes. The multiallelic sites will be deleted from this dataframe to continue the analysis with bialleleic sites only.
    
    Delete multiallelic sites
    Let's start droping the multiallelic. The criteria compares if the reference allele in the SNPs resulting from mapping against PSG1 pseudogenome is different form the reference and the alternative alleles from the mapping against PSG2 pseudogenome. This is later repeated for the alternative allele. Note that the SNPs mapped against PSG1 pseudogenome have been pre-filtered for biallelic sites for all the samples.
    
    Collect the indexes where the reference allele in PSG1 is different from both alleles in the other mappings or the alternative allele in PSG1 is different from both alleles in the PSG2 mapping.
    """

    # Call the function to clean the multiallelic sites and drop the dfs with this condition:
    multi_index = []
    df_bi = multiallelic(df, samples, PSGs, multi_index)

    # Make a temporary folder where you can copy temporary files for backups
    # os.mkdir("temp")
    os.chdir("temp")
    df_bi.to_csv('Biallelic_SNPs.csv')
    os.chdir(PATH)

    # From now on we work with df_bi standing for "dataframe_biallelic_SNPs"

    """
    COMPARE GENOTYPES AND MAKE THE AVERAGE OF COUNTS
    ------------------------------------------------
    
    Now let's make the average on the counts for the reference and the alternative alleles. For each sample, select the 
    genotype of the reference allele established by the reads mapped against PSG1 and compare it to both alternative and 
    reference genotypes of the values in the mapping with PSG2 genome. Then proceed to make the average of the counts. The 
    possible cases are the next: 
    
    1.- Both homozigots for the reference allele: In that case set the same genotype to both alleles and sum all the 
    counts from both sites, making the average for two individuals. 
    
    2.- Reference and alternative allele in the SNP from the PSG1 pseudogenome are placed in the same position for the SNP 
    from the PSG2 pseudogenome. Set the genotype to the reference and alternative alleles to the ones from the SNP in PSG1 
    and make the average of the counts for reference and alternative genome 
    
    3.- The reference allele in the reference pseudogenome is called as alternative allele in the alternative 
    pseudogenome. In that case, set the genotype to the reference allele and make the average of the counts for the 
    corresponding variant. 
    
    4.- While the mapping against the reference genome resulted in an homozygot for the reference allele, the mapping 
    against the alternative genome resulted in an heterozygot with the alternative allele in the alternative position. 
    The genotype of the reference allele will be set from genotype of the reference allele in the reference mapping and 
    the alternative allele will be set from the genotype of the alternative allele in the alternative mapping. The counts 
    of each allele will be set accordingly to this distribution, thus dividing the counts by two individuals in the 
    reference allele. 
    
    5.- The mapping against the reference genome resulted in an homozygot for the reference allele and the mapping 
    against the alternative genome is heterozygot. In that case, the genotype of the alternative allele is set as the 
    reference allele on the mapping against the alternative genome. The genotypes must be set accordingly, 
    being the reference allele the one found in the reference genome and the alternative allele the one found in the 
    alternative genome. The average of the counts will follow this pattern, being the average in the erference allele. 
    
    6.- The mapping against the reference genome produced an heterozygot with reference and alternative alleles. The 
    mapping against the alternative pseudogenome produced an homozygot for the reference allele. In that case, 
    reference and alternative allele must be set according to the reference pseudogenome and the counts must be averaged 
    for both individuals in the alternative allele. 
    
    7.- The mapping against the reference genome produced an heterozygot with reference and alternative alleles. The 
    mapping against the alternative pseudogenome produced an homozygot for the alternative allele. In that case, 
    reference and alternative allele must be set according to the reference pseudogenome and the counts must be averaged 
    for both individuals in the alternative allele. 
    
    This analysis is repeated for each sample and will assign reference and alternative alleles independently. The 
    assignation of reference and alternative alleles uniformly in all the samples will be developed in a further step of the 
    workflow.
    
    """

    df_average, cols = sample_average(df_bi, samples, PSGs)

    # Drop the columns that are not needed
    df_average.drop(df_average.columns[cols], axis=1, inplace=True)

    os.chdir("temp")
    df_average.to_csv("Average_SNPs.csv")
    os.chdir(PATH)

    """
    DELETE SNPs WITH AD<10
    ----------------------
    This new function cleans the SNPs that do not have enough counts and are considered possible bad reads
    """

    df_AD10 = AD10(df_average, samples)

    # Write the temporary file
    os.chdir("temp")
    df_AD10.to_csv("AD10_SNPs.csv")
    os.chdir(PATH)

    """
    ASSIGN REFERENCE AND ALTERNATIVE ALLELES UNIFORMLY IN ALL SAMPLES
    -----------------------------------------------------------------
    Reference and alternative alleles have been assigned by sample. This function will correct it and will assign the 
    reference and alternative alleles according to the pattern established in sample 1.
    """

    df_uni = genotype(df_AD10, samples)

    # Write the temporary file
    os.chdir("temp")
    df_uni.to_csv("Uniform_SNPs.csv")
    os.chdir(PATH)

    """ 
    ELIMINATE MONOALLELIC EXPRESSION
    --------------------------------
    The monoallelic expression (MAE) can result from homozygots as well as imprinted genes. In order to distinguish each
    case, we need to know the genotype. However, that's not possible for the amount of SNPs that we are working with.
    Therefore we drop all the SNPs whose allele frequency is 0 or 1.
    
    If the samples are taken from different organs from the same individual, the system needs a file called 
    "Samples_MAE.csv" where the sample name from an individual are in the same line and there are two columns, 
    one for each tissue we need to compare. The name of the columns must be in plural.
    """

    df_af = frequencies(df_uni, samples)

    df_mae = MAE(df_af, samples)

    df_mae.to_csv("SNPs_ready.csv")


if __name__ == '__main__':
    main()

# DATA SCIENCE FOR SNPs

# AUTHOR
Aurora Campo, ayla.bcn@gmail.com

# HARDWARE/SOFTWARE REQUIREMENTS
    x86-64 compatible processors
    Python 3.7



# INSTRUCTIONS
In this git you can find the scripts for managing SNP counts retrieved from GATK after mapping against two pseudogenomes.

The order of the files is:
	ASE_data_wrangling.py 
	QC.py
	Stats.py

Input files are csv files including reads and genotypes for reference and alternative alleles separate in columns for each sample.

Output will be Allelic Imbalance and Allele Specific Expression for experimental groups determined by Chi-square, binomial test and Fisher exact test.


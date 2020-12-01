# NbdLamprey
Pipeline to process rapture genotypes from a raw vcf file to filtered SNP sets for Nb/Ns estimates in COLONY and NeEstimator. Contains a test data set to work through all of the scripts before scaling up for a large population. There is a markdown tutorial with more details about scripts as well as output summaries for the test data set.

## Guide to navigating repository
### Folder overviews:
Input: initial inputs for the pipeline (pre-processed by VCFtools)
Summary Statistics: Outputs from the Summary Statistics scripts
SNPsets: Outputs from the SNP_set_creation script
Software Outputs: Outputs from Colony to be used for cohorts and Nb calculations
Output: Final outputs for the pipeline
Homebrew: contains all the user-created functions that are used in the pipeline


### Script overviews:
1. Generating summary statistics for the input genotype file
2. Generating a set of SNPs for COLONY for each location
3. Length-based Aging Models for each location
4. Generating COLONY files and NeEstimator SNP set for all cohorts
5. Calculating Nb and Ns for all cohorts
6. Additional Scripts to create figures found in the manuscript


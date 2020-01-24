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
Scripts are a combination of bash command line code (summarized in Rmarkdown file) and R scripts
1. VCFtools_processing - bash
2. Summary Statistics - R
2a. ***Figures and Tables for Summary Statistics - R
3. SNP_set_creation - R
4. Cohort_determination - R
5. Nb_Ns_calculation - R
S1. PCA - R


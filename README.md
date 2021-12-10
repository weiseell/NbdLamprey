# NbdLamprey
Pipeline to process RAD-capture genotypes from a vcf file to filtered SNP sets for Nb/Ns estimates in COLONY and NeEstimator. Contains a test data set to work through all of the scripts before scaling up for a large population. There is a markdown tutorial with more details about scripts as well as output summaries for the test data set.

## Guide to navigating repository
### Folder overviews:
Input: initial inputs for the pipeline (pre-processed by VCFtools)
Aging_Models: Data inputs for length mixture models, and output containing cohort assignments
Summary_Stats: Outputs from the Summary Statistics scripts
SNPsets: Outputs from the SNP_set_creation script
Software_outputs: Outputs from Colony to be used for cohorts and Nb calculations
Homebrew: contains all the user-created functions that are used in the pipeline. Each function has comments at the top describing each function input for easier use.
Figures: Folder containing .tiff files of final figures generated for the manuscript Weise et al. in review

### Script overviews:
1. Generating summary statistics for the input genotype file
1A. Additional script to test for SNPs deviant from Hardy-Weinburg expectations across collections
2. Generating a set of SNPs for COLONY for each location
3a. Mixture models to generate cohort assignment using length data for each location
3b. Combining results from script 3a with Colony clusters generated from the output of script 2 to correct mixture models
4. Generating COLONY files and NeEstimator SNP set for all cohorts
5. Calculating Nb and Ns for all cohorts
6. Additional Scripts to create figures found in the manuscript
S1. Extra figures and summary calculations to further summarize data by populations
S2. PCA analysis to identify potential native lamprey individuals sequenced with sea lamprey (if any are present). Known sequences for brook and northern lamprey, along with known sea lamprey, are used to calibrate the PCA.


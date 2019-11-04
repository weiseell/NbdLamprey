# NbdLamprey
Pipeline to process rapture genotypes from a raw vcf file to filtered SNP sets for Nb/Ns estimates in COLONY and NeEstimator

## Guide to navigating repository
### Folder overviews:
Input: contains initial data and intermediate data frame outputs used as inputs in later scripts
Output: will contain final file outputs of the pipeline as well as several figures throughout pipeline that are generated
Homebrew: contains all the user-created functions that are used in the pipeline


### Script overviews:
Scripts are a combination of bash command line code (summarized in Rmarkdown file) and R scripts
0. Summary of bioinformatics pipeline to process raw Illumina reads into a .vcf file with called genotypes (scripts from Sard et al 2019 rad genotyping pipeline)
1. VCFtools_processing
2. vcf_prep
3. ld_filtering
4. quality_filtering
5. COLONY_file_creation
6. GENEPOP_file_creation
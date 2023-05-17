# Genetic diversity varies with species traits and latitude in predatory soil arthropods (Myriapoda: Chilopoda)
---



## Summary of contents

The folders in this repository contain raw data and R scripts to reproduce the results in our manuscript, where we investigate the correlates of intra-specific genetic diversity in centipedes using species traits and biogeography as explanatory variables. The scripts can be used to read in the raw data with sequence and coordinate information, perform species-wise sequence alignments, calculate alignment summary statistics, create input files for beta regression and perform sample size sensitivity analysis. Figures 2 to 5 and Table 1 in the main manuscript, and all the figures and tables in the Supplementary Information can be reproduced using the data and scripts provided.


## Input data and folder structure

To reproduce the analysis, download the folders to a suitable directory. Change the working directory in the beginning of each R script within the 'scripts' folder and run the R code sequentially based on the file name.

### 1. scripts 
This folder contains the R code required for data processing, analysis and visualization. The scripts are numbered sequentially as code01 to code05 based on the order in which they should be run.

### 2. data_raw
This folder contains all the raw data files needed for analysis.

The sequence and coordinate information is present in 'sheet4_popgen_database_analysis_clean_11Apr23.csv' and the species trait information is present in 'centipedes_life_history_11Apr23.csv'.
    
### 3. data
This folder contains all the processed data files generated using the R scripts.

This includes the files 'input_betareg_no_introductions_11Apr23.csv' and 'input_betareg_11Apr23.csv', which contain the input data for beta regression analysis.
    
### 4. figures
This folder contains all the main and supplementary figures generated using the R scripts.

### 5. results
This folder contains all the tables generated using the R scripts.

## Analysis code

The R scripts present in the 'scripts' folder are:
### 1. code01_sequence_alignment_11Apr23.R
This code is to clean and process raw data, query accession numbers on GenBank, perform multiple sequence alignment for each species and save these files to disk.

### 2. code02_phylogatr_sequences_11Apr23.R
This code is to processes phylogatr data, compare our database with the phylogatr database, look at common and different accession numbers between the two and save results.

### 3. code03_sequence_summary_11Apr23.R
This code is to calculate sequence summaries from the multiple sequence alignment of each species, combine species trait and biogeography information and create input files for statistical analysis.

### 4. code04a_analysis_plots_no_introductions_19Apr23.R
This code is to prepare input data for analysis, run beta regression models, perform model diagnostics, correct for spatial autocorrelation in model residuals, obtain bootstrapped model coefficients, test for phylogenetic signal in model residuals and create figures and tables. 

### 5. code04b_analysis_plots_11Apr23.R
Same analysis pipeline as above but using data including synanthropic introductions. 

### 6. code05_sample_size_sensitivity_11Apr23.R
This code is to resample sequence alignments for each species iteratively across a range of sample sizes, calculate the variance in estimated genetic diversity for each sample size and run beta regression models using the optimal sample size.

## Sharing and access information

### Please use the following citation for the use of the data compilation or analysis code:

Bharti, D. K., Pawar, P. Y., Edgecombe, G. D., & Joshi, J. (2023). Genetic diversity varies with species traits and latitude in predatory soil arthropods (Myriapoda: Chilopoda). Global Ecology and Biogeography.

### Links to other publicly accessible locations of the data:
  * https://datadryad.org/stash/share/Z8tG0lFyS-DBIMgLbAQqF7FNIs1n9pwciblORMsI1M0

### Detailed information regarding data sources is available in the Supplementary Information associated with our manuscript:
  * Appendix S1: Mitochondrial COI accession numbers and associated sequence coordinates used for data analysis
  * Appendix S2: Centipede species traits and biogeographic variables used for data analysis

Supplementary Information is available at https://doi.org/10.5281/zenodo.7940353



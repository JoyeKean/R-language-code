# R-language-code
Mendelian randomization of the original code
# 🚪Project Overview： Project Overview: This codebase is used to analyse the causal relationship between cerebrospinal fluid (CSF) metabolites and ovarian disease based on Mendelian randomization (MR) method, and is planned to include data preprocessing, MR analysis, and visualization scripts. However, our team does not have the programming ability for the time being, the content of this writing is obtained with the help of artificial intelligence and temporary learning, we are actively contacting the relevant technical staff, and will continue to improve the content, and strive to run smoothly as soon as possible~.
# 📦 Dependencies and Installation：
1. R language environment
R version: ≥ 4.0.0
Dependent Package Installation: Run the following code to install the required packages: install.packages(c("TwoSampleMR", "MendelianRandomization", "ggplot2", "dplyr", "readr"))
2. Data Preparation
Input file:
exposure： All metabolites are available in the GWAS database (https://www.ebi.ac.uk/gwas/), with the code range from ebi-a-GCST90026002 to ebi-a-GCST90026289.
outcome: The GWAS data of four common ovarian diseases and infertility among women of childbearing age were sourced from the Finngen database (https://www.finngen.fi/en/).
#The following are the sources of the ending factors
phenotype	database	object of study	cases	controls	Gwas-id
Polycystic ovarian syndrome	FinnGen	EUR	34388	195922	finn-b-E4_POCS
Salpingitis and oophoritis	FinnGen	EUR	6814	226439	finn-b-N14_SALPHOOPH
Ovarian dysfunction	FinnGen	EUR	2309	218970	finn-b-E4_OVARDYS
Benign neoplasm of ovary	FinnGen	EUR	6323	248295	finn-b-CD2_BENIGN_OVARY
Female infertility, cervigal, vaginal, other or unspecified origin	FinnGen	EUR	14407	119468	finn-b-N14_FIOTHNAS
Female infertility	FinnGen	EUR	16720	119468	finn-b-N14_FEMALEINFERT
# 🍵Instructions for installing and using the R package, including examples:
The following are the methods for installing and running the main R packages in this study.
1. TwoSampleMR
#To install the TwoSampleMR package in R, follow these steps.
1) Install the devtools package: install.packages(‘devtools’)
2) Load the devtools package: library(devtools)
3) Install the TwoSampleMR package from GitHub: devtools::install_github(‘MRCIEU/TwoSampleMR’)
#The main purpose of this package is to help the user to extract genetic information on instrumental variables (SNPs) from two independent samples, to perform causal inference analyses, and to assess the causal relationship between exposure and outcome.
Example:
4) Load the package using:library(TwoSampleMR)
2. ggplot2
Install devtools: install.packages(‘devtools’)
library(ggplot2)
#Create graphs with ggplot(data, aes(x, y)) + geom_xxx().
#ggplot2 can plot scatterplots, histograms, histograms, regression curves, and more.
3. foreach
install.packages(‘foreach’)
library(foreach)
#The main function of doParallel is to compute the causal estimation of multiple SNPs in parallel, which improves the computation speed and is especially suitable for MR analysis of large-scale GWAS data. Combined with doParallel, we can compute IVW, MR-Egger, Weighted Median estimation of multiple SNPs at the same time, which greatly improves the efficiency of analysis.
Example: # Load the necessary packages
library(MendelianRandomisation)
library(foreach)
library(doParallel)
#Create a parallel computing environment
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
#Assume we have MR input data from multiple SNPs
beta_exposure <- c(0.1, 0.05, 0.02, 0.15, 0.07)
beta_outcome <- c(0.3, 0.2, 0.1, 0.4, 0.25)
se_exposure <- c(0.02, 0.03, 0.01, 0.02, 0.04)
se_outcome <- c(0.05, 0.06, 0.04, 0.05, 0.07)
#Parallel computation of IVW, MR-Egger and Weighted Median estimates
mr_results <- foreach(i = 1:length(beta_exposure), .combine = rbind) %dopar% {
 #Create the MR input object
  mr_input <- mr_input(bx = beta_exposure[i], bxse = se_exposure[i],
                        by = beta_outcome[i], byse = se_outcome[i])
 #Calculate the IVW estimate
  ivw_result <- mr_ivw(mr_input)
 #Calculate MR-Egger estimate
  egger_result <- mr_egger(mr_input)
 #Calculate Weighted Median estimate
  weighted_median_result <- mr_median(mr_input)
  #Return the result
  c(ivw_result@Estimate, ivw_result@Pvalue,
    egger_result@Estimate, egger_result@Pvalue, weighted_median_result
    weighted_median_result@Estimate, weighted_median_result@Pvalue)}
#Shut down the parallel computing environment
stopCluster(cl)
#Set the column names and print the results
colnames(mr_results) <- c(‘IVW_Estimate’, ‘IVW_P’. 
                          ‘Egger_Estimate’, “Egger_P”. 
                          ‘WeightedMedian_Estimate’, “WeightedMedian_P”)
print(mr_results)
4. Mendelian Randomisation
install.packages(‘MendelianRandomisation’)
library(MendelianRandomisation)
Example:
#Load packages
library(MendelianRandomisation)
#Assuming you already have data related to the exposure (X) and the ending (Y), stored in data frame format
#Example data (effects of exposure and ending variables)
exposure <- c(0.1, 0.05, 0.02) # exposure effect
outcome <- c(0.3, 0.2, 0.1) # outcome effect
se_exposure <- c(0.02, 0.03, 0.01) # standard error of exposure effect
se_outcome <- c(0.05, 0.06, 0.04) # Standard error of the ending effect
#Construct a data frame
data <- data.frame(exposure, outcome, se_exposure, se_outcome)
#Perform MR analysis
mr_results <- mr_ivw(data)
#View results
print(mr_results)
5. data.table
install.packages(‘data.table’)
library(data.table)
Example:
#Load the data.table package
library(data.table)
#Create a data.table object
dt <- data.table(ID = 1:5, Name = c(‘A’, ‘B’, ‘C’, ‘D’, ‘E’), Age = c(23, 45, 34, 29, 56))
#View data
print(dt)
#Filter data: select rows with Age > 30
dt[Age > 30]
#Add a column: calculate Age plus 10 for each person
dt[, Age_plus_10 := Age + 10]
 #View the updated data.table
print(dt)
# 🚀 Quick Start Guide
1. Running a complete analysis
Execute the main script main_analysis.R to automate data loading, MR analysis, and results saving: source("scripts/main_analysis.R")
2. Step-by-step operation
Step 1: Data pre-processing
source("scripts/01_data_preprocessing.R")
Functions: cleaning data, standardizing metabolite concentrations, screening instrumental variables.
Output file
Step 2: MR analysis
source("scripts/02_mr_analysis.R")
Methods: inverse variance weighting (IVW), MR-Egger regression, weighted median method. Output file
Step 3: Visualisation source("scripts/03_visualization.R")
Generate charts: forest diagrams
# 📊 Interpretation of results
Key output files
mr_results.csv: contains the following:
method: MR analysis method (IVW, MR-Egger, etc.).
or: ratio ratio.
se: standard error.
pval: p-value.
Example results:
method or se pval
IVW 1.35 0.12 0.003
MR-Egger 1.28 0.18 0.08
Chart description
Forest plot: shows OR and confidence intervals for different MR methods.
# 📂 Code structure
├── data/ # raw data
│ ├── exposure_data.csv
│ └── outcome_data.csv
├── scripts/ # R scripts
│ ├── 01_data_preprocessing.
│ ├── 02_mr_analysis.R
│ └── 03_visualisation.R
├── results/ # Analysis results (OR, p-value, etc.)
├── figures/ # Graph output
└── README.md # This file
# ❓ Frequently Asked Questions (FAQ)
Q1：‘Package not found’ when running.
Solution: Reinstall the dependent packages, make sure to use install.packages() command.
Q2：Wrong data path
Solution: Check if exposure_data.csv is located in . /data/ directory, or modify the file path in the script.
Q3：Insufficient strength of tool variable (F-statistic <10)
Solution: Adjust SNP screening thresholds (e.g. lower p-value thresholds) in 01_data_preprocessing.R.
# 📜 Licences and citations
#📧 Contact
Maintainer:
Technical cooperation: 

# BART with Targeted Smoothing: Analysis
This repository contains code to replicate figures and tables in the paper "BART with Targeted Smoothing: An analysis of patient-specific stillbirth risk", by Starling et al (2018).

Note that many of these scripts use the **fastbart** R package for comparing tsBART with regular BART.  This package is not publicly available, but can be provided on request.  

The following are instructions for the workflow to replicate the analysis.

## Dataset.
The raw data is publicly available via the National Center for Health Statistics at the following locations:  

* 2005 data: ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/cohortlinkedus/LinkCO05US.zip  
* 2006 data: ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/cohortlinkedus/LinkCO06US.zip  

We provide several items to support replication of our process: 

* A 'data dictionary' for our dataset, which includes a link to the original data source and its documentation.     
    + Our data dictionary is located at **/data/04_stillbirth-example-data/data-dictionary.pdf**.    
    + The original dataset and documentation are located at https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm#Period_Linked        
* A sample stillbirth dataset, to illustrate the structure and format of our data set post cleaning and processing.      
    + The sample dataset is located at **data/perinatal-mortality/stillbirth-data-example.csv**. 
    + Note that this dataset will not replicate results from the paper; it is intended only to illustrate the data structure.  
* The scripts we use to clean, case control sample, and analyze the stillbirth data are located in at **code/00_stillbirth-data-processing**.  

## Data preparation and cleaning.
To use the scripts to run our analysis, the user must download the raw data files.

1. A stata script at **/Data/000-stata-script-to-create-new-variables.txt** processes the data, creating variables used in the analysis.  The user must save the output of this script in a .csv file at **/Data/BIRTH_DEATH_DATA_2004-2006.csv**.  The workflow then proceeds as follows.

2. Run the script **/code/00_stillbirth-data-processing/00a-stillbirth-process-raw-data.R** to clean the raw data and create variables for stillbirth analysis.  

* This file takes the input raw dataset **/Data/BIRTH_DEATH_DATA_2004-2006.csv**.  
    + This dataset must be saved by the user, by downloading the original publicly accessible data, and assembling in a .csv file.  
* This file outputs the processed dataset **/Data/stillbirth-clean.csv**.    

3. Run the script **/code/00_stillbirth-data-processing/00b-stillbirth-create-case-control-samples.R** to create the case-control sample of the stillbirth data.

* This file takes the input file **/Data/stillbirth-clean.csv** (created in the previous step).  
* This file outputs one file, 
**/data/perinatal-mortality/stillbirth-data-casecontrol-50-for-paper.csv**.  This file is used in the remainder of the analysis.
 

## Code execution for analysis and simulation.
The code for both stillbirth analysis and the simulation can be run by the following command.  

**bash tsbart-analysis-coderun.sh**    

This script runs several R scripts, all located in the **./code** subfolder.  These scripts can also be executed individually.

* **table-01.R**   
    + Generates Table 1 for cohort characteristics.
* **table-02.R**  
    + Repeats second simulation on 100 datasets to calculate log-likelihood.
    + Creates Table 2.    
* **figure-02-and-figure06.R**  
    + Creates Figure 2 and Figure 6 (in Appendix A1).
    + Figure 2 illustrates results of the second simulation, comparing tsBART and BART.
    + Figure 6 is example of tuning expected number of crossings.    
* **table-03-a_generate-data.R**   
* **table-03-b_one-example-run.R**   
    + These two scripts illustrate how data sets for the first simulation are generated and fit.  
    + This code does not replicate Table 3; Table 3 is the result of repeating this process on 500 datasets for each weighting scenario to calculate MSE and point-wise coverage for each method and weighting scenario combination.  
* **figure-03.R** 
    + Illustrates hazard data generated for first simulation, comparing coverage across methods.  
* **figure-04-and-table-04.R**  
    + Fits models to up-sampled stillbirth data. 
    + Creates Figure 4 and Table 4.
    + The file ending in "-single-iteration.R" illustrates the fit for a single dataset instead of all five runs.
* **figure-05.R**
    + Fits models to stillbirth data.
    + Creates Figure 5 to show estimated stillbirth risk curves for each method.
    
## Helper functions.
There are several helper functions in the analysis that are used for cleaning data, assembling tables, etc.  These are not included in the **tsbart** package, but are available in **/code/helper-functions/**.  (These are sourced by other scripts as needed.)
# Creates case-control samples for the stillbirth data for
# Figures 5, 6, 7.
# Also creates smaller case-control samples for the
# benchmarking in Table 2.

# Case-control sampling for the benchmarking section (part of Table 2).
# Creates a case-control sample, where proportion of cases is multiplied
# by 50.  Creates sample size of n=1500.


#===================================================================
# User inputs.
#===================================================================
source('./code/helper-functions/stillbirth-testtrain-split.R')

# Filepath for stillbirth data.
# Uses cleaned data from output script 05_01-stillbirth-process-data.R
perinatal_source = paste0(getwd(),'/Data/stillbirth-clean.csv')
       
div_prop_by = 50

# Output paths/names for each case-control sample.
outfile_cc  = paste0(getwd(),'/data/perinatal-mortality/stillbirth-data-casecontrol-50-for-paper.csv')

#===================================================================
# Workspace prep
#===================================================================

# Install fastbart package from local directory.
library(devtools)
library(data.table)
library(tidyverse)


#===================================================================
# Read data
#===================================================================

# Load functions for re-leveling perinatal data and for test-train split.
source(paste(getwd(),'code/00_stillbirth-data-processing/stillbirth-dataformat-functions.R',sep='/'))

# Read both data sets and adjust factor levels.
data = fix_levels(as.data.frame(fread(perinatal_source, sep=',', header=T, stringsAsFactors = T)))

####################################################################
# Calculate case-control sample sizes.
####################################################################

#===================================================================
# For Figures 5-6-7.
#===================================================================

# Note: Use all fd cases, and adjust sample sizes for each gest age to reflect proportions.
ga_tab = as.data.frame(data %>% group_by(ga) %>% summarize('n' = n(), 'fd' = sum(fd==1)))
ga_tab$p = ga_tab$fd / ga_tab$n
ga_tab$p_adj = ga_tab$p * div_prop_by
ga_tab$n_adj = ga_tab$fd / ga_tab$p_adj

# Set sample sizes for both fd=1 and fd=0 for each gestational age.
ga_tab$fd1 = ga_tab$fd
ga_tab$fd0 = round(ga_tab$n_adj - ga_tab$fd)

ga_tab$divide_prop_by = div_prop_by

####################################################################
# Case-control sampling.
####################################################################


# Create vector of indices to keep.
keeps = NULL

# Loop through each gestational age, and sample fd=1 and fd=0 cases in appropriate size.
for(i in 1:nrow(ga_tab)){
   temp_ga = ga_tab$ga[i]
   n1 = ga_tab$fd1[i]
   n0 = ga_tab$fd0[i]
   
   keeps1 = sample(data %>% filter(ga==temp_ga, fd==1) %>% collect %>% .[["id"]], size=n1, replace=F)
   keeps0 = sample(data %>% filter(ga==temp_ga, fd==0) %>% collect %>% .[["id"]], size=n0, replace=F)
   
   keeps = c(keeps, keeps1, keeps0)
}

data_fd_cc = data %>% filter(id %in% keeps)

# Clean up variables - get rid of weeks, switch to a variable that represents baseline ga (34).
data_fd_cc = data_fd_cc %>% dplyr::rename(gest_age34 = weeks, gest_age = ga)


#===================================================================
# Add test/train indicators to files.
#===================================================================

# Stillbirth.
fd_split = testtrain(data_fd_cc, 'fd', 0.3)

data_fd_cc$testtrain = 'train'
data_fd_cc$testtrain[which(data_fd_cc$id %in% fd_split$test_id)] = 'test'


#===================================================================
# Write output to csv files.
#===================================================================

# Stillbrith and stillbirth expanded.
fwrite(data_fd_cc,outfile_cc)


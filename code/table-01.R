#######################################################################################
##  DESCRIPTION:
##
## This script summarizes information about the dataset, and creates the Table 1 to
## display cohort characteristics.  Note that the full dataset is not publicly available;
## this script is an example, but will not recreate Table 1.
#######################################################################################

###########################################################################
# WORKSPACE PREP
###########################################################################

# Clean workspace.
rm(list=ls())

library(data.table)  # For fast csv read/write.
library(tidyverse)
library(xtable)

# ggplot settings
theme_set(theme_bw(base_size=13, base_family='Helvetica'))

# Set output directories.
out.fig = paste0(getwd(), '/output-figures')
out.csv = paste0(getwd(), '/output-files/')

# Load my custom functions.
source(paste0(getwd(),'/code/helper-functions/create-table1.R'))

###########################################################################
# MISSING DATA AND NUMBERS OF RECORDS
###########################################################################

# These values come from the following unix command to count lines in a file.
# find . -type f -exec wc -l {} \;
nraw = 5476562-1
nclean = 4553869-1

# Number of records removed.
nraw - nclean


###########################################################################
# TABLE 1
###########################################################################

df = as.data.frame(fread(paste0(getwd(),'/data/perinatal-mortality/stillbirth-data-example.csv')))
colvar = 'hirisk'

table1 = cohortchars(df,colvar)

print(xtable(table1),
      include.rownames=FALSE,
      hline.after = c(0,1,6,7,11,12,14,15,17,18))


###########################################################################
library(tidyverse)

# How many stillbirth records in dataset? Prevalence per 1000 births.
df %>% summarise('total_sb' = sum(fd), 'total_b' = n(), 'prev per 1000' = round(sum(fd) / n() * 1000,2))

# How many stillbirth records in dataset for hr and lr groups? Prevalence per 1000 births.
df %>% filter(hirisk==1) %>% summarise('prev per 1000' = round(sum(fd) / n() * 1000,2))
df %>% filter(hirisk==0) %>% summarise('prev per 1000' = round(sum(fd) / n() * 1000,2))

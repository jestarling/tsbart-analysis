# Create Table 1.
Rscript ./code/table-01.R

# Create Table 2.
Rscript ./code/table-02.R

# Create Figure 2 and Figure 6.
Rscript ./code/figure-02-and-figure06.R

# Create Figure 3.
Rscript ./code/figure-03.R

# Create Figures 4 and Table 4.
Rscript ./code/figure-04-and-table-04.R

# Create Figure 5.
Rscript ./code/figure-05.R


### For the first simulation study (Table 2), 
### this work was done in parallel on the TACC Stampede2
### supercomputer.  The following scripts create one dataset
### and illustrate model fit and calculations.

# To create one dataset:
Rscript ./code/table-02-a_generate-data.R

# To illustrate coverage analysis:
Rscript ./code/table-02-b_one-example-run.R





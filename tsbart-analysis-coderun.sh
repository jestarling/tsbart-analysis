# Create Table 1.
Rscript ./code/table-01.R

# Create Figure 2.
Rscript ./code/figure-02.R

# Create Figure 3 and Figure 6.
Rscript ./code/figure-03-and-figure06.R

# Create Figure 4.
Rscript ./code/figure-04.R

# Create Figures 5 and Table 4.
Rscript ./code/figure-05-and-table-04.R

# Create Table 3.
Rscript ./code/table-03.R

### For the first simulation study (Table 2), 
### this work was done in parallel on the TACC Stampede2
### supercomputer.  The following scripts create one dataset
### and illustrate model fit and calculations.

# To create one dataset:
Rscript ./code/table-02-a_generate-data.R

# To illustrate coverage analysis:
Rscript ./code/table-02-b_one-example-run.R





# This script shows the calculations for coverage and MSE
# using a single dataset and one scenario, as an illustration.

# Requires the fastbart package; available on request.

# Set filename.
file = "df_linear_interaction.csv"
ntree=200; nburn=500; nsim=5000

###########################################################################
# WORKSPACE PREP
###########################################################################

# Install and load funbart package.
library(devtools)
install_github("jestarling/tsbart")
library(tsbart)
library(dbarts)
library(splines)
library(data.table)  # For fast csv read/write.
library(tidyverse)   # For tidy data/plotting.
library(mgcv)        # For penalized splines

# Load helper functions.
source('./code/helper-functions/stillbirth-functions-modelfitutils.R')
source('./code/helper-functions/stillbirth-functions-cvutils.R')
source('./code/helper-functions/stillbirth-functions-testtrain-split.R')
source('./code/helper-functions/stillbirth-functions-binll.R')

###########################################################################
# DATA LOAD AND PREP
###########################################################################

data = as.data.frame(fread(paste0('./data/simulations/', file)))

# 2. Set the column names of covariates to include in the models.
xcols = paste0('x',1:10)

# 3. Test-train splits.
train = na.omit(data %>% filter(testtrain=='train'))
test = data %>% filter(testtrain=='test')


###########################################################################
# Fit models.
###########################################################################
   
# tsbart with tuned hyperparameter.
cv_tsb = tsbart_fit_util(train, test, xcols, yvar='y',tvar='t', ec=1.2, ntree=ntree, nburn=nburn, nsim=nsim, init_upper=1.65, init_lower=-1.65)

# tsbart with default hyperparameter 1
cv_tsb_default = tsbart_fit_util(train, test, xcols, yvar='y',tvar='t', ec=1, ntree=ntree, nburn=nburn, nsim=nsim, init_upper=1.65, init_lower=-1.65)

# vanilla bart
cv_vb = bart_fit_util(train, test, yvar='y', tvar='t', xcols=c(xcols,'t'), m, burn, nsim, init_upper=1.65, init_lower=-1.65)

# Fit splines model with linear interaction.
cv_sp1 = spline_fit_util_sim(train, test)

# Fit splines model with basis interaction.
cv_sp2 = spline_fit_util_sim(train, test, spline_interactions=TRUE)

# Fit p-splines model.
cv_pensp = penspline_fit_util_sim(train, test)

# Save test fits.
test = extractFits(test, cv_tsb, cv_vb, cv_sp1, cv_sp2, cv_pensp, cv_tsb_default, adjust=F)


###########################################################################
# Assess coverage.
###########################################################################

# Calculate whether each observation in each iteration is covered.
test$covg_tsb = ifelse((test$haz_fun >= test$phat_oos_tsb_lb) &
                          (test$haz_fun <= test$phat_oos_tsb_ub), 1, 0)

test$covg_tsb_default = ifelse((test$haz_fun >= test$phat_oos_tsb_default_lb) &
                                  (test$haz_fun <= test$phat_oos_tsb_default_ub), 1, 0)

test$covg_vb = ifelse((test$haz_fun >= test$phat_oos_vb_lb) &
                         (test$haz_fun <= test$phat_oos_vb_ub), 1, 0)

test$covg_sp1 = ifelse((test$haz_fun >= test$phat_oos_sp1_lb) &
                          (test$haz_fun <= test$phat_oos_sp1_ub), 1, 0)

test$covg_sp2 = ifelse((test$haz_fun >= test$phat_oos_sp2_lb) &
                          (test$haz_fun <= test$phat_oos_sp2_ub), 1, 0)

test$covg_pensp = ifelse((test$haz_fun >= test$phat_oos_pensp_lb) &
                            (test$haz_fun <= test$phat_oos_pensp_ub), 1, 0)

# Average coverage for each method.
covg = as.data.frame(matrix(0, nrow=6, ncol=2))
colnames(covg) = c('method','coverage')
covg$method = c('tsb','tsb_default','vb','sp1','sp2','pensp')

for(i in 1:length(covg$method)){
   method = covg$method[i]
   
   col_idx = grepl(paste0("covg_",method), colnames(test))
   col_idx[which(colnames(test) %in% c('iteration','haz_fun'))] = TRUE
   temp = test[,col_idx]
   
   varname = paste0("covg_",method)
   covg$coverage[i] = mean(temp[,2])
}

covg$scenario = gsub("_test","",file)

###########################################################################
# MSE for each method.
###########################################################################

covg$mse = NULL

for(i in 1:length(covg$method)){
   method = covg$method[i]
   
   col_idx = grepl(paste0("phat_oos_",method), colnames(test))
   col_idx[which(colnames(test) %in% c('iteration','haz_fun'))] = TRUE
   temp = test[,col_idx]
   
   temp$sqerr = (temp[,1] - temp[,2])^2
   covg$mse[i] = round(mean(temp$sqerr),6)
}

covg

write.csv('./output-files/table-02-example.csv',row.names=F)

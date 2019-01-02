#######################################################################################
##  DESCRIPTION:
##
##  This script analyzes the stillbirth dataset, comparing the following models:
##    1. tsbART
##    2. vanilla BART (using sparapani method for correct conditional probabilities)
##    3. conditional logistic regression with splines (preliminary analysis)
##    4. treed gaussian processes (TGP), with only time included in trees
#
##  The script does the following:
##    1. Fits each model to the training dataset.
##    2. Calculates out-of-sample logloss on the test dataset.
##    3. Calculates estimates for a patient panel, for plotting purposes.
# 
##  Steps 1-2 are done k=5 five times for each dataset, to get cross-validated 
##  out-of-sample log-loss.
#######################################################################################

###########################################################################
# WORKSPACE PREP
###########################################################################

rm(list=ls())

# Set trees/burns/draws for tsbart and bart.
ntree=200; nburn=1000; nsim=10000
ntree=20; nburn=20; nsim=20
#-------------------------------------------------------------------------
# Install tsbart and fastbart packages.
#-------------------------------------------------------------------------

library(devtools)    # For installing packages via github.
install_github('jestarling/tsbart')
library(tsbart)

# Install fastbart package from local directory. (Vanilla BART first, for propensity scores).
install.packages('./fastbart-path/fastbart_1.0.tar.gz', type='source', repos=NULL)
library(fastbart)

# Install packages for splines and TGP.
library(splines)

# Load other libraries.
library(data.table)  # For fast csv read/write.
library(gridExtra)   # For plotting multi-pane ggplots.
library(tidyverse)   # For tidy data/plotting.
library(mgcv)        # For penalized splines
library(magrittr)
library(cowplot)
library(ggthemes)
library(viridis)

# ggplot settings
theme_set(theme_bw(base_size=13, base_family='Helvetica'))

# Set output directories.
out.fig = './output-figures/'
out.csv = './output-files/'

# Load my custom functions.
source('./code/helper-functions/stillbirth-functions-modelfitutils.R')
source('./code/helper-functions/stillbirth-functions-cvutils.R')
source('./code/helper-functions/stillbirth-functions-testtrain-split.R')
source('./code/helper-functions/stillbirth-functions-binll.R')

###########################################################################
# DATA LOAD AND PREP
###########################################################################

#-------------------------------------------------------------------
# 1. Test/train data from case-control sampled obstetrics dataset.

# Read csv file and expand data frame for correct conditional probabilities.
data = as.data.frame(fread(paste0(getwd(),'/data/perinatal-mortality/stillbirth-data-casecontrol-50-for-paper.csv')))

# Drop the testtrain column already constructed.  we will draw new test/train splits for each cv fold.
data = tsbart::survPrep(data, 'gest_age34', 'fd'); data$gest_age = data$gest_age34+33

#-------------------------------------------------------------------
# 2. Set the column names of covariates to include in the models.

xcols = c('diabetes','htn', 'otherRisk','momAge','momEthnicity','multiparous','induce','male','wtgainQ','birthwtQ')

#-------------------------------------------------------------------
# 3. Test-train split.
temp = testtrain(data, test_pct=.2)
test = temp$test; test$testtrain = 'test'
train = temp$train; train$testtrain = 'train'

#-------------------------------------------------------------------
# 3. Read in new patient panel.
newpts = read.csv(paste0(getwd(), '/data/perinatal-mortality/perinatal-ptpanel.csv'),stringsAsFactors=F)
newpts$id = newpts$id + 5000000  #Panel IDs start at 5000001.
ptpanel = tsbart::makePredGrid(newpts, 'gest_age34', sort(unique(data$gest_age34)))

# Add 'fd' as placeholder, even though we don't know it for patient panel.
# This is so that we can join the test and panel dataframes together.
ptpanel$fd=0
ptpanel$testtrain = 'panel'

# Add patient panel to end of training data.
#    (Want models to predict on both the test data and the new patient panel,
#    for plotting and OOS log-loss.)

test = dplyr::bind_rows(ptpanel[c(xcols, 'id', 'gest_age34','fd','testtrain')],
                        test[c(xcols, 'id', 'gest_age34','fd','testtrain')])

test$testtrain = c(rep('panel', nrow(ptpanel)), rep('test',nrow(test) - nrow(ptpanel)))


###########################################################################
# Fit models.
###########################################################################

# tsbart with previously selected tuned hyperparameter of .1
cv_tsb = tsbart_fit_util(train, test, xcols, ec=.1, ntree=ntree, nburn=nburn, nsim=nsim)

# tsbart with default hyperparameter 1
cv_tsb_default = tsbart_fit_util(train, test, xcols, ec=1, ntree=ntree, nburn=nburn, nsim=nsim)

# vanilla bart
cv_vb = bart_fit_util(train, test, xcols=c(xcols,'gest_age34'), ntree=ntree, nburn=nburn, nsim=nsim)

# Fit splines model with linear interaction.
cv_sp1 = spline_fit_util(train, test)

# Fit splines model with basis interaction.
cv_sp2 = spline_fit_util(train, test, spline_interactions=TRUE)

# Fit p-splines model.
cv_pensp = penspline_fit_util(train, test)

# Extract fits.
test = extractFits(test, cv_tsb, cv_vb, cv_sp1, cv_sp2, cv_pensp, cv_tsb_default, adjust=T)
ptpanel = test %>% filter(testtrain=='panel')


###########################################################################
# Plotting
###########################################################################

ptpanel$gest_age = ptpanel$gest_age34 + 33
y_limits = c(0,6)
ids_newplt = c(5000001, 5000002, 5000004, 5000027, 5000032)

row1 = panplt_tsb_combos = ggplot((ptpanel %>% filter(id %in% ids_newplt)), aes(x=gest_age, y=phat_oos_tsb), colour='black') +
   geom_line(aes(y=phat_oos_vb), colour='grey70') +
   geom_line(aes(y=phat_oos_sp1), colour='grey70') +
   geom_line(aes(y=phat_oos_sp2), colour='grey70') +
   geom_line(aes(y=phat_oos_pensp), colour='grey70') +
   geom_ribbon(aes(x=gest_age, ymin=phat_oos_tsb_lb, ymax=phat_oos_tsb_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line() +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   coord_cartesian(ylim=y_limits) +
   labs(x='',y='tsBART') +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines")) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') 

row2 = ggplot((ptpanel %>% filter(id %in% ids_newplt)), aes(x=gest_age, y=phat_oos_vb), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey70') +
   geom_line(aes(y=phat_oos_sp1), colour='grey70') +
   geom_line(aes(y=phat_oos_sp2), colour='grey70') +
   geom_line(aes(y=phat_oos_pensp), colour='grey70') +
   geom_ribbon(aes(x=gest_age, ymin=phat_oos_vb_lb, ymax=phat_oos_vb_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line() +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='BART') +
   coord_cartesian(ylim=y_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines")) +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85')

row3 = ggplot((ptpanel %>% filter(id %in% ids_newplt)), aes(x=gest_age, y=phat_oos_sp1), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey70') +
   geom_line(aes(y=phat_oos_vb), colour='grey70') +
   geom_line(aes(y=phat_oos_sp2), colour='grey70') +
   geom_line(aes(y=phat_oos_pensp), colour='grey70') +
   geom_ribbon(aes(x=gest_age, ymin=phat_oos_sp1_lb, ymax=phat_oos_sp1_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line() +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='Splines 1') +
   coord_cartesian(ylim=y_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines"))  +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85')

row4 = ggplot((ptpanel %>% filter(id %in% ids_newplt)), aes(x=gest_age, y=phat_oos_sp2), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey70') +
   geom_line(aes(y=phat_oos_vb), colour='grey70') +
   geom_line(aes(y=phat_oos_sp1), colour='grey70') +
   geom_line(aes(y=phat_oos_pensp), colour='grey70') +
   geom_ribbon(aes(x=gest_age, ymin=phat_oos_sp2_lb, ymax=phat_oos_sp2_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line() +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='Splines 2') +
   coord_cartesian(ylim=y_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines"))  +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85')

row5 = ggplot((ptpanel %>% filter(id %in% ids_newplt)), aes(x=gest_age, y=phat_oos_pensp), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey70') +
   geom_line(aes(y=phat_oos_vb), colour='grey70') +
   geom_line(aes(y=phat_oos_sp1), colour='grey70') +
   geom_line(aes(y=phat_oos_sp2), colour='grey70') +
   geom_ribbon(aes(x=gest_age, ymin=phat_oos_pensp_lb, ymax=phat_oos_pensp_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line() +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='P-splines') +
   coord_cartesian(ylim=y_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines"))  +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85')


library(grid)
grid.newpage()
all = grid.arrange(row1, row2, row3, row4, row5, 
             ncol=1,
             bottom='Gestational age (Wks)',
             left='Risk of stillbirth per 1,000 remaining pregnancies',
             top='       Patient 1                   Patient 2                   Patient 3                 Patient 4                Patient 5')

ggsave(paste0(out.fig,'figure-05.pdf'), all,
       width=8, height=8, units='in', dpi=300, limitsize=TRUE)

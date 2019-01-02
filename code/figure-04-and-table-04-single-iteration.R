###########################################################################
# WORKSPACE PREP
###########################################################################

rm(list=ls())

# Set trees/burns/draws for tsbart and bart.
ntree=200; nburn=1000; nsim=10000  

#-------------------------------------------------------------------------
# Install tsbart and fastbart packages.
#-------------------------------------------------------------------------

library(devtools)    # For installing packages via github.
library(tsbart)
library(fastbart)
library(splines)

# Load other libraries.
library(data.table)  # For fast csv read/write.
library(gridExtra)   # For plotting multi-pane ggplots.
library(tidyverse)   # For tidy data/plotting.
library(mgcv)        # For penalized splines
library(viridis)
library(ggthemes)

# ggplot settings
theme_set(theme_bw(base_size=13, base_family='Helvetica'))

# Load my custom functions.
source('./code/helper-functions/stillbirth-functions-modelfitutils.R')
source('./code/helper-functions/stillbirth-functions-cvutils.R')
source('./code/helper-functions/stillbirth-functions-testtrain-split.R')
source('./code/helper-functions/stillbirth-functions-binll.R')
source('./code/helper-functions/ggtheme-publication.R')

###########################################################################
# DATA LOAD AND PREP
###########################################################################

#-------------------------------------------------------------------
# 1. Test/train data from case-control sampled obstetrics dataset.

# Read csv file and expand data frame for correct conditional probabilities.
data = as.data.frame(fread(paste0(getwd(),'/data/perinatal-mortality/stillbirth-data-casecontrol-50-for-paper.csv')))
data = upsample(data, n=1000, sb_pct = .5)
data = tsbart::survPrep(data, 'gest_age34', 'fd'); data$gest_age = data$gest_age34+33

#-------------------------------------------------------------------
# 2. Set the column names of covariates to include in the models.
xcols = c('diabetes','htn', 'otherRisk','momAge','momEthnicity','multiparous','induce','male','wtgainQ','birthwtQ')

#-------------------------------------------------------------------
# 3. Test-train split.
temp = testtrain(data, test_pct=.2)
test = temp$test; test$testtrain = 'test'
train = temp$train; train$testtrain = 'train'

###########################################################################
# Fit models.
###########################################################################

# tsbart with previously tuned hyperparameter of .1.
cv_tsb = tsbart_fit_util(train, test, xcols, ec=5, ntree=ntree, nburn=nburn, nsim=nsim, init_upper=.01, init_lower=qnorm(.05))

# tsbart with default hyperparameter 1
cv_tsb_default = tsbart_fit_util(train, test, xcols, ec=1, ntree=ntree, nburn=nburn, nsim=nsim, init_upper=.01, init_lower=qnorm(.05))

# vanilla bart
cv_vb = bart_fit_util(train, test, xcols=c(xcols,'gest_age34'), ntree=ntree, nburn=nburn, nsim=nsim, init_upper=.01, init_lower=qnorm(.05))

# Fit splines model with linear interaction.
cv_sp1 = spline_fit_util(train, test)

# Fit splines model with basis interaction.
cv_sp2 = spline_fit_util(train, test, spline_interactions=TRUE)

# Fit p-splines model.
cv_pensp = penspline_fit_util(train, test)

# Add out of sample fit info to dataframe.
test = extractFits(test, cv_tsb, cv_vb, cv_sp1, cv_sp2, cv_pensp, cv_tsb_default)

###########################################################################
# Binomial log-likelihoods (per person) & Weekly misclassification error.
###########################################################################
thrsh = .5

binll_oos = binll_per_person(test, method='tsb', thresh=thrsh)

for(meth in c('tsb_default','vb','sp1','sp2','pensp')){
   binll_oos = rbind.data.frame(binll_oos, binll_per_person(test, meth, thresh=thrsh))
}

# Display results.
binll_oos %>% filter(gest_age34 =='all') %>% select(method,logl)
binll_oos %>% filter(gest_age34 !='all')


# Separate logl into its own table.
bll_tbl = binll_oos %>% 
   filter(gest_age34 != "all") %>%
   select(method, gest_age34, logl) %>%
   spread(key=gest_age34, value=logl)
colnames(bll_tbl)[2:10] = 34:42

write.csv(binll_oos %>% filter(gest_age34 =='all') %>% select(method,logl), './output-files/table-04.csv', row.names=F)

#--------------------------------------------------------------------------
# Save and plot results.
#--------------------------------------------------------------------------

binll_rel = binll_oos %>% filter(gest_age34 !='all')
binll_rel$bll_relto_tsb = NA

for(m in unique(binll_rel$method)){
   binll_rel$bll_relto_tsb[which(binll_rel$method==m)] = binll_rel$logl[which(binll_rel$method==m)] / binll_rel$logl[which(binll_rel$method=='tsb')]
}

ll_plt = ggplot(binll_rel, aes(x=gest_age, y=bll_relto_tsb, colour=method)) + 
   geom_line(size=.75) + 
   labs(x='Gestational age (Wks)', y='Out of sample logl relative to tsBART') +
   scale_colour_Publication()+ theme_Publication()

ggsave('./output-figures/figure-04.pdf', ll_plt,
       width=8, height=8, units='in', dpi=300, limitsize=TRUE)


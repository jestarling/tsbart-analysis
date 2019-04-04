###########################################################################
# WORKSPACE PREP
###########################################################################

rm(list=ls())

# Average over k reps.
kreps=5

# Set trees/burns/draws for tsbart and bart.
ntree=200; nburn=1000; nsim=10000  

#-------------------------------------------------------------------------
# Workspace setup.
#-------------------------------------------------------------------------

library(devtools)    # For installing packages via github.
library(tsbart)
library(dbarts)
library(splines)

# Load other libraries.
library(data.table)  # For fast csv read/write.
library(gridExtra)   # For plotting multi-pane ggplots.
library(tidyverse)   # For tidy data/plotting.
library(mgcv)        # For penalized splines
library(viridis)
library(ggthemes)

# ggplot settings
theme_set(theme_bw(base_size=16, base_family='Helvetica'))

# Load my custom functions.
source('./code/helper-functions/stillbirth-functions-modelfitutils.R')
source('./code/helper-functions/stillbirth-functions-cvutils.R')
source('./code/helper-functions/stillbirth-functions-testtrain-split.R')
source('./code/helper-functions/stillbirth-functions-binll.R')
source('./code/helper-functions/ggtheme-publication.R')

###########################################################################
# DATA LOAD AND PREP
###########################################################################

for(k in 1:kreps){
   
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
   
   if(k==1){
      binll_oos = binll_per_person(test, method='tsb', thresh=thrsh)
      
      for(meth in c('tsb_default','vb','sp1','sp2','pensp')){
         binll_oos = rbind.data.frame(binll_oos, binll_per_person(test, meth, thresh=thrsh))
      }
   } else{
      
      binll_oos_k = binll_per_person(test, method='tsb', thresh=thrsh)
      
      for(meth in c('tsb_default','vb','sp1','sp2','pensp')){
         binll_oos_k = rbind.data.frame(binll_oos, binll_per_person(test, meth, thresh=thrsh))
      }
      
      binll_oos[,3:6] = binll_oos[,3:6] + binll_oos_k[,3:6] 
   }
}

# Average over iterations.
binll_oos[,3:6] = binll_oos[,3:6]/kreps

# Separate logl into its own table.
bll_tbl = binll_oos %>% 
   filter(gest_age34 != "all") %>%
   select(method, gest_age34, logl) %>%
   spread(key=gest_age34, value=logl)
colnames(bll_tbl)[2:10] = 34:42

tbl_4 = binll_oos %>% filter(gest_age34 =='all') %>% select(method,logl)
write.csv(tbl_4, './output-files/table-04.csv', row.names=F)

#--------------------------------------------------------------------------
# Save and plot results.
#--------------------------------------------------------------------------

### -----------
# Using saved file.
binll_oos = read.csv('./output-files/binll_oos_upsampled_binll-by-week.csv')
colnames(binll_oos)[-1] = 34:42
binll_oos = binll_oos %>% gather('gest_age','logl','34':'42')
binll_oos$gest_age = as.numeric(binll_oos$gest_age)
binll_oos$gest_age34 = binll_oos$gest_age - 33
### -----------

binll_rel = binll_oos %>% filter(gest_age34 !='all')
binll_rel$bll_relto_tsb = NA

binll_rel$labels = ifelse(binll_rel$method=="pensp", "P-splines",
                          ifelse(binll_rel$method=="sp1", "Splines 1",
                                 ifelse(binll_rel$method=="sp2", "Splines 2",
                                        ifelse(binll_rel$method=="tsb", "tsBART (tuned)",
                                               ifelse(binll_rel$method=="tsb_default", "tsBART (default)",
                                                      "BART")))))

for(m in unique(binll_rel$method)){
   binll_rel$bll_relto_tsb[which(binll_rel$method==m)] = binll_rel$logl[which(binll_rel$method==m)] / binll_rel$logl[which(binll_rel$method=='tsb')]
}

ll_plt = ggplot(binll_rel, aes(x=gest_age, y=bll_relto_tsb, colour=labels, linetype=labels)) + 
   geom_line(size=.75) + 
   labs(x='Gestational age (Wks)', y='Out of sample relative log-loss') +
   scale_colour_colorblind(name="") + scale_linetype_manual(name='', values=c(6,5,3,2,4,1)) + theme_Publication() +
   theme(legend.key.size=unit(.75,"cm"))

ll_plt

ggsave('./output-figures/figure-04.pdf', ll_plt,
      width=5, height=5, units='in', dpi=300, limitsize=TRUE)

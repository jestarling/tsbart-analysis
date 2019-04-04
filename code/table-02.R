# Creates data used in Table 1's simulation study.

#===================================================================
# User inputs.
#===================================================================

# fastbart filepath. Example is:
fastbart_path = './fastbart'

#===================================================================
# Workspace prep.
#===================================================================

library(data.table) # For fast csv read.
library(tidyverse)
library(devtools)
install_github('jestarling/tsbart')
library(tsbart)

source(paste0(getwd(), '/code/helper-functions/create-table4.R'))

#===================================================================
# Read all file names.
#===================================================================

files <- list.files(paste0(getwd(), "/data/supp-simulations-table-04/"))

# Define function to compare each dataset's model fits.
tsbartCompare = function(filepath, filename, probit=FALSE, out.csv){
   
   require(tidyverse)
   require(data.table)
   require(tsbart)
   require(dbarts)
   
   # Check that csv.out exists before proceeding, if it does not already.
   dir.create(out.csv)
   
   #===================================================================
   # Read data.
   #===================================================================
   
   sim = fread(paste0(filepath, filename))
   
   # Set up train and test data sets.
   train = sim %>% filter(train=='train')
   test = sim %>% filter(train=='test')
   
   # Extract values.
   xx = train[,which( substr(colnames(train),1,1) == "x")]
   x_pred =  test[,which( substr(colnames(test),1,1) == "x")]
   
   ti = train$ti;     ti_pred = test$ti      # Time points for each obs.
   fx = train$fx;     fx_pred = test$fx      # True underlying function, without epsilon noise.
   y = train$y;       y_pred = test$y        # Vector of responses.
   
   yobs=NULL; yobs_pred=NULL; offset=NULL
   
   if(probit){
      yobs = train$yobs;  yobs_pred = test$yobs
      offset = offset$train[1]
   } 
   
   #############################################################################################
   ###   1. Vanilla BART.
   #############################################################################################
   
   # Install fastbart package from local directory.
   library(dbarts)
   
   # Calibrate BART's error variance a la CGM 2010 (method 2)
   df = data.frame(xx, ti, y)
   lmf = lm(y~.,df)
   sighat = sigma(lmf) #sd(y) #summary(lmf)$sigma
   
   # Hyperparameters
   nu = 3
   sigq = .9
   qchi = qchisq(1.0-sigq,nu)
   lambda = (sighat*sighat*qchi)/nu
   
   # Decision rule cutpoints for each predictor.  
   # Note: Not using makeCutpoints() function here, since each xj has one fixed value in this sim.
   # Need to set up grids individually to cover all possible values for each xj.
   cutpoints = makeCutpoints(cbind.data.frame(xx,ti))
   
   # Fit vanilla BART model.
   m = 200; burn = 200; nd = 1000
   # fit_v = bartRcppClean(y_ = y, x_ = t(cbind.data.frame(xx,ti)), # obsvervations must be in *columns*
   #                       xpred_ = t(cbind.data.frame(x_pred,'ti' = ti_pred)), #Dpred,#[,-2,drop=F],
   #                       xinfo_list = cutpoints,
   #                       burn, nd, m,
   #                       lambda, nu, kfac=2,
   #                       paste0(getwd(),"/trees.txt"),
   #                       RJ=FALSE)
   
   xcols = c(xcols,which(colnames(train)=='ti'))
   fit_v = bart(x.train=train[xcols], y.train=train$y, x.test=test[xcols],ntree=ntree,ndpost=nsim,nskip=nburn)
   
   print('Vanilla BART completed.')
   
   #############################################################################################
   ###   1. tsBART.
   #############################################################################################
   
   # Remove tsBART package.
   detach("package:tsbart", unload=TRUE)
   
   # Install and load tsBART package.
   library(tsbart)
   
   #===================================================================
   # BART setup for model fitting.  (Hyperparameters, etc.)
   #===================================================================
   
   # Calibrate BART's error variance a la CGM 2010 (method 2)
   df = data.frame(xx, y)
   lmf = lm(y~.,df)
   sighat = sigma(lmf) #sd(y) #summary(lmf)$sigma
   
   # Hyperparameters
   nu = 3
   sigq = .9
   qchi = qchisq(1.0-sigq,nu)
   lambda = (sighat*sighat*qchi)/nu
   
   # Decision rule cutpoints for each predictor.  
   # Note: Not using makeCutpoints() function here, since each xj has one fixed value in this sim.
   # Need to set up grids individually to cover all possible values for each xj.
   cutpoints = makeCutpoints(xx)
   
   #=====================================================================
   #=== Evaluate optimal expected number of crossings by selecting
   #=== one with lowest WAIC.
   #=====================================================================
   
   # Candidate values for expected number of crossings.
   ecross_candidates = seq(.5,8,by=.5)
   
   head(ti_pred)
   
   # Evaluate optimal number of crossings.
   if(probit){
      ecrossTune = tuneEcross(ecross_candidates,
                              m=200, nd=50, burn=50,
                              y_ = y, x_ = t(xx) , xpred_ = t(x_pred),# obsvervations must be in *columns*
                              t_=ti , tpred_ =ti_pred, xinfo_list=cutpoints, 
                              lambda, nu, treef_name_ = paste0(getwd(),'trees.txt'), offset=offset, yobs=yobs, probit=TRUE)
   } else{
      ecrossTune = tuneEcross(ecross_candidates,
                              m=200, nd=50, burn=50,
                              y_ = y, x_ = t(xx) , xpred_ = t(x_pred),# obsvervations must be in *columns*
                              t_=ti , tpred_ =ti_pred, xinfo_list=cutpoints, 
                              lambda, nu, treef_name_ = paste0(getwd(),'trees.txt'))
   }
   
   exp_cross = ecrossTune$ecross_opt
   
   print('Ecross tune completed.')
   
   #===================================================================
   # Fit tsBART model.
   #===================================================================
   
   m = 200; burn = 200; nd = 1000
   
   if(probit){
      
      # Calculate offset
      phat = mean(yobs)
      offset = qnorm(phat)
      
      fit = tsbartProbit(y_ = y,
                         yobs_ = yobs,
                         x_ = t(xx), # obsvervations must be in *columns*
                         t_ = ti, # new input for time points; default 1.
                         xpred_ = t(x_pred), #Dpred,#[,-2,drop=F],
                         tpred_ = ti_pred,
                         xinfo_list = cutpoints,
                         burn = burn, nd = nd, m = m,
                         lambda=lambda, nu=nu, ecross = exp_cross,offset=offset,
                         treef_name_ =  paste0(getwd(),'_trees.txt'))
      
   } else{
      fit = tsbartFit(y_ = y, x_ = t(xx), # obsvervations must be in *columns*
                      t_ = ti, # new input for time points; default 1.
                      xpred_ = t(x_pred), #Dpred,#[,-2,drop=F],
                      tpred_ = ti_pred,
                      xinfo_list = cutpoints,
                      burn = burn, nd = nd, m = m,
                      lambda=lambda, nu=nu, ecross = exp_cross,
                      treef_name_ =  paste0(getwd(),'_trees.txt'))
   }
   
   print('tsBART completed.')
   
   #############################################################################################
   ###   3. Comparison between tsbart and vanilla.
   #############################################################################################
   
   #===================================================================
   # In-sample
   #===================================================================
   
   # tsBART
   check_is = checkFit(y=y, 
                       mcmcdraws = fit$mcmcdraws,
                       sig = fit$sigma,
                       probit,
                       doWaic=TRUE,
                       yobs)
   
   # Vanilla BART
   check_is_v = checkFit(y=y, 
                         mcmcdraws = fit_v$yhat.train,
                         sig = fit_v$sigma,
                         probit,
                         doWaic=TRUE,
                         yobs)
   
   #===================================================================
   # Out-of-sample
   #===================================================================
   
   # Functional BART
   check_oos = checkFit(y = y_pred,
                        mcmcdraws = fit$mcmcdraws_oos,
                        sig = fit$sigma,
                        probit,
                        doWaic=TRUE,
                        yobs)
   
   # Vanilla BART
   check_oos_v = checkFit(y = y_pred,
                          mcmcdraws = fit_v$yhat.test,
                          sig = fit$sigma,
                          probit,
                          doWaic=TRUE,
                          yobs=yobs)
   
   #===================================================================
   # Assemble data frame for comparisons.
   #===================================================================
   
   # Set up ggplot data frame.
   log_dens = c(check_is$logdensity_mcmc, check_is_v$logdensity_mcmc, check_oos$logdensity_mcmc, check_oos_v$logdensity_mcmc)
   data_set = c(rep('in-sample',nd*2), rep('out-of-sample',nd*2))
   model = c(rep('functional',nd), rep('vanilla',nd), rep('functional',nd), rep('vanilla',nd))
   lldf = cbind.data.frame('logdensity' = log_dens, 'dataset' = data_set, 'model' = model)
   
   #############################################################################################
   ###   4. Arrange and output data frame.
   #############################################################################################
   
   #===================================================================
   # Begin with training data set.
   #===================================================================
   
   df_tr = train
   
   df_tr$pmean_f = fit$pmean
   df_tr$pmean_v = colMeans(fit_v$postfit)
   
   # Lower and upper bounds.
   df_tr$lb_f = apply(fit$mcmcdraws, 2, function(x) quantile(x,.025))
   df_tr$ub_f = apply(fit$mcmcdraws, 2, function(x) quantile(x,.975))
   
   df_tr$lb_v = apply(fit_v$postfit, 2,function(x) quantile(x,.025))
   df_tr$ub_v = apply(fit_v$postfit, 2,function(x) quantile(x,.975))
   
   # Full MCMC.
   df_tr = cbind.data.frame(df_tr, 'funMCMC'=as.data.frame(t(fit$mcmcdraws)))
   df_tr = cbind.data.frame(df_tr, 'vMCMCM'=t(fit_v$postfit))
   
   #===================================================================
   # Assemble test set.
   #===================================================================
   
   df_te = test
   
   df_te$pmean_f = fit$pmean_oos
   df_te$pmean_v = colMeans(fit_v$postpred)
   
   # Lower and upper bounds.
   df_te$lb_f = apply(fit$mcmcdraws_oos, 2, function(x) quantile(x,.025))
   df_te$ub_f = apply(fit$mcmcdraws_oos, 2, function(x) quantile(x,.975))
   
   df_te$lb_v = apply(fit_v$postpred, 2,function(x) quantile(x,.025))
   df_te$ub_v = apply(fit_v$postpred, 2,function(x) quantile(x,.975))
   
   # Full MCMC.
   df_te = cbind.data.frame(df_te, 'funMCMC'=as.data.frame(t(fit$mcmcdraws_oos)))
   df_te = cbind.data.frame(df_te, 'vMCMCM'=t(fit_v$postpred))
   
   
   #===================================================================
   # Combine test and train data frames for output.
   #===================================================================
   
   # Combine test and train output file
   out = rbind(df_tr, df_te)
   
   # Remove tsBART package.
   detach("package:tsbart", unload=TRUE)
   
   #===================================================================
   # Write logl output files for the log-density info.
   #===================================================================
   
   fwrite(lldf, paste0(out.csv, gsub('.csv','',filename),'_logl.csv'), row.names=FALSE)
   
   return(NULL)
   
}

#####################################################################
# Fit models for all data sets in directory above.
#####################################################################

# source is location of csv files.  
# out.csv is location to output result files, which are used to assemble table.
source = paste0(getwd(), "/data/supp-simulations-table-04/")
out.csv = paste0(getwd(),"/data/supp-simulations-table-04/results/")

# Iterate through files and generate comparisons.
counter = 0

for(f in files){
   
   run = funbartCompare(filepath=source, filename=f, probit=FALSE, out.csv)
   
   counter = counter + 1
   print(paste0('SIMULATION ',counter,' of ',length(files),' completed.'))
   
}

#####################################################################
# Assemble table.
#####################################################################

#===================================================================
# Read all file names.
#===================================================================

# Keep only files with logl in name.
files <- grep('logl', list.files(out.csv), value=TRUE)

# List of unique scenarios.
scenarios = unique(gsub( "-rep.*$", "", files))

#===================================================================
# Table assembly
#===================================================================

table = as.data.frame(matrix(NA,nrow=0,ncol=5))
colnames(table) = c('in-sample BART', 'in-sample tsBART','out-of-sample BART', 'out-of-sample tsBART','p-value')

# Loop through scenarios.
for(i in 1:length(scenarios)){
   
   # Grab scenario i.
   sc = paste0(scenarios[i],'-')
   
   # Read all files into a list.
   fs = files[which(substr(files, 1,nchar(sc))==sc)]
   
   all.files <- paste0(out.csv,fs)
   mylist <- lapply(all.files, fread)
   
   # Convert list to data frame and add scenario name.
   mydata <- rbindlist(mylist)
   mydata$scenario = rep(gsub('_logl.csv','',fs), each=4000)
   
   # Calculate means for each scenario.
   meanstbl = mydata %>% group_by(dataset, model) %>% summarize('ll' = mean(logdensity))
   
   table[i,] = c(meanstbl[2,3], meanstbl[1,3], meanstbl[4,3], meanstbl[3,3])
   rownames(table)[i] = substr(sc, 17, nchar(sc)-1)
   
   # #------------------------------------
   # # In-sample p-value.
   # #------------------------------------
   # temp_fun = mydata %>% 
   #    filter(dataset=="in-sample", model=="functional") %>% 
   #    dplyr::group_by(scenario) %>%
   #    summarize('logdens'=mean(logdensity))
   # 
   # temp_van = mydata %>% 
   #    filter(dataset=="in-sample", model=="vanilla") %>% 
   #    dplyr::group_by(scenario) %>%
   #    summarize('logdens'=mean(logdensity))
   # 
   # temp_p = wilcox.test(x=temp_fun$logdens, y=temp_van$logdens)$p.value
   # table[i,3] = temp_p
   
   #------------------------------------
   # Out-of-sample p-value.
   #------------------------------------
   temp_fun = mydata %>% 
      filter(dataset=="out-of-sample", model=="functional") %>% 
      dplyr::group_by(scenario) %>%
      summarize('logdens'=mean(logdensity))
      
   temp_van = mydata %>% 
      filter(dataset=="out-of-sample", model=="vanilla") %>% 
      dplyr::group_by(scenario) %>%
      summarize('logdens'=mean(logdensity))
   
   temp_p = wilcox.test(x=temp_fun$logdens, y=temp_van$logdens)$p.value
   table[i,5] = temp_p
   
}

# Rounding and p-value formatting.
table[,-5] = round(table[,-5],2)
table[,5] = round(table[,5],3)
table[,5] = as.character(table[,5])
table[,5][which(table[,5]=="0")] = "<0.001"

order = c(5,8,6,7,9,12,10,11,1,4,2,3)
table = table[order,]

# Add p and n info to table, instead of keeping it in the row names.
p = rep(c(4,8,20),each=4)
n = rep(c(100,500,1000,2500),times=3)
table = cbind.data.frame(p,n,table)



fwrite(table,"./output-files/table-02.csv")

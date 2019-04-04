# Function to run functional BART and vanilla BART on data set, and return
# model fits and log-densities for comparison.

# Reads simulated BART file and runs simulation.

# Returns list containing two data frames: model fits and log-likelihood info.


funbartCompare = function(filepath, filename, probit=FALSE, out.csv){
   
   require(tidyverse)
   require(data.table)
   require(funbart)
   require(fastbart)
   
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
   library(fastbart)
   
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
   fit_v = bartRcppClean(y_ = y, x_ = t(cbind.data.frame(xx,ti)), # obsvervations must be in *columns*
                         xpred_ = t(cbind.data.frame(x_pred,'ti' = ti_pred)), #Dpred,#[,-2,drop=F],
                         xinfo_list = cutpoints,
                         burn, nd, m,
                         lambda, nu, kfac=2,
                         paste0(getwd(),"/trees.txt"),
                         RJ=FALSE)
   
   print('Vanilla BART completed.')
   
   #############################################################################################
   ###   1. Functional BART.
   #############################################################################################
   
   # Remove funbart package.
   detach("package:fastbart", unload=TRUE)
   
   # Install and load funbart package.
   library(funbart)
   
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
   # Fit BART model.
   #===================================================================
   
   m = 200; burn = 200; nd = 1000
   
   if(probit){
      
      # Calculate offset
      phat = mean(yobs)
      offset = qnorm(phat)
      
      fit = bartRcppProbit(y_ = y,
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
      fit = bartRcppClean(y_ = y, x_ = t(xx), # obsvervations must be in *columns*
                          t_ = ti, # new input for time points; default 1.
                          xpred_ = t(x_pred), #Dpred,#[,-2,drop=F],
                          tpred_ = ti_pred,
                          xinfo_list = cutpoints,
                          burn = burn, nd = nd, m = m,
                          lambda=lambda, nu=nu, ecross = exp_cross,
                          treef_name_ =  paste0(getwd(),'_trees.txt'))
   }
   
   print('Functional BART completed.')
   
   #############################################################################################
   ###   3. Comparison between functional and vanilla.
   #############################################################################################
   
   #===================================================================
   # In-sample
   #===================================================================
   
   # Functional BART
   check_is = checkFit(y=y, 
                       mcmcdraws = fit$mcmcdraws,
                       sig = fit$sigma,
                       probit,
                       doWaic=TRUE,
                       yobs)
   
   # Vanilla BART
   check_is_v = checkFit(y=y, 
                         mcmcdraws = fit_v$postfit,
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
                          mcmcdraws = fit_v$postpred,
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
   
   # Remove funbart package.
   detach("package:funbart", unload=TRUE)
   
   #===================================================================
   # Write logl output files for the log-density info.
   #===================================================================
   
   fwrite(lldf, paste0(out.csv, gsub('.csv','',filename),'_logl.csv'), row.names=FALSE)
   
   return(NULL)
   
}


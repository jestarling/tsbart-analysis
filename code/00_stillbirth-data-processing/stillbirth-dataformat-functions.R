# OB FUNCTIONS

library(Rcpp)
library(RcppArmadillo)

#===================================================================
# Function for test/train split.
#===================================================================

testtrain = function(df, response='nd', test_pct = .3){
   
   ### Stratified sampling based on all combos of fd or nd and diab/htn and gest age.
   ### Requires an 'id' column identifying unique obs.
   
   require(dplyr)    # For left_join.
   
   # Set up grid of groups from which to sample.
   grps = expand.grid('status' = c(0,1),
                      'diab_htn' = factor(c('Neither','Diabetes','Htn','Both'), 
                                          levels = c('Neither', 'Diabetes', 'Htn','Both')),
                      'gest_age' = sort(unique(df$gest_age)))
   
   grps$strata = as.numeric(rownames(grps))
   colnames(grps)[1] = paste(response)
   
   # Pull in strata to data frame.
   df = left_join(df, grps, by=colnames(grps)[1:3])
   
   # Get rid of strata which do not exist in the data set.
   strata_nodata = grps$strata[-which(grps$strata %in% df$strata)]
   grps = grps %>% filter(!strata %in% strata_nodata)
   
   # Sample sizes of test set. 
   test_sampsize = round(xtabs(~df$strata) * test_pct,0)
   
   # loop through strata.
   test_id = NULL
   train_id = NULL
   
   for(i in sort(unique(df$strata))){
      
      temp = df %>% filter(strata==i)
      
      # If not enough rows available, allocate them to train.
      if(nrow(temp) * test_pct < 1){
         temp_train_id = temp$id
         temp_test_id = NULL
      } else{
         temp_test_id = sample(temp$id, size = test_sampsize[which(names(test_sampsize)==i)], replace=F)
         temp_train_id = temp$id[-which(temp$id %in% temp_test_id)]
      }
      
      test_id = c(test_id, temp_test_id)
      train_id = c(train_id, temp_train_id)
   }
   
   # Assemble test/train data sets and return.
   
   train = df %>% filter(id %in% train_id)
   test = df %>% filter(id %in% test_id)
   
   # Return output
   return(list(#'train' = train,
      #'test' = test,
      'train_id' = train_id,
      'test_id' = test_id))
}

#***********************************************************************
#***********************************************************************
#*** Functions for Data Processing, Test/Train Split      **************
#***********************************************************************
#***********************************************************************

#=======================================================================
#=== Utility functions for formatting data set.
#=======================================================================

#-----------------------------------------------------------------------
# Sets levels ordering for factor variables.
#-----------------------------------------------------------------------

fix_levels = function(data){
   
   if('momEthnicity' %in% colnames(data)){
      data$momEthnicity = factor(data$momEthnicity, levels=c('White NonHisp','Black NonHisp', 'Hispanic', 'Other'))
   }
   
   if('diab_htn' %in% colnames(data)){
      data$diab_htn = factor(data$diab_htn, levels = c('Neither','Diabetes','Htn','Both'))
   }
   
   if('agegrp' %in% colnames(data)){
      data$agegrp = factor(data$agegrp, levels = c('<20','20-24','25-29','30-34','35-39','40-44','45+'))
   }
   
   # Fix levels for birthwtQ and wtgainQ so that [.25-.75) is baseline level.
   if('birthwtQ' %in% colnames(data)){
      base_level = '[0.25,0.75)'
      base_idx = which(levels(data$birthwtQ)==base_level)
      other_idx = (1:length(levels(data$birthwtQ)))[-base_idx]
      data$birthwtQ = factor(data$birthwtQ, levels = c(levels(data$birthwtQ)[c(base_idx, other_idx)]))
   }
   
   if('wtgainQ' %in% colnames(data)){
      base_level = '[0.25,0.75)'
      base_idx = which(levels(data$wtgainQ)==base_level)
      other_idx = (1:length(levels(data$wtgainQ)))[-base_idx]
      data$wtgainQ = factor(data$wtgainQ, levels = c(levels(data$wtgainQ)[c(base_idx, other_idx)]))
   }
   
   return(data)
}

#-----------------------------------------------------------------------
# Drops rows from the data frame df, where variables specified in vars are 99 or 'Unknown' (case sensitive).
#-----------------------------------------------------------------------

drop_unknown = function(df,vars){
   
   # For dropping all unused levels in df.
   require(gdata)
   
   # Drop rows and unused levels.
   subset = df[!rowSums(df[,vars,drop=FALSE]==99 | df[,vars,drop=FALSE]=='Unknown' ) > 0, ]
   subset = drop.levels(subset)
   return(subset)
}

#-----------------------------------------------------------------------
# Drops rows from the data frame df, where variables specified in vars are 99 or 'Unknown' (case sensitive).
#-----------------------------------------------------------------------

drop_NA = function(df,vars){
   
   # For dropping all unused levels in df.
   require(gdata)
   
   # Drop rows and unused levels.
   subset = df[!rowSums(is.na(df[,vars,drop=FALSE])) > 0, ]
   subset = drop.levels(subset)
   return(subset)
}

#-----------------------------------------------------------------------
# Utilities to repeat rows in a data frame; ie for setting up xnew_full grid for 
# new patient panel.
#-----------------------------------------------------------------------

rep.row <- function(r, n){
   # r = data frame.
   # n = number of reps for each row in data frame.
   require(plyr)
   colwise(function(x) rep(x, n))(r)
}

ga.grid = function(df, incr){
   # df = data frame of new patients
   # incr = increment for weeks (ga) to create grid.  (.25 = divide grid into quarter-weeks)
   
   # Returns df_full; df, expanded into a grid where obs are repeated at each ga incr from 0 to 8 weeks.
   # Add weeks variable.
   wks = seq(0,8,by=incr)
   n_wks = length(unique(wks))
   
   # If id is not a column in df, add it.
   if(! 'id' %in% colnames(df)){
      df$id = 1:nrow(df)
   }
   
   # Add repeated rows for each week in grid.
   df_full = rep.row(df,n_wks)
   df_full$weeks = rep(wks, each=nrow(df))
   
   return(df_full)
}

#=======================================================================
#===   Split test/train data based on stratified sample   ==============
#=======================================================================

testtrain = function(data_nd, data_fd, test_pct = .3){
   
   ### Stratified sampling based on all combos of fd/nd and diab_htn.
   
   
   # Library for left_join.
   require(dplyr)
   
   # Combine nd and fd into one temp variable for convenience.
   data_nd$status = as.factor(ifelse(data_nd$nd==0 & data_nd$fd==0, 'surv',
                                     ifelse(data_nd$nd==1,'nd','fd')))
   
   grps = expand.grid('status' = as.factor(c('surv','nd','fd')), 
                      'diab_htn' = factor(c('Neither','Diabetes','Htn','Both'), 
                                          levels = c('Neither', 'Diabetes', 'Htn','Both')))
   grps$strata = rownames(grps)
   
   # Match two columns to get strata group (out of possible 12 groups).
   data_nd$strata = left_join(data_nd, grps, by = c('status','diab_htn'))$strata
   
   # Sample size of test set.
   test_sampsize = round(xtabs(~data_nd$strata) * test_pct,0)
   
   # Sample ids which are part of train data set.
   data_nd = data_nd[order(data_nd$strata),]
   id_test = strata(data_nd, stratanames='strata', size=test_sampsize, method = 'srswor')$ID_unit 
   
   # Split test and train data sets.
   test_nd = data_nd   %>% filter(id %in% id_test)
   train_nd = data_nd  %>% filter(!(id %in% id_test))
   test_fd = data_fd   %>% filter(id %in% id_test)
   train_fd = data_fd  %>% filter(!(id%in% id_test))
   
   # Return output
   return(list('test_nd' = test_nd,
               'train_nd' = train_nd,
               'test_fd' = test_fd,
               'train_fd' = train_fd))
}


#=======================================================================
#=== Utility function to create quantiles conditional on each gestational age.
#=======================================================================

conditQuantile = function(df, var, prbs = c(0,.25,.75,1), condit = "ga"){
   # Creates conditional quantiles for variable var, using probs, conditional on condit.
   # Writes csv file key for levels/percentiles for reference.
   # Data is from data frame df.
   
   ### Example of cuts function to show how binning works.
   ### test = rep(1:10, each=2)
   ### cut(test, breaks=c(1,3,7,10), include.lowest=T, right=F)
   
   # Initialize key matrix and vector to store results.
   key = matrix(NA, nrow = length(unique(data[,paste(condit)])), ncol = length(prbs)-1, dimnames = list(sort(unique(data[,paste(condit)])),prbs[-length(prbs)]))
   qtiles = rep(NA, nrow(data))
   vals = rep(NA, nrow(data))
   
   # Colnames for key.
   for(i in 1:(length(prbs)-1)){
      if(i == length(prbs)-1){
         colnames(key)[i] = paste0('[', prbs[i], ',', prbs[i+1], ']')
      } else{
         colnames(key)[i] = paste0('[', prbs[i], ',', prbs[i+1], ')')
      }
   }
   
   # Loop through gestational ages.
   for(ga in sort(unique(data[,paste(condit)]))){
      
      # Create quantile breakpoints.
      ga_idx = which(data$ga == ga)
      qs = quantile(data[ga_idx,paste(var)], probs = prbs)
      
      # Assign poins to quantile.  Breaks give ranges, including upper and lower limits.  0-.2, .2-.4, ..., .8-1.
      temp_vals = cut(data[ga_idx,paste(var)], breaks = qs, include.lowest = T, right = F)
      temp = cut(data[ga_idx,paste(var)], breaks = qs, include.lowest = T, right = F, labels = colnames(key))
      
      qtiles[ga_idx] = as.character(temp)
      vals[ga_idx] = as.character(temp_vals)
      
      # Assign breakpoints to key.
      key[ga - min(sort(unique(data[,paste(condit)]))  )+1,] = levels(temp_vals)
   }
   
   # Output csv key.
   write.csv(key, paste0('Data/DataForModels/quantile_key_',var,'_condit_on_',condit,'.csv'))
   
   # Return results.
   return(list('key' = key, 'vals' = vals, 'qtiles' = qtiles))
   
}


#***********************************************************************
#***********************************************************************
#*** Functions for Plotting Output of Model (Coefs, Patient Panel)   ***
#***********************************************************************
#***********************************************************************


#===================================================================================
#===   Function for plotting and outputting coefficient confidence interval info  ==
#===================================================================================

myCI.plot = function(model, out.csv, out.fig, fname, exponentiate = FALSE, wd = 16, ht = 12, xlab = 'Regression Coefficients'){
   #-----------------------------------------------------------
   # FUNCTION:    For plotting confidence intervals for fixed effects, and 
   #              outputting results to csv file. 
   #-----------------------------------------------------------
   # INPUTS:      model = a model object.
   #              out.csv = Path used for saving csv files.
   #              out.fig = Path used for saving figures.
   #              fname = a file name, used for both plots. Ex: "myfile"
   #              exponentiate = FALSE for logit plots, TRUE for odds ratio plots.
   #               wd = width of output pdf
   #               ht = height of output pdf
   #-----------------------------------------------------------
   # OUTPUTS:     Function does not return output.
   #              Saves pdf and csv files to specified folder.
   #-----------------------------------------------------------
   
   require(plotrix) 		   #For plotting CIs
   
   #-----------------------------------------------------------
   # Set up CI data frame for model.
   #-----------------------------------------------------------
   
   # Output coefficients which do not vary by state, including baseline intercept.
   coefs.fixed = summary(model)$coefficients
   
   sig = ifelse(coefs.fixed[,4] < 0.001, '***',
                ifelse(coefs.fixed[,4] < 0.01, '**',
                       ifelse(coefs.fixed[,4] < 0.05, '*',
                              ifelse(coefs.fixed[,4] < 0.1, '+', ''))))
   coefs.fixed = cbind.data.frame(coefs.fixed, sig)
   
   # Set up CI bounds for plotting.
   CI = data.frame('lb' = coefs.fixed[,1] - 1.96 * coefs.fixed[,2],
                   'ub' = coefs.fixed[,1] + 1.96 * coefs.fixed[,2],
                   'hat' = coefs.fixed[,1],
                   'sig' = coefs.fixed[,5],
                   row.names = rownames(coefs.fixed))
   
   # If exponentiate=TRUE, exponentiate the confidence intervals.
   if(exponentiate==TRUE) {
      CI[,1:3] = exp(CI[,1:3])
   }
   
   #-----------------------------------------------------------
   # Confidence Intervals
   #-----------------------------------------------------------
   
   # Plot confidence intervals.
   pdf(paste(out.fig, fname,".pdf", sep=''), width=wd, height=ht, family='ArialMT')
   par(oma = c(5.5,4,0,4), mar = c(16, 4, 4, 2) + 0.1, ps=12)
   
   # Set up dynamic title.
   main.title = ''
   
   # Main plot.
   plotCI(x = CI[-1,3], li = CI[-1,1], ui = CI[-1,2],
          xaxt='n',
          xlab = '',
          ylab = '95% Confidence Interval',
          main = main.title,
          cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
   
   # Legend
   par(xpd=T)
   legend('bottomright',c('+     p-value < 0.1', 
                          '*      p-value < 0.05',
                          '**    p-value < 0.01',
                          '***   p-value < 0.001'), cex = 1, inset=c(0,-.65))
   par(xpd=F)
   
   # Red line at zero; or at one, if exponentiated (multiplicative effect).
   if(exponentiate==FALSE){
      abline(h=0,lty=2,col='red',lwd=2) 
   } else{
      abline(h=1,lty=2,col='red',lwd=2) 
   }
   
   # Add axis labels.
   axis(1, at = 1:(nrow(CI)-1), labels=FALSE, CI$sig[-1], cex = 1.5)
   
   # X-axis label names.
   labnames = trimws(paste(rownames(CI)[-1], CI$sig[-1]), 'right')
   
   # Add x-axis labels for individual plots.
   text(seq(1,nrow(CI)-1,by=1), 
        y = par("usr")[3] - .1,
        labels = labnames, 
        srt = 80, pos = 2, xpd = TRUE, cex = 1)
   
   # Add x-axis for page.
   mtext(paste(xlab), side=1, line=3, outer=TRUE, cex = 1.2)
   dev.off() #End plot.
   
   #-----------------------------------------------------------
   # CSV Output
   #-----------------------------------------------------------
   write.csv(CI,paste(out.csv, fname, '.csv',sep=''),na='') 
   
} #End myCI.plot function.

#=======================================================================
#===   Optimal gestational age generator   =============================
#=======================================================================

optimalga = function(ptinfo, nd_model, fd_model, incr=.25, pred_only=FALSE){
   #-----------------------------------------------------------
   # FUNCTION:    Returns two data frames:
   #              1. ptinfo: Contains optimal ga, with bounds, for each patient.
   #                         If ptinfo contains a gestational age for each new pt, also contains
   #                         predicted values for phat_fd and phat_nd at that gestational age.
   #              2. curve_info: Contains an array risk curve information for all patients.
   #                          Dimensions are: (ga, risk curve columns, patient)
   #                          Columns for each patient's matrix are:
   #                          nd_curve, nd_se, nd_lb, nd_ub, fd_curve, fd_se, nd_lb, nd_ub.
   #                          Note that these are not true CIs; more like Bayesian credible intervals.
   #-----------------------------------------------------------
   # INPUTS:      ptinfo = Data frame of patient characteristics.
   #              nd_model = Model object for neonatal death.
   #              fd_model = Model object for fetal death.
   #              incr = increments for weeks (ga)
   #              pred_only = if pred_only = false, only output predict() objects for new data set; do not return grid of curves.
   #-----------------------------------------------------------
   # OUTPUTS:     optga = Optimal gestational age for each patient.
   #              optga_lb = Lower bound for plausible optimal ga for each patient.
   #              optga_ub = Upper bound for plausible optimal ga for each patient.
   #              plots = List of plot objects for each patient.
   #-----------------------------------------------------------
   # DEPENDENCIES: optimalga_cpp_util()
   #-----------------------------------------------------------
   
   # Expand grid of patient info to include all week increments for prediction.
   baseline_ga = 34
   wks = seq(0,8,by=incr)
   ptinfo_full = ga.grid(ptinfo, incr)
   pt_ids = ptinfo_full$id
   
   # Calculate predicted values and standard errors for new patients for only observed values.
   # Do this only if a gestational age is provided in pinfo data frame.
   phat_nd = rep(NA,nrow(ptinfo))
   phat_fd = rep(NA,nrow(ptinfo))
   
   if("weeks" %in% colnames(ptinfo)){
      phat_nd = predict(nd_model, newdata=ptinfo, type='response', se.fit=F)
      phat_fd = predict(fd_model, newdata=ptinfo, type='response', se.fit=F)
   }
   
   # Calculate predicted values and standard errors for new patients for full grid, for curve info.
   newpt_nd = predict(nd_model, newdata=ptinfo_full, type='response', se.fit=T)
   newpt_fd = predict(fd_model, newdata=ptinfo_full, type='response', se.fit=T)
   
   if(pred_only==T){
      return(list('pred_nd' = phat_nd, 'pred_fd' = phat_fd))
   }
   
   # Estimates.
   line_fd = newpt_fd$fit
   line_nd = newpt_nd$fit
   
   # Standard errors.
   se_fd = newpt_fd$se.fit
   se_nd = newpt_nd$se.fit
   
   # Lower and upper bounds for CIs.
   lb_fd = line_fd - 1.96*se_fd
   ub_fd = line_fd + 1.96*se_fd
   lb_nd = line_nd - 1.96*se_nd
   ub_nd = line_nd + 1.96*se_nd
   
   # Number of new patients.
   n_new_pts = nrow(ptinfo)
   
   out = optimalga_cpp_util(line_fd, line_nd, se_fd, se_nd, n_new_pts, wks, baseline_ga, pt_ids)
   
   dimnames(out$curve_info) = list(wks + baseline_ga, 
                                   c('nd_curve','nd_se','nd_lb','nd_ub','fd_curve','fd_se','fd_lb','fd_ub'),
                                   paste0('Patient',1:n_new_pts))
   
   #---------------------------------------------
   # Return output
   #---------------------------------------------
   
   # Info for each patient in the input. (Curve_info, in contrast, is over the entire GA grid.)
   ptinfo = cbind.data.frame('opt_ga' = out$opt_ga, 
                             'opt_ga_lb' = out$opt_ga_lb,
                             'opt_ga_ub' = out$opt_ga_ub,
                             'phat_nd' = phat_nd,
                             'phat_fd' = phat_fd)
   
   return(list('ptinfo' = ptinfo, 'curve_info' = out$curve_info))
   
}

#=======================================================================
#==   Plot predicted curves for new patients (patient panel).
#=======================================================================

print_panel = function(optimalga_output, fname, nr=4, nc=4, ht=11, wd=8.5){
   #-----------------------------------------------------------
   # FUNCTION:    For plotting predicted values and CIs for
   #              panel of new (predicted) patients.
   #-----------------------------------------------------------
   # INPUTS:      optimalga_output = the output object of the optimalga function.
   #              fname = filename
   #              nr, nc = number of rows, columns for panel on each pg
   #-----------------------------------------------------------
   # OUTPUTS:     Saves pdf file to specified folder.
   #-----------------------------------------------------------   
   
   # Plot predicted outcomes for grid of new patients.
   pdf(paste('R/Figs/',fname,'.pdf',sep=''),height=ht,width=wd)
   
   par(mfrow=c(nr,nc), oma = c(4.5,3,3,4) + .05)
   ptlabel_cex = 1
   
   # y-axis details.
   ylims = c(0,.006)
   yaxis_at = seq(ylims[1],ylims[2],length=7)
   
   # Number of new patients, and weeks grid.
   n_new_pts = length(dimnames(optimalga_output$curve_info)[[3]])
   wks = as.numeric(dimnames(optimalga_output$curve_info)[[1]])
   
   #-----------------------------------------------------------
   # Loop through patients.
   #-----------------------------------------------------------
   
   # Loop through patients.
   for(i in 1:n_new_pts){
      
      # Extract values from input object.
      line_nd_idx = optimalga_output$curve_info[,'nd_curve',i]
      line_fd_idx = optimalga_output$curve_info[,'fd_curve',i]
      
      lb_nd_idx = optimalga_output$curve_info[,'nd_lb',i]
      lb_fd_idx = optimalga_output$curve_info[,'fd_lb',i]
      
      ub_nd_idx = optimalga_output$curve_info[,'nd_ub',i]
      ub_fd_idx = optimalga_output$curve_info[,'fd_ub',i]
      
      opt_ga    = optimalga_output$ptinfo$opt_ga[i]
      opt_ga_lb = optimalga_output$ptinfo$opt_ga_lb[i]
      opt_ga_ub = optimalga_output$ptinfo$opt_ga_ub[i]
      
      #---------------------------------------------
      # Plotting
      #---------------------------------------------
      
      # Plot lines for estimated coefficients. Calls new plot.
      plot(wks, line_fd_idx, type='l', ylim=ylims, col='grey10', lty=5,
           main='', yaxt = 'n', xaxs="i", yaxs="i",xlab='', ylab = '',cex.main=.75)
      
      # Main title for each plot.
      mtext(paste0('Patient ', i), cex=ptlabel_cex, line=2)
      mtext(paste0(opt_ga,' wks'),cex=ptlabel_cex, line=1)
      
      
      # Shaded rectangle for plausible gestational age bounds.
      rect(opt_ga_lb, ylims[1], opt_ga_ub, ylims[2], col =  adjustcolor("grey70", alpha.f = 0.3), border=NA)
      
      # Add dotted line at zero and y-axis labels.
      abline(h=0,lty=3,col='blue')
      axis(side=2, at=yaxis_at, labels=round(yaxis_at,3))
      
      # Vertical line at optimal gestational age.
      abline(v=opt_ga, col='firebrick3',lty=2,lwd=1)
      
      # Plot CIs.
      polygon(x = c(wks, rev(wks)), y = c(ub_fd_idx, rev(lb_fd_idx)), col =  adjustcolor("grey40", alpha.f = 0.30), border = NA) # deepskyblue4
      polygon(x = c(wks, rev(wks)), y = c(ub_nd_idx, rev(lb_nd_idx)), col =  adjustcolor("grey40", alpha.f = 0.30), border = NA) # darkorange3
      
      # Plot lines for estimated coefficients. 
      lines(wks, line_fd_idx, type='l', col='grey10', lwd=2, lty=5)
      lines(wks, line_nd_idx, type='l', col='grey30', lwd=2)
      
      
      # Text to add gestational age range.
      # mtext(paste0('(',opt_ga_lb,' wks - ',opt_ga_ub,' wks)'), cex=ptlabel_cex, line=0)
      mtext(paste0('(',opt_ga_lb,' - ',opt_ga_ub,')'), cex=ptlabel_cex, line=0)
      
      # Overall title and legend at bottom of each page.
      if(i%%(nr*nc)==0) {
         
         mtext("Optimal Delivery to Minimize Stillbirth & Neonatal Death: Patient Panel", outer = TRUE, cex = 1.5)
         
         mtext(text="Gestational age (weeks)",side=1,line=0,outer=TRUE)
         mtext(text="Stillbirth and Neonatal Death Risks",side=2,line=0,outer=TRUE)
         
         legend('bottom', inset=-.75 ,legend=c('fetal death','neonatal death','optimal delivery'),
                col=c('grey10','grey30','firebrick3'),
                lwd=c(2,2,2),lty=c(5,1,3),bty="n",xpd=NA, cex=1.35)
      }
   }
   
   dev.off()
   
}

#****************************************************************************
#****************************************************************************
#*** Helper functions for post-hoc analysis, Mandujano comparison, etc   ****
#****************************************************************************
#****************************************************************************

#===========================================================================================
#=== Compute the binomial log-likelihood function for a specific parameter value.   ========
#===========================================================================================

binom.loglik = function(n, p, x){
   # n = vector of sample sizes (n_1, ..., n_n).
   # y = vector of counts (y_1, ..., y_n).
   # p = vector of probs (p_1, ..., p_n).
   out = dbinom(x=x, size=n, prob=p, log=T)
   return(sum(out))
}

#===========================================================================================
#===   Helper function to generate Mandujano table 2 for a given collapsed data set   ======
#===========================================================================================

mand_tblgen_collapsed = function(data){
   
   require(tidyverse)
   
   #-----------------------------------------------------------------------
   # Table 2 Count Section
   #-----------------------------------------------------------------------
   
   tbl = as.data.frame(matrix(NA,nrow=9,ncol=8))
   
   tbl = data  %>% 
      group_by('ga' = as.factor(ga)) %>%
      summarise('n_ga' = sum(nd_total), 'lb' = sum(lb), 'fd' = sum(fd_yes), 'nd' = sum(nd_yes))
   
   # Adjust tf so that it is a running count of all fetuses.
   tbl$tf = rep(0,9)
   
   tbl$tf[1] = sum(tbl$n_ga)
   for(i in 2:nrow(tbl)){
      tbl$tf[i] = tbl$tf[i-1] - tbl$lb[i-1] - tbl$fd[i-1]
   }
   
   tbl = tbl[, -2]
   
   #-----------------------------------------------------------------------
   # Table 2 Analysis Section
   #-----------------------------------------------------------------------
   
   # Populate cum % fd.
   tbl$cumfd = NA
   
   for (i in 1:nrow(tbl)){
      tbl$cumfd[i] = round(sum(tbl$fd[1:i]) / sum(tbl$fd) * 100, 2)
   }
   
   # Populate ndr.
   tbl$ndr = round(tbl$nd / (tbl$lb / 1000), 2)
   
   # Populate fdrre.
   deliveries = tbl$lb + tbl$fd
   tbl$fdrre = NA
   
   for (i in 1:nrow(tbl)){
      
      # Numerator: sum of fetal death for ga >= current ga i.
      num = sum(tbl$fd[i:nrow(tbl)])
      
      # Denominator.  If i=last row of table, deliveries only.  Else deliveries at current ga + remaining total fetuses for ga i+1.
      denom = ifelse(i==nrow(tbl), deliveries[i], deliveries[i] + tbl$tf[i+1]  )
      tbl$fdrre[i] = round(num/denom * 1000, 2)
   }
   
   #Divide by 1000 for ndr and fdrre.
   tbl$fdrre = tbl$fdrre / 1000
   tbl$ndr = tbl$ndr / 1000
   
   # Return output table.
   return(tbl)
}


#===========================================================================================
#===   Helper functions to expand the predictions on the collapsed data set for       ======
#===   creating histograms (for assessment of range of optga values)                  ======
#===========================================================================================

#-----------------------------------------------------------------------
# Utility function to get expanded vector of characteristics for each patient, 
# for plotting histograms conditional on that characteristic.
#-----------------------------------------------------------------------

expandchar = function(df=NULL, var, counts){
   if(!is.null(df)){
      temp = rep(df[,paste(var)], counts)
   } else{
      temp = rep(var, counts)
   }
   return(temp)
}

#-----------------------------------------------------------------------
# Utility function for plotting histogram
#-----------------------------------------------------------------------

pt_hist = function(opt_ga_exp, char_exp=NULL, level=NULL, title=NULL, brks=20){
   
   # Creates histogram for optimal gestational ages for level of char with specified title.
   idx = NULL
   
   if(is.null(char_exp)){
      idx = 1:length(opt_ga_exp)
   } else{
      idx = which(char_exp==level)
   }
   
   hist(opt_ga_exp[idx], breaks=brks, col='grey', main=title, xlab='Optimal Gestational Age (wks)', xlim=c(34,42), xaxt='n')
   axis(1, at=34:42, labels=34:42)
   
   mtext(paste0('Mean ', round(mean(opt_ga_exp[idx]),2), ', SD ', round(sd(opt_ga_exp[idx]),2), ', Range (', min(opt_ga_exp[idx]), ', ', max(opt_ga_exp[idx]), ')'  ))
   
   abline(h=0)
}



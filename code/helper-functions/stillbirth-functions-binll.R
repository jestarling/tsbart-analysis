########################################################################
# Binomial log-likelihood.  For all methods, assesses quality of fit.
########################################################################

# Calculates individual log-likelihood contributions.
# Used by fit utility functions.
binloglik = function(yobs, probs){
   
   # Catch probs too close to 0.
   probs[which(probs<0.000001)] = 0.000001
   
   # Calculate contribution for each obs.
   #ll = yobs * probs + (1-yobs) * (1-probs)
   ll = dbinom(yobs, size=1, probs, log=TRUE)
   #binll = mean(ll)
   
   # return ll for each indivdual obs, also, in case of subgroup analysis needs.
   return(ll)
   #return(list('binll' = binll, 'binll_vec' = ll))
}

########################################################################
# Binomial log-likelihood on average, and by week.
########################################################################

binll_per_person = function(df, method='tsb', thresh=.5){
   
   # thresh is phat threshold for calling something a 1 vs 0.  0.5 is reasonable to start.
   require(tidyverse)
   
   # For weekly binll, take week l's average of the per-person binll for all people with pregnancies delivering at week l, or continuing past week l.
   n_wks = length(unique(df$gest_age34))
   wks = unique(df$gest_age34)
   
   col_idx = grepl(paste(method), colnames(df)); col_idx[which(colnames(df) %in% c('id','gest_age34'))] = TRUE
   df = df[,col_idx]
   
   df$phat_il = df[,3]
   df$yhat_il = ifelse(df$phat_il > thresh, 1, 0)
   df$yobs_il = test$fd
   df$binll = df[,6]
   df$fp = ifelse(df$yhat_il==1 & df$yobs_il==0, 1, 0)
   df$fn = ifelse(df$yhat_il==0 & df$yobs_il==1, 1, 0)
   df$misclass = ifelse(df$fp==1 || df$fn==1, 1, 0)
   
   # Df to hold logl and misclass error.
   logl = as.data.frame(matrix(NA, nrow=n_wks+1, ncol=6))
   colnames(logl) = c('method','gest_age34','logl','merr','fp','fn')
   logl$method = method
   
   for(l in 1:n_wks){
      temp_df = as.data.frame(df %>% filter(gest_age34 == wks[l]))
      logl$logl[l] = mean(temp_df$binll)
      logl$merr[l] = mean(temp_df$yhat_il!=temp_df$yobs_il)
      logl$fp[l] = mean(temp_df$yhat_il==1 & temp_df$yobs_il==0)
      logl$fn[l] = mean(temp_df$yhat_il==0 & temp_df$yobs_il==1)
      logl$gest_age34[l] = wks[l]
   }
   
   # For overall binll, take the average of the per-person binll.  (Since each person has difference # of time points.)
   logl$gest_age34[n_wks+1] = 'all'
   
   # Overall binll.
   binll_person = as.data.frame(df %>% dplyr::group_by(id) %>% summarise('binll' = sum(binll)))
   logl$logl[n_wks+1] = mean(binll_person$binll)
   
   # Add gest age.
   logl$gest_age = NA
   logl$gest_age[1:n_wks] = as.numeric(logl$gest_age34[1:n_wks]) + 33
   
   return(logl)
}

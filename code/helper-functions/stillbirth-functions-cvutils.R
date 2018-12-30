
########################################################################
# Extract oos fit data for patient panel.
########################################################################

extractFits = function(test, cv_tsb, cv_vb, cv_sp1, cv_sp2, cv_pensp, cv_tsb_default = NULL, adjust=F){
   
   #=======================================================================
   # Set parameters for adjusting the predicted oos probabilities.
   #=======================================================================
   
   padjust = 1
   rescale = 1
   
   if(adjust==T){
      padjust=50
      rescale=1000
   }
   
   #=======================================================================
   # Add probabilities to test dataframe.  Adjust if indicated.
   #=======================================================================
   
   # tsBART.
   test$phat_oos_tsb = cv_tsb$phat_oos / padjust * rescale
   test$phat_oos_tsb_lb = cv_tsb$phat_oos_lb / padjust * rescale
   test$phat_oos_tsb_ub = cv_tsb$phat_oos_ub / padjust * rescale
   test$phat_oos_tsb_binll = cv_tsb$binll_oos
   
   # tsBART with default settings.
   if(!is.null(cv_tsb_default)){
      test$phat_oos_tsb_default = cv_tsb_default$phat_oos / padjust * rescale
      test$phat_oos_tsb_default_lb = cv_tsb_default$phat_oos_lb / padjust * rescale
      test$phat_oos_tsb_default_ub = cv_tsb_default$phat_oos_ub / padjust * rescale
      test$phat_oos_tsb_default_binll = cv_tsb_default$binll_oos 
   }
   
   # BART.
   test$phat_oos_vb = cv_vb$phat_oos / padjust * rescale
   test$phat_oos_vb_lb = cv_vb$phat_oos_lb / padjust * rescale
   test$phat_oos_vb_ub = cv_vb$phat_oos_ub / padjust * rescale
   test$phat_oos_vb_binll = cv_vb$binll_oos
   
   # sp1.
   test$phat_oos_sp1 = cv_sp1$phat_oos / padjust * rescale
   test$phat_oos_sp1_lb = cv_sp1$phat_oos_lb / padjust * rescale
   test$phat_oos_sp1_ub = cv_sp1$phat_oos_ub / padjust * rescale
   test$phat_oos_sp1_binll = cv_sp1$binll_oos
   
   # sp2.
   test$phat_oos_sp2 = cv_sp2$phat_oos / padjust * rescale
   test$phat_oos_sp2_lb = cv_sp2$phat_oos_lb / padjust * rescale
   test$phat_oos_sp2_ub = cv_sp2$phat_oos_ub / padjust * rescale
   test$phat_oos_sp2_binll = cv_sp2$binll_oos
   
   # pensp.
   test$phat_oos_pensp = cv_pensp$phat_oos / padjust * rescale
   test$phat_oos_pensp_lb = cv_pensp$phat_oos_lb / padjust * rescale
   test$phat_oos_pensp_ub = cv_pensp$phat_oos_ub / padjust * rescale
   test$phat_oos_pensp_binll = cv_pensp$binll_oos
   
   
   #=======================================================================
   # Add labels and return panel.
   #=======================================================================
   
   # Extract right two chars of ids for labels.
   RIGHT = function(x,n){
      substring(x,nchar(x)-n+1)
   }
   
   test$label = ifelse(test$testtrain=="panel",
                       paste0('Patient ', RIGHT(test$id,2)),
                       paste0('Patient', test$id))
   
   return(test)
}

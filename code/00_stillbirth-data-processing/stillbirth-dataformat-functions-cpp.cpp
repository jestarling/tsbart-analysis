#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;
using namespace Rcpp;


//==========================================================================
//===   concat_data   ======================================================
//==========================================================================
// [[Rcpp::export]]
mat concat_data(mat df, int time_colnum, int response_colnum){
   // Returns a concatenated version of the data set, with rows repeated up as in 
   // Sparapani, NP Surv Analysis using BART; each row (obs) is repeated
   // at each timevar until (binary) response event, with fd=0 at times prior to fd event.)
      
   // Example: If time = gestational age, response=fd, then if fd occurs at ga=36, then
   // row appears with ga=34 and fd=0, ga=35 and fd=0, and ga=36 and fd=1. 
   
   //--------------------------------------------------------------------
   // Initialization
   //--------------------------------------------------------------------
   
   // Add ID column to matrix.
   int n = df.n_rows;               // Number of rows in data set.
   vec id = linspace(1,n,n);
   df.insert_cols(df.n_cols, id);
   
   // Initialize data structures.
   mat df_i = df.row(0);               // A 1-row data frame to hold an individual row of df.
   int time_i = 0;                    // Holds the time variable at each i.
   
   // Initialize key variables.
   int t_idx = time_colnum - 1;        // Col number (zero indexed) for time variable.
   int y_idx = response_colnum - 1;    // Col number (zero indexed) for response variable. 
   int baseline = min(df.col(t_idx));  // The min timepoint present in the data set.
   
   // Calculate number of rows needed for new matrix then initialize to zeros.
   int size_new = sum(df.col(t_idx) - baseline + 1) - n;
   mat df_new = zeros(size_new, df.n_cols);
   int k = 0; //Pointer for filling df_new.

   //--------------------------------------------------------------------
   // Loop through rows of data frame.
   //--------------------------------------------------------------------
   
   for(int i=0; i<n; i++){
       
       // Progress report.
       Rcout<< "Rows remaining to process: " << n-i<<std::endl;
       
       // Extract ith row and value of the time variable.
       df_i = df.row(i);
       time_i = df_i(t_idx);
       
       // Check if time_i > baseline.  If yes, set up additional rows.
       if(time_i > baseline){
         
          // Loop through new rows to add.
          for(int j = 0; j < time_i - baseline; j++){
             df_new.row(k) = df_i;     // Add new row.
             df_new(k,y_idx) = 0;      // Set response to 0.
             df_new(k,t_idx) = baseline + j; // Set time variable.
             k = k+1;
          }
       } // End if statement.
   } // End data frame loop.
   
   //--------------------------------------------------------------------
   // Concatenate original df with df_new.
   //--------------------------------------------------------------------
   
   Rcout << "Beginning matrix join" << std::endl;
   mat concat = join_cols(df,df_new);
   
   return(concat);
}

//=======================================================================
//===   Optimal gestational age generator   =============================
//=======================================================================
// [[Rcpp::export]]
List optimalga_cpp_util(vec line_fd, vec line_nd, 
                        vec se_fd, vec se_nd, 
                        double n_new_pts, vec wks, int baseline_ga, vec pt_ids){
   // Utility function to calculate optimal gestational age for each patient; speeds up optimalga function.

   // Calculate lower and upper bounds.
   vec lb_fd; vec ub_fd;
   vec lb_nd; vec ub_nd;
   
   lb_fd = line_fd - 1.96*se_fd;
   ub_fd = line_fd + 1.96*se_fd;     
   lb_nd = line_nd - 1.96*se_nd;
   ub_nd = line_nd + 1.96*se_nd;


   // Initialize data structures to hold optimal ga, optimal ga bounds, and curve info.
   vec opt_ga = zeros(n_new_pts);
   vec opt_ga_lb = zeros(n_new_pts);
   vec opt_ga_ub = zeros(n_new_pts);

    // Initialize cube (array) to hold curve info for each patient.
    cube curve_info = cube(wks.size(), 8, n_new_pts);


    // Initialize temp vector to hold incides for each patient.
    uvec idx;

    // Initialize temp vectors to hold line info.
    vec line_nd_idx;
    vec line_fd_idx;

    vec lb_nd_idx;
    vec lb_fd_idx;

    vec ub_nd_idx;
    vec ub_fd_idx;

    uvec nd_leq_fd;
    uvec fd_leq_nd;

   //-----------------------------------------------------------
   // Loop through patients.
   //-----------------------------------------------------------

   // Loop through patients.
   for(int i=0; i<n_new_pts; i++){

      // Housekeeping
      Rcpp::Rcout << "Processing patient " << i+1 << " of " << n_new_pts << endl;

      idx = find(pt_ids == i+1);

      //---------------------------------------------
      // Set up optimal gestational age estimate and boundaries.
      //---------------------------------------------

      // Calculate optimal gestational age:
      //   If nd curve > fd curve for all gas, use min nd.
      //   If fd curve > nd curve for all gas, use min fd.
      //   Else calculate ga at which lines cross (nd becomes <= fd).

      // Precache some values.
      line_nd_idx = line_nd(idx);
      line_fd_idx = line_fd(idx);
      
      lb_nd_idx = lb_nd(idx);
      lb_fd_idx = lb_fd(idx);

      ub_nd_idx = ub_nd(idx);
      ub_fd_idx = ub_fd(idx);

      nd_leq_fd.reset();
      fd_leq_nd.reset();

      nd_leq_fd = find(line_nd_idx <= line_fd_idx);
      fd_leq_nd = find(line_fd_idx <= line_nd_idx);

      //---------------------------------------------
      // Calculate optimal gestational age.
      //---------------------------------------------

       if(nd_leq_fd.is_empty() ){
               // nd always > fd; ga at which nd is minimized.
               opt_ga(i) = as_scalar(wks(find(line_nd_idx==min(line_nd_idx))));

      } else if(fd_leq_nd.is_empty() ){
               // fd always > nd; ga at which fd is minimized.
               opt_ga(i) = as_scalar(wks(find(line_fd_idx==min(line_fd_idx))));

      } else{
               // Gestational age at which lines cross (nd becomes <= fd)
               opt_ga(i) = wks(min(find(line_nd_idx <= line_fd_idx)));

      }
      
      //---------------------------------------------
      // Calculate optimal gestational age lower and upper ranges.
      //---------------------------------------------
         
      if(nd_leq_fd.is_empty() | fd_leq_nd.is_empty()){
         
         // If nd>fd for all, or fd>nd for all, we know exactly where minimizes; the opt_ga.
         opt_ga_lb(i) = opt_ga(i);
         opt_ga_ub(i) = opt_ga(i);
     

      } else{
         
         // GA lower limit:  max ga where lb_nd >= fd_ub. Bounded at 0 (34 weeks).
         uvec lb_test = find(lb_nd_idx >= ub_fd_idx);

         if(lb_test.is_empty()){
            opt_ga_lb(i) = 0;
         } else{
            opt_ga_lb(i) = wks(max(find(lb_nd_idx >= ub_fd_idx)));
         }

         // GA upper limit: max ga where ub_nd >= fd_lb.  Bounded at 8 (42 weeks).
         uvec ub_test = find(ub_nd_idx >= lb_fd_idx);
         if(ub_test.is_empty()){
            opt_ga_ub(i) = 8;
         } else{
            opt_ga_ub(i) = wks(max(find(ub_nd_idx >= lb_fd_idx)));
         }

      }

      //---------------------------------------------
      // Store curve info for patient i.
      //---------------------------------------------
     //  'nd_curve','nd_se','nd_CI_lb','nd_CI_ub','fd_curve','fd_se','fd_CI_lb','fd_CI_ub')

      curve_info(span(), span(0), span(i)) = line_nd_idx;
      curve_info(span(), span(1), span(i)) = se_nd(idx);
      curve_info(span(), span(2), span(i)) = lb_nd_idx;
      curve_info(span(), span(3), span(i)) = ub_nd_idx;
      
      curve_info(span(), span(4), span(i)) = line_fd_idx;
      curve_info(span(), span(5), span(i)) = se_fd(idx);
      curve_info(span(), span(6), span(i)) = lb_fd_idx;
      curve_info(span(), span(7), span(i)) = ub_fd_idx;
      

    } // End patient loop.

   //---------------------------------------------
   // Return output
   //---------------------------------------------

   return Rcpp::List::create(
      _["opt_ga"]    = opt_ga + baseline_ga,
      _["opt_ga_lb"] = opt_ga_lb + baseline_ga,
      _["opt_ga_ub"] = opt_ga_ub + baseline_ga,
      _["curve_info"]= curve_info
   ) ;

}

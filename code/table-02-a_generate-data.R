##########################################################
# Simulate datasets for tsBART simulation study.
##########################################################

rm(list=ls())
library(mosaic)
library(foreach)

##########################################################
# Function definitions.
##########################################################

sigmoid = function(x) 1/{1+exp(-x)}

# tuning the power in my "kick-up" function
rho = 1 + log(0.1)/log(0.75)

# baseline risk function
f1 = function(t) {
	0.075*(t)
}

# risk function with a big kick-up at t = 0.75
# risk at t = 1 is 5x over baseline
f2 = function(t) {
	0.75*pmax(0.75, t)^rho
}


# Display these two functions and the kick-up point.
curve(f2, col='red', ylim = c(0, 1))
curve(f1,add=TRUE, col='black')
abline(v=0.75)

##########################################################
# Next, we set up various interpolating functions.
# These will always return something between 0 and 1.
# The return value is intended to be the weight on f2, i.e. the function with a kick-up.
# In all cases, x is assumed to be 10-dimensional:
# First five dimensions matter, next five don't
#
# Three types of interpolating functions:
# 1. Linear (wx_linear).
# 2. Linear with interaction (wx_linear_interaction).
# 3. Nonlinear with interaction (wx_nonlinear_interaction).
##########################################################

#---------------------------------------------------------
# 1. Linear without interaction.
#---------------------------------------------------------

# Define linear function.
wx_linear = function(x) {
	psi =  5*(x[1] - x[2] + x[3] - x[4])
	sigmoid(psi)
}

# Simulate one random observation, with 10-dim x.
wx_linear(runif(10))

# Simulate random x's to calibrate level.
t_grid = seq(0,1,by=0.01)
test1 = do(50)*{
	x = runif(10)
	w = wx_linear(x)
	0.25*x[5] + w*f2(t_grid) + (1-w)*f1(t_grid) #.75 from board is absorbed in f1, f2.
}

# Plots hazard functions.
plot(t_grid, .test1[1,], ylim=c(0,1), type='l',
	col=rgb(0.1,0.1,0.1,0.1))
for(i in 1:50) {
	lines(t_grid, test1[i,],col=rgb(0.1,0.1,0.1,0.1))
}

#---------------------------------------------------------
# 1. Linear with interaction.
#---------------------------------------------------------

# There are four basic shapes:  start hi/low, kick up/or not.

## Linear single-index with an interaction
# With an x1-x2 and an x3-x4 interaction term.
wx_linear_interaction = function(x) {
	psi = 5*(x[1] - x[2]) + 5*(x[1]-0.5)*(x[2] - 0.5) + 5*(x[3] - x[4]) + 5*(x[3]-0.5)*(x[4] - 0.5)
	sigmoid(psi)
}

# Simulate one random observation, with 10-dim x.
wx_linear_interaction(runif(10))

# Simulate random x's to calibrate level.
test2 = do(50)*{
	x = runif(10)
	w = wx_linear_interaction(x)
	0.25*x[5] + w*f2(t_grid) + (1-w)*f1(t_grid)
}

# Plots hazard functions.
plot(t_grid, test1[1,], ylim=c(0,1), type='l',
	col=rgb(0.1,0.1,0.1,0.1))
for(i in 1:50) {
	lines(t_grid, test2[i,],col=rgb(0.1,0.1,0.1,0.1))
}

#---------------------------------------------------------
# 13. Nonlinear with interaction.
#---------------------------------------------------------

# Same four basic shapes, but more complicated function of x.

## Non-linear single-index with an interaction
# Involves maxes, bc they are hard to learn.
wx_nonlinear_interaction = function(x) {
	psi = 5*max(x[1:2]) - 5*max(x[3:4])
	sigmoid(psi)
}

# Simulate one random observation, with 10-dim x.
wx_nonlinear_interaction(runif(10))

# Simulate random x's to calibrate level.
test3 = do(50)*{
	x = runif(10)
	w = wx_nonlinear_interaction(x)
	0.25*x[5] + w*f2(t_grid) + (1-w)*f1(t_grid)
}

# Plots hazard functions.
plot(t_grid, test1[1,], ylim=c(0,1), type='l',
	col=rgb(0.1,0.1,0.1,0.1))
for(i in 1:50) {
	lines(t_grid, test3[i,],col=rgb(0.1,0.1,0.1,0.1))
}

#----------------------------------------------------------
##### Now how to deflate these in a way to expand them to hazard function.

makeHazard = function(N, P, t_grid= seq(0.1, 1, by=0.1), wx_fun = "wx_linear"){
   
   # create the design matrix and calculate the weight function
   X = matrix(runif(N*P, 0, 1), nrow=N)
   w = apply(X,1, paste(wx_fun)) # alpha on board.
   
   print(w)
   
   # calculate the hazard function at each time point
   hazard = foreach(i=1:N, .combine='rbind') %do% {
      0.25*X[i,5] + w[i]*f2(t_grid) + (1-w[i])*f1(t_grid)
   }
   
   # overall survival probability for each row.  knock back hazards by multiplicative factor of m.
   # may be different for each of the three functions.  want roughly 50%.
   if(wx_fun=='wx_linear'){ m = 0.3 }
   if(wx_fun=='wx_linear_interaction'){ m = 0.3 }    
   if(wx_fun=='wx_nonlinear_interaction'){ m = 0.3 } 
   
   # Avoid rescaling.
   m=1
   
   hazard_rescaled = m*hazard
   surv_prob = apply(hazard_rescaled, 1, function(x) {prod(1-x)})
   hist(surv_prob)
   mean(surv_prob)
   
   # Calculate whether event actually happened at each time.
   # Move through columns of hazard_rescaled matrix.
   events = matrix(NA, nrow=N, ncol=P)
   events[,1] = rbinom(N, size=1, prob=hazard_rescaled[,1])

   for(t in 2:ncol(events)){
      no_event_yet = which(events[,t-1]!=1 & !is.na(events[,t-1]))
      events[no_event_yet ,t] = rbinom(length(no_event_yet), 
                                              size=1, 
                                              prob=hazard_rescaled[no_event_yet,t])
   }
   
   # Format as long dataframe.  Each row is an obs, with (y_i, t_i, haz_it, x_i vector.)
   colnames(hazard_rescaled) = t_grid
   colnames(events) = t_grid
   
   df = as.data.frame(matrix(0, nrow=N * length(t_grid), ncol=15))
   colnames(df) = c('id','t','y','haz_fun','overall_surv_prob', paste0('x',1:10))
   
   row = 1
   
   for(i in 1:N){
      
      for(t in 1:length(t_grid)){
         
         temp_id = i
         temp_t = t #t_grid[t]
         temp_y = events[i,t]
         temp_haz = hazard_rescaled[i,t]
         temp_surv = surv_prob[i]
         temp_X = X[i,]
         
         df[row,] = c(temp_id, temp_t, temp_y, temp_haz, temp_surv, temp_X)
         row = row + 1
      }
      
   }
   
   # Add test/train indicators for each ID, with a 80/20 split.
   tt = sample(1:N, size=N*.8, replace=F)
   df$testtrain = ifelse(df$id %in% tt, 'train','test')
   
   return(df)

}


#----------------------------------------------------------
# Create datasets with N=1000.
#----------------------------------------------------------
N = 1000
P = 10

df_lin = makeHazard(N, P, wx_fun = 'wx_linear')
df_lin_int = makeHazard(N, P, wx_fun = 'wx_linear_interaction')
df_nonlin_int = makeHazard(N, P, wx_fun = 'wx_nonlinear_interaction')

write.csv(df_lin, './data/simulations/df_linear.csv', row.names=F)
write.csv(df_lin_int, './data/simulations/df_linear_interaction.csv', row.names=F)
write.csv(df_nonlin_int, './data/simulations/df_nonlinear_interaction.csv', row.names=F)

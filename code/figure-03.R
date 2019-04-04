##########################################################
# Plot hazard functions for the simulated data study.
##########################################################

library(mosaic)
library(foreach)
library(tidyverse)
library(ggplot2)
library(gridExtra)

source('./code/helper-functions/ggtheme-publication.R')
theme_set(theme_bw(base_size=16, base_family='Helvetica'))

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

# Define three wx functions:
wx_linear = function(x) {
	psi =  5*(x[1] - x[2] + x[3] - x[4])
	sigmoid(psi)
}

wx_linear_interaction = function(x) {
   psi = 5*(x[1] - x[2]) + 5*(x[1]-0.5)*(x[2] - 0.5) + 5*(x[3] - x[4]) + 5*(x[3]-0.5)*(x[4] - 0.5)
   sigmoid(psi)
}

wx_nonlinear_interaction = function(x) {
   psi = 5*max(x[1:2]) - 5*max(x[3:4])
   sigmoid(psi)
}


# Simulate random x's to calibrate level.
n=15
t_grid = seq(0,1,by=0.01)

test1 = do(n)*{
	x = runif(10)
	w = wx_linear(x)
	0.25*x[5] + w*f2(t_grid) + (1-w)*f1(t_grid) #.75 from board is absorbed in f1, f2.
}

test2 = do(n)*{
	x = runif(10)
	w = wx_linear_interaction(x)
	0.25*x[5] + w*f2(t_grid) + (1-w)*f1(t_grid)
}

test3 = do(n)*{
   x = runif(10)
   w = wx_nonlinear_interaction(x)
   0.25*x[5] + w*f2(t_grid) + (1-w)*f1(t_grid)
}

# Plots hazard functions.
df = as.data.frame(matrix(0,nrow=n*length(t_grid)*3, ncol=4))
colnames(df) = c('method','id','t', 'haz')

df$method = rep(c('linear','linear_interaction','nonlinear_interaction'), each=n*length(t_grid))
df$label = rep(c('Linear','Linear with interaction','Nonlinear with interaction'), each=n*length(t_grid))
df$id = rep(1:n, times=3*length(t_grid))
df$t = rep(rep(t_grid, each=n), times=3)
df$haz = c(unlist(test1), unlist(test2), unlist(test3))

plt = ggplot(df, aes(x=t, y=haz, colour=factor(id))) + 
   geom_line(size=.5) + 
   facet_wrap(~label,ncol=3) + 
   scale_colour_manual(values=rep('black',100), guide=F) + 
   theme_Publication() + 
   labs(x='t', y='Hazard function')
plt
ggsave('./output-figures/figure-03.pdf', plt, height=3, width=9, dpi=300)   

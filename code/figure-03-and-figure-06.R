# For creating Figures 5 and 6.
# Note: requires tsbart package to reproduce regular BART code.
# tsbart package can be obtained by request from Jared S. Murray.

rm(list=ls())

#===================================================================
# User Inputs
#===================================================================

# tsbart filepath. Example is:
tsbart_path = getwd()

# Simulation name to execute.  
filename = "/data/supp-simulations/tsb-sumofcosines-simdata-T8-p4-n100.csv"  

#===================================================================
# Workspace prep
#===================================================================

# Install tsbart package from local directory.
library(devtools)
install_github("jestarling/tsbart")
library(mosaic)
library(tidyverse)
library(gridExtra)
library(tsbart)
library(fastbart)

# ggplot settings
theme_set(theme_bw(base_size=13, base_family='Helvetica'))

# Output directories.
out.fig = paste0(getwd(),'/output-figures/')


##################################################################################
### Loop through each simulation file, if an individual name is not specified.
##################################################################################

# Optimize expected number of crossings over user-defined grid.
exp_cross = NULL

#===================================================================
# Read data.
#===================================================================

# Read csv file.
sim = read.csv(paste0(getwd(),'/',filename))
sim = read.csv("/Users/jennstarling/Stats/UTAustin/Research/tsbart-analysis/data/supp-simulations/tsb-sumofcosines-simdata-T8-p4-n100.csv")

# Set up train and test data sets.
train = sim %>% filter(train=='train')
test = sim %>% filter(train=='test')

# Extract values from data frame.
y = train$y;       y_pred = test$y        # Vector of responses.
ti = train$ti;     ti_pred = test$ti      # Time points for each obs.
fx = train$fx;     fx_pred = test$fx      # True underlying function, without epsilon noise.

xcols = which(substr(colnames(train),1,1)=="x")
xx = train[xcols]; x_pred = test[xcols]     # Matrix of predicted covariates.


#############################################################################################
###   1. tsBART.
#############################################################################################

#=====================================================================
#=== Evaluate optimal expected number of crossings by selecting
#=== one with lowest WAIC.
#=====================================================================

ecross_candidates = seq(.5,5,by=.5)
# Evaluate optimal number of crossings.
ecrossTune = tuneEcross(ecross_candidates,
                        y=y, tgt=ti, tpred=ti_pred,
                        x=xx, xpred=x_pred, nburn=50, nsim=50, ntree=200)

# Set expected number of crossings.
exp_cross = ecrossTune$ecross_opt

# Output figure for appendix.
ggsave(paste0(out.fig,"figure-06.pdf"),
       ecrossTune$waic_plot, height=5, width=7, dpi=300)

#===================================================================
# Fit BART model.
#===================================================================
ntree = 200; nburn = 2000; nsim = 2000

fit = tsbart(y=y, tgt=ti, tpred=ti_pred, x=xx, xpred=x_pred, 
      nburn, nsim, ntree, ecross=exp_cross, use_fscale=T)

#===================================================================
# Plot out of sample curves over time.
# Note: not included in plot function bc need faceting for x's.
#===================================================================

# Set up data frames for plotting out of sample curves over time.
ggdf_pred = cbind.data.frame(
   'id' = test$id,
   't' = ti_pred,
   'y' = y_pred,
   'fx' = test$fx,
   'pmean' = apply(fit$mcmcdraws_oos, 2,function(x) mean(x)),
   'lb' = apply(fit$mcmcdraws_oos, 2,function(x) quantile(x,.025)),
   'ub' = apply(fit$mcmcdraws_oos, 2,function(x) quantile(x,.975)),
   x_pred
)

ggdf_pred$label = paste0('Patient ',ggdf_pred$id)

# Select a few IDs to include.
id_sel = c(5,41,79,98)

# Set plotting limits.
ylims = c(min(ggdf_pred$lb)-.5, max(ggdf_pred$ub)+.5)

# Plot out-of-sample fit for each curve over time, for combos of covariates x.
p1 <- ggplot(ggdf_pred %>% filter(id %in% id_sel)) +
   geom_line(aes(y = pmean, x = t, colour='pmean'), linetype=1, size=1, stat="identity") +
   geom_line(aes(y=fx, x=t, colour='fx'), linetype=2, size=1, stat="identity",alpha=.85) +
   geom_ribbon(aes(ymin=lb,ymax=ub, x=t),alpha=0.2) +
   facet_wrap(~label,ncol=4) + 
   scale_colour_manual(values=c("firebrick3","grey10","grey40"), 
                       labels=c('True function value','Estimated function value'),
                       name=' ') +
   scale_fill_manual(values="grey80") +
   guides(colour=FALSE) +
   coord_cartesian(ylim=ylims) +
   labs(x='',
        y='y = f(x,t)',
        title = 'tsBART out of sample fit') 


#############################################################################################
###   2. Vanilla BART.
#############################################################################################

# Re-create x and xpred to include the time covariate.
xx = cbind.data.frame(xx,ti)
x_pred = cbind.data.frame(x_pred,'ti' = ti_pred)

# Calibrate BART's error variance a la CGM 2010 (method 2)
df = data.frame(xx, y)
lmf = lm(y~.,df)
sighat = sigma(lmf) 

# Hyperparameters
nu = 3
sigq = .9
qchi = qchisq(1.0-sigq,nu)
lambda = (sighat*sighat*qchi)/nu

# Update cutpoints to include t.
cutpoints = makeCutpoints(xx)

# Fit BART model.
fit_v = bartRcppClean(y_ = y, x_ = t(xx), # obsvervations must be in *columns*
                                xpred_ = t(x_pred), #Dpred,#[,-2,drop=F],
                                xinfo_list = cutpoints,
                                nburn, nsim, ntree,
                                lambda, nu, kfac=2,
                                paste0(getwd(),"/vanilla-bart-trees.txt"), # saves trees - can delete later
                                RJ=FALSE)

#===================================================================
# Plot vanilla out of sample curves over time.
#===================================================================

# Set up data frames for plotting out of sample curves over time.
ggdf_pred_v = cbind.data.frame(
   't' = ti_pred,
   'id' = test$id,
   'y' = y_pred,
   'fx' = sim %>% filter(train=='test') %>% select(fx),
   'pmean' = colMeans(fit_v$postpred),
   'lb' = apply(fit_v$postpred, 2,function(x) quantile(x,.025)),
   'ub' = apply(fit_v$postpred, 2,function(x) quantile(x,.975)),
   x_pred
)

ggdf_pred_v$label = paste0('Patient ',ggdf_pred_v$id)


# Plot out-of-sample fit for each curve over time, for combos of covariates x.
p2 <- ggplot(ggdf_pred_v %>% filter(id %in% id_sel)) +
   geom_line(aes(y = pmean, x = t, colour='pmean'), linetype=1, size=1, stat="identity") +
   geom_line(aes(y=fx, x=t, colour='fx'), linetype=2, size=1, stat="identity") +
   geom_ribbon(aes(ymin=lb,ymax=ub, x=t),alpha=0.2) +
   scale_colour_manual(values=c("firebrick3","grey10","grey40"), 
                       labels=c('True function value','Estimated function value'),
                       name=' ') +
   scale_fill_manual(values="grey80") +
   facet_wrap(~label, ncol=4) +
   coord_cartesian(ylim=ylims) +
   labs(x='Time (t)',
        y='y = f(x,t)',
        title = 'BART out of sample fit') +
   guides(colour=FALSE)

p2

#############################################################################################
###   3. Out-of-sample Comparison
#############################################################################################

# OOS RMSE
sqrt(mean((ggdf_pred$y - ggdf_pred$pmean)^2))
sqrt(mean((ggdf_pred_v$y - ggdf_pred_v$pmean)^2))

# Out of sample log-likelihoods for posterior means.
ll_oos = dnorm(ggdf_pred$y, ggdf_pred$pmean, 1, T)
ll_oos_v = dnorm(ggdf_pred_v$y, ggdf_pred_v$pmean, 1, T)

# Set up ggplot data frame.
lldf = cbind.data.frame('ll' = c(ll_oos, ll_oos_v),
                        'model' = c(rep('tsb',length(ll_oos)), 
                                    rep('vb', length(ll_oos_v))),
                        't' = rep(ggdf_pred$t,times=2))
# Summarize by time point.
lldf_t = lldf %>% dplyr::group_by(t,model) %>% summarise('ll'=mean(ll))

# Overall means.
tsb_meanll = round(mean(ll_oos),2)
vb_meanll = round(mean(ll_oos_v),2)

# Create plot.
llplt = ggplot(lldf_t, aes(x=ll, fill=model)) + 
   geom_density(alpha=.6,na.rm=T) +
   scale_fill_manual(values = c('grey10','grey60'), 
                     labels=c(paste0('tsBART (',tsb_meanll,')'), paste0('BART (',vb_meanll,')')), 
                     name=' ') +
   labs(x = 'Log likelihood', y = 'Density') +
   theme(strip.text.x = element_text(size = 13, colour = "black", angle = 0)) +
   theme(legend.position='top') +
   scale_x_continuous(expand=c(0,0))


################################################################
# Assemble grid of plots.
################################################################

# Define layout.
lay <- rbind(c(1,3),c(2,3))

# Set up panel.
panel = grid.arrange(grobs = list(p1,p2,llplt), 
             layout_matrix = lay, widths=c(2,1))

# Save plot.
ggsave(paste0(out.fig,"figure-03.pdf"), 
       plot(panel), height=6, width=12, dpi=300)



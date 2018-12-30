#**********************************************************************
# This script takes the raw National Health Center data and processes
# it.  Processing steps include: selecting variables for analysis,
# creating new variables where necessary, and removing missing and 
# NA rows.
#
# Reads a raw file called BIRTH_DEATH_DATA_2004-2006.csv, so the user
# should use this filename for their raw file, which they can assemble 
# from the publicly available data,
#**********************************************************************

# Clean workspace.
rm(list=ls())

#=======================================================================
# User inputs: file path to BIRTH_DEATH_DATA_2004-2006.csv/
#=======================================================================

# Example is filled in - set your own path.
raw_data_path = paste0(getwd(),'/Data/BIRTH_DEATH_DATA_2004-2006.csv')
output_path = paste0(getwd(),'/Data/stillbirth-clean.csv')

#=======================================================================
# Workspace prep.
#=======================================================================

library(data.table)  # For fast csv read/write.
library(dplyr)       # For re-leveling factor.
library(Rcpp)
library(RcppArmadillo)

# Load helper functions.
source(paste0(getwd(),'/code/00_stillbirth-data-processing/stillbirth-dataformat-functions.R'))
sourceCpp(paste0(getwd(),'/code/00_stillbirth-data-processing/stillbirth-dataformat-functions-cpp.cpp'))

#=======================================================================
#===   Read data   =====================================================
#=======================================================================

dropcols = c(278:339)
data = as.data.frame(fread(raw_data_path, sep=',', header=T, drop=dropcols))
rm(dropcols)

#=======================================================================
# Covariate setup: gestational age
#=======================================================================

# Limit analysis to pregnancies delivered between 34 and 42 weeks' gestation.
data = data[which(data$gestation_detail >= 34 & data$gestation_detail <= 42),]

# Gestational age.
ga = data$gestation_detail

#=======================================================================
# Set up response variables
#=======================================================================

# Neonatal death.
nd = ifelse(!is.na(data$neon_death), data$neon_death, ifelse(data$neonataldeath==2,0,1) )
table(nd, useNA='always')
#xtabs(~nd + ga, addNA=T)

# Fetal death.
fd = data$fetal_death
table(fd, useNA='always')
#xtabs(~fd + ga, addNA=T)

# Live birth info.
livebirths = data$livebirth
#table(livebirths, useNA='always')


#=======================================================================
#===   Set up covariates   =============================================
#=======================================================================

#-----------------------------------------------------------------------
# Whether labor was induced.
#-----------------------------------------------------------------------

induce = ifelse(data$induction_labr_u==2,0,data$induction_labr_u)

#-----------------------------------------------------------------------
# High-risk indicator covariates. 
#  (Anemia, Cardiac Disease, Lung Disease, Diabetes Mellitus,
#  Hemoglobinopathy, Chronic HTN, Renal Disease, Rh Sens)
#-----------------------------------------------------------------------

# Anemia
anemia = pmin(data$anemia_u, data$ab_anemia_u,  na.rm=TRUE)

# Cardac disease
cardiac = data$cardiac_u

# Diabetes
diabetes = data$diabetes_u

# Hemoglobinopathy
hemoglob = data$hemoglobinopathy_u

# Hypertension.   (Adding prepregnancy_ht does not add any new info.)
htn = data$chronic_ht_u

# Lung disease.
lung = data$lung_disease_u

# Renal disease.
renal = data$renal_disease_u

# Rh sensitization.
rhsens = data$rh_sensitizatn_u

# High risk indicator.
hirisk = pmin(anemia, cardiac, diabetes, hemoglob, htn, lung, renal, rhsens, na.rm=T)

#-----------------------------------------------------------------------
# Weight Gain
#-----------------------------------------------------------------------
wtGain = data$wtgain

#-----------------------------------------------------------------------
# Mother age
#-----------------------------------------------------------------------
momAge = data$mother_age

#----------------------------------------------------------------------
# Create variable for finer discretization of age variable.  
# Age buckets of every 5 years.  No separate bucket for 45-50 and then 50+, as no nd or fd cases in 50+.
#----------------------------------------------------------------------
agegrp = ifelse(momAge >= 45, "45+",
                ifelse(momAge >= 40, "40-44",
                       ifelse(momAge >= 35, "35-39",
                              ifelse(momAge >= 30, "30-34",
                                     ifelse(momAge >=25, "25-29",
                                            ifelse(momAge >= 20, "20-24", "<20"))))))


#-----------------------------------------------------------------------
# Mother ethnicity
#-----------------------------------------------------------------------
momEthnicity = ifelse(data$mother_racehisp_origin %in% c(1:5), 'Hispanic',
                      ifelse(data$mother_racehisp_origin==6, 'White NonHisp',
                             ifelse(data$mother_racehisp_origin==7, 'Black NonHisp','Other')))

momEthnicity = factor(momEthnicity, levels=c('White NonHisp', 'Black NonHisp', 'Hispanic', 'Other'))


#-----------------------------------------------------------------------
# Prior pregnancies.
#-----------------------------------------------------------------------
# Prior live births.
priorLive = ifelse(data$prior_live > 0, 1, 0)

# Prior stillbirths.
priorDead = ifelse(data$prior_dead > 0, 1, 0)

# Prior terminated pregnancies.
priorTerms = ifelse(data$prior_terminations > 0, 1, 0)

#-----------------------------------------------------------------------
# Parity
#-----------------------------------------------------------------------
multiparous = ifelse(data$total_birth_order==99, 'Unknown', ifelse(data$total_birth_order==1, 0,1))

#-----------------------------------------------------------------------
# Male infant
#-----------------------------------------------------------------------
male = as.factor(data$sex_infant)

#-----------------------------------------------------------------------
# Birth Weight
#-----------------------------------------------------------------------
birthwt = data$birth_wt

#-----------------------------------------------------------------------
# Birth Year (for test/train splitting) (or year of fetal death)
#-----------------------------------------------------------------------
year = ifelse(!is.na(data$birth_yr), data$birth_yr, data$fetal_death_yr)



#=======================================================================
#===   Assemble new data set   =========================================
#=======================================================================

# convert to factors.
nd = as.factor(nd)
fd = as.factor(fd)
induce = as.factor(induce)
hirisk = as.factor(hirisk)
anemia = as.factor(anemia)
cardiac = as.factor(cardiac)
diabetes = as.factor(diabetes)
hemoglob = as.factor(hemoglob)
htn = as.factor(htn)
lung = as.factor(lung)
renal = as.factor(renal)
rhsens = as.factor(rhsens)


# Assemble full data set.
ob = data.frame('nd' = nd,
                'fd' = fd,
                'livebirth' = livebirths,
                'ga' = ga,
                'induce' = induce,
                'hirisk' = hirisk,
                'priorLive' = priorLive,
                'priorDead' = priorDead,
                'priorTerms' = priorTerms,
                'wtGain' = wtGain,
                'momAge' = momAge,
                'agegrp' = agegrp,
                'momEthnicity' = momEthnicity,
                'anemia' = anemia,
                'cardiac' = cardiac,
                'diabetes' = diabetes,
                'hemoglob' = hemoglob,
                'htn' = htn,
                'lung' = lung,
                'renal' = renal,
                'rhsens' = rhsens,
                'multiparous' = multiparous,
                'male' = male,
                'birthwt' = birthwt,
                'year' = year)


# Re-level so that 0=no instead of 2=no.
levels(ob$hirisk)[2] = "0"
levels(ob$anemia)[2] = "0"
levels(ob$cardiac)[2] = "0"
levels(ob$diabetes)[2] = "0"
levels(ob$hemoglob)[2] = "0"
levels(ob$htn)[2] = "0"
levels(ob$lung)[2] = "0"
levels(ob$renal)[2] = "0"
levels(ob$rhsens)[2] = "0"
levels(ob$male)[2] = "0"

summary(ob)


#==========================================================================
# Omit NA values and missing values, and save a clean copy of the data set.
#==========================================================================

# Remove NAs.
ob.clean = na.omit(ob)

# Remove high risk unknowns, induce.  Drop unused levels.
ob.clean = drop_unknown(ob.clean, colnames(ob.clean))

# Output full data set, including NA data.
#fwrite(ob.clean, 'Data/data_infantdeath.csv')

######################################################################
# Part 2: Add quantile variables, clean up levels, etc.
######################################################################

data = fix_levels(ob.clean)

# Add 'weeks' continuous variable to data set.
data$weeks = as.numeric(data[,"ga"]) - as.numeric(baseline_week)

# Add 'other significant risks' - all hirisk excl diabetes and htn.
data$otherRisk = ifelse(data$anemia ==1 | data$cardiac ==1 | data$hemoglob ==1 | data$lung ==1 | data$renal ==1 | data$rhsens ==1, 1, 0)

# Add categorical variable for neither-diabetes-htn-both.
data$diab_htn = factor(ifelse(data$diabetes==0 & data$htn==0, 'Neither',
                              ifelse(data$diabetes==1 & data$htn==0, 'Diabetes',
                                     ifelse(data$diabetes==0 & data$htn==1, 'Htn', 'Both'))),
                       levels = c('Neither','Diabetes','Htn','Both'))

# Add id variable.
data$id = 1:nrow(data)

#----------------------------------------------------------------------
# Create quintiles for weight gain (conditional on gestational age).
#----------------------------------------------------------------------
prbs = c(0,.10,.25,.75,.90,1)

wtgainQ_out = conditQuantile(data, var='wtGain', prbs, condit='ga')

data$wtgainQ = as.factor(wtgainQ_out$qtiles)
data$wtgainQ_vals = as.factor(wtgainQ_out$vals)


#----------------------------------------------------------------------
# Quintiles for birth weight, at each gestational age level.  (conditional on gestational age: birthwt | ga)
#----------------------------------------------------------------------

prbs = c(0,.10,.25,.75,.90,1)

birthwtQ_out = conditQuantile(data, var='birthwt', prbs, condit='ga')

data$birthwtQ = as.factor(birthwtQ_out$qtiles)
data$birthwtQ_vals = as.factor(birthwtQ_out$vals)

#----------------------------------------------------------------------

# Output cleaned dataset.
fwrite(data, output_path)
print('Neonatal death csv complete.')

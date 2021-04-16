## Test model simulations with standardized inputs and outputs for model fitting with grid search to Variability experiment 1 run in January 2021
## 15/04/2021 

rm(list=ls())

#----------------------------------------------------------
# set dirs and load required scripts
#----------------------------------------------------------
# # When sourcing this script
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)

# When running in console
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

#----------------------------------------------------------
source ("../utilities/utilities.R")
source ("../utilities/joint.utilities.R")
source ("../utilities/condprob.utilities.R")
source ("../utilities/cgm.R")
source ("../models/mutsampler.R")
source ("../models/param_variability.R")

source ("../models/BMS.R")
source ("../models/beta_variability_model.R")




#----------------------------------------------------------
# Testing inputs
#----------------------------------------------------------


# Inputs to for simulations for fitting to Variability experiment 1 data, we want inputs to be causal strengths and baserates
# we used common cause network, with causal strengths 75%, and baserates for all variables 50%. (implying strength of background causes for effects is 20%)

var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['Y', 'X1'] = ms['Y','X2'] = .75
bs = c (.5, .5, .5)
  


#----------------------------------------------------------
# Beta Var model
#----------------------------------------------------------

conc <- 10
beta.sims <- betavar.pdf.samps(ms, bs, conc, nSamps=10000)
beta.sims



#----------------------------------------------------------
# BMS
#----------------------------------------------------------

# input here should be ms, bs, chainlength, betavar, nchains
# output should be dataframe with response distributions
# we employ a middle step where chains are generated, as these can be reused for different chainlengths and betavars



chains <- genchainsMSclean(len=50, nchains=10000, ms=ms, bs=bs) #generate chains

meanlen <- 10
betavar <- 1


BMS.sims <- chainstoresps.BMS(chains, meanlen=meanlen, betavar=betavar) #takes a bit
BMS.sims

#----------------------------------------------------------
# Par Var model
#----------------------------------------------------------



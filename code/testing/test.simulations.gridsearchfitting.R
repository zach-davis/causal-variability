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
source ("../fitting/fitting_utilities.R")

source ("../models/BMS.R")
source ("../models/beta_variability_model.R")


load("../fitting/datexp1_clean.RData") 
#datlong is Df with cleaned data

#----------------------------------------------------------
# Testing inputs
#----------------------------------------------------------


# Inputs to simulations for fitting to Variability experiment 1 data, we want inputs to be causal strengths and baserates
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


plot(smoothkernel(beta.sims[,13]))
plot(density(beta.sims[,13]), main=colnames(beta.sims)[13]) # test plot


PDA.SLL(smoothkernel(beta.sims[,13]), c(.5,.5))


#----------------------------------------------------------
# BMS
#----------------------------------------------------------

# input here should be ms, bs, chainlength, betavar, nchains
# output should be dataframe with response distributions
# we employ a middle step where chains are generated, as these can be reused for different chainlengths and betavars

meanlen <- 10
betavar <- 1

BMS.sims <- genchainsMSclean(len=50, nchains=100, ms=ms, bs=bs) %>% #generate chains
  chainstoresps.BMS(meanlen=meanlen, betavar=betavar) #takes a bit
  
plot(density(BMS.sims[,4]), main=colnames(BMS.sims)[4]) # test plot




#----------------------------------------------------------
# Par Var model
#----------------------------------------------------------
# bs = c (.5, .5, 0) #cant do this, all ahve to be nonzero?
parvar.sims <- param_variability(ms = ms,
                       bs = bs,
                       ms_conc = 10,
                       bs_conc = 10,
                       nSamples = 100)


parvar.sims <- Parvarrespdistr(parvar.sims$chainjoints)

plot(density(Parvarrespdistr(parvar.sims$chainjoints)[,11]))


#----------------------------------------------------------
# Guess fnc (to be moved to other script)
#----------------------------------------------------------
#needs to count how many responses there are, and then add a percentage of .5 resps to that 
guess_funIK <- function(simresps, guess_prob) {
  #simsresps is vector of simulated responses
  nresps <- length(simresps) #nr of resps
  nrguesses <- (nresps*guess_prob)/(1-guess_prob) #nr guesses that need to be added (so that guesses are guess_prob of total)
  resps <- append(simresps, rep(.5, nrguesses))
  return(resps)
}

respvec1 <- Parvarrespdistr(parvar.sims$chainjoints)[,11]
plot(density(guess_funIK(respvec1, .60)))



#----------------------------------------------------------
# Stimulus encoding noise model
#----------------------------------------------------------
#a probability of misreading a variable value in the stimulus. 
# 2 possible responses if 1 variable given
# 4 possible responses if 2 variables in stimulus.


#inputs: Option1, to do for all inferences at once inputs:
# inference type (maybe multiple vars?)
# probability of misreading a variable value
# number of responses to generate

#this is formatting consistent with condprob functions, needed for regex wrangling
datlong$inftype2 <- factor(datlong$inftype, 
                           levels = c("P(Xi|Y=1,Xj=1)", "P(Xi|Y=1)", "P(Xi|Y=1,Xj=0)", "P(Y|Xi=1,Xj=1)", "P(Y|Xi=1)", "P(Y|Xi=1,Xj=0)"),
                           labels = c("X1|Y==1 & X2==1", "X1|Y==1", "X1|Y==1 & X2==0", "Y|X1==1 & X2==1", "Y|X1==1", "Y|X1==1 & X2==0"))


# Should work on every inference type!?
SimEncodeNoise <- function(ms, bs, inferencetype, misprob, nSamps=100) {
  #inferentype needs to be of format 'X1|Y==0 & X2==0', ie 
  #misprob is probability of misreading a value
  normresps <- normativeresps(ms, bs)
  correctans <- normresps[inferencetype] 
  resps <- rep(NA, nSamps)
  
  #here with regex get change the possible inference types??
  
  
  
  
  
  
  if (inftype==0){
    resps <- rep(correctans, nSamps)
  }
  
}

#hardcoded to only work on datlong data set experiment 1 variability. 
SimEncodeHard <- function(ms, bs, inferencetype, misprob, nSamps=100) {
  #inferentype needs to be of format 'X1|Y==0 & X2==0'
  #misprob is probability of misreading a value
  normresps <- normativeresps(ms, bs)
  resps <- rep(NA, nSamps)
  
  if (inftype=="P(Xi|Y=0,Xj=1"){
    
  }
  
  if (inftype=="P(Xi|Y=0,Xj=1"){
    
  }
  
}




#----------------------------------------------------------
# Stimulus coherence guess model
#----------------------------------------------------------



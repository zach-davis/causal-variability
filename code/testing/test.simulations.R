## Test model simulations with standardized inputs and outputs

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
source ("../utilities/cgm.R")
source ("../models/mutsampler.R")
source ("../models/param_variability.R")

source ("../models/BMS.R")
source ("../models/beta_variability_model.R")


#######################
# TODO
# set up outputs to include all input parameters !!!!!!!! Parameter var model still needs to be done.


#----------------------------------------------------------
# Testing inputs
#----------------------------------------------------------

# input joint directly (IK naming convention: states 1-8 go from 111-000)
testjoint <- data.frame(         X1=c(1,1,1,1,0,0,0,0),
                                 Y=c(1,1,0,0,1,1,0,0),
                                 X2=c(1,0,1,0,1,0,1,0),
                                 p=c(9,3,1,3,3,1,3,9)/32)  
rownames(testjoint)<- c("111", "110", "101", "100", "011", "010", "001", "000")
testjoint


# input using baserates and causal strengths
var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['X1','Y'] = ms['Y','X2'] = .5
bs = c (.5, .75, .75)
testjoint = joint.cgm.generic3(ms, bs) # joint.cgm.generic3 follows IK convention, starting at 111, and output lists with joint, ms, bs, normresps
testjoint





#----------------------------------------------------------
# Beta Var model
#----------------------------------------------------------

respdistr.betavar <- betavariability.sim(joint=testjoint, concentration=100, nSamps=100) #immediately outputs inferences and not joint, since the joint used is normative.


#----------------------------------------------------------
# Par Var model
#----------------------------------------------------------

a <- param_variability(ms=ms, 
                       bs=bs, 
                       ms_conc=.0,
                       bs_conc=.0,
                       nSamples=500)
a$meanjoint #THIS IS INCORRECT, SOMETHING GOING WRONG? THIS SHOULD BE THE NORMATIVE JOINT with noise=0, or at least somewhat resembling it? its output though are assymmetric even if input is not.

respdistr.parvar <- genrespdistr(a, betavar=0) #when using genrespdistr on param_variability model BETAVAR NEEDS TO BE ZERO!!!!!

#----------------------------------------------------------
# ZD mutation sampler
#----------------------------------------------------------

test.neighbors = neighbors.of.joint (testjoint$joint)
testmsjoint1 <- mutsampler2(joint=testjoint, chainLen=10, nChains=15, neighbors=test.neighbors)
respdistr.ZDMS <- genrespdistr(testmsjoint1, betavar=0) #betavar=1 is uniform prior
respdistr.ZDMS


#----------------------------------------------------------
# IK mutation sampler
#----------------------------------------------------------


res <- genchainsMSpoislen(meanlen=10, nchains=12, bias=0.5, joint=testjoint) #generate chains with lengths according poisson distr
#res <- genchainsMS(len=10, nchains=15, bias=0.5, joint=testjoint) #generate chains
testmsjoint2 <- genjoint(res) #calculate mean joint, and joint per chain, and vector chainlengths
respdistr.IKMS <- genrespdistr(testmsjoint2, betavar=0) # function  going from chainjoints + chainlens to response distributions.
respdistr.IKMS






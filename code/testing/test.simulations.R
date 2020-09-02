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
bs = c (.5, .25, .25)
testjoint = joint.cgm.generic2(ms, bs) # joint.cgm.generic2 follows IK convention, starting at 111
testjoint

# normative inferences
normresps <- list()
normresps[[1]]<-condproballBay(testjoint, 0, 10) #chainlength (3rd arg) is irrelevant when betavar(2nd arg)=0
normresps <- preddistrBay(1, normresps, infnames=0) #need to make this in DF then apply preddistrBay to obtain same order as model predictions
normresps

#----------------------------------------------------------
# Beta Var model
#----------------------------------------------------------

respdistr.betavar <- betavariability.sim(joint=testjoint, concentration=100, nSamps=1000) #immediately outputs inferences and not joint, since the joint used is normative.


#----------------------------------------------------------
# Par Var model
#----------------------------------------------------------

a <- param_variability2(ms=ms, 
                       bs=bs, 
                       ms_sd=.01,
                       bs_sd=.01,
                       nSamples=500)
a$meanjoint #THIS IS INCORRECT, SOMETHING GOING WRONG? THIS SHOULD BE THE NORMATIVE JOINT with noise=0, or at least somewhat resembling it? its output though are assymmetric even if input is not.

respdistr.parvar <- genrespdistr(a, betavar=0) #when using genrespdistr on param_variability model BETAVAR NEEDS TO BE ZERO!!!!!

#----------------------------------------------------------
# ZD mutation sampler
#----------------------------------------------------------

test.neighbors = neighbors.of.joint (testjoint)
testmsjoint1 <- mutsampler2(joint=testjoint, chainLen=10, nChains=15, neighbors=test.neighbors)
respdistr.ZDMS <- genrespdistr(testmsjoint1, betavar=0) #betavar=1 is uniform prior
respdistr.ZDMS


#----------------------------------------------------------
# IK mutation sampler
#----------------------------------------------------------

#res <- genchainsMSpoislen(meanlen=10, nchains=12, bias=0.5, joint=testjoint) #generate chains with lengths according poisson distr
res <- genchainsMS(len=10, nchains=15, bias=0.5, joint=testjoint) #generate chains
testmsjoint2 <- genjoint(res) #calculate mean joint, and joint per chain, and vector chainlengths
respdistr.IKMS <- genrespdistr(testmsjoint2, betavar=0) # function  going from chainjoints + chainlens to response distributions.
respdistr.IKMS



#----------------------------------------------------------
# Fnc to generate plots for each inference based on two sets of response distributions to compare
#----------------------------------------------------------

library(ggplot2)

# fnc generates 27 density plots to compare simulated response distributions.
compare.resp.distr.plots <- function(respdistr1, respdistr2, normresps){
        #function to input two response inference distributions and plot all densities to compare
        plots.compare.allinfs <- list()
        for (inf in 1:27){
                plots.compare.allinfs[[inf]] <- local({
                        inf <- inf
                        distr1 <- respdistr1[,inf]
                        distr2 <- respdistr2[,inf]
                        normresp <- normresps[[inf]]
                        titleinf <- names(respdistr1)[inf]
                        ks.stat <- ks.test(distr1, distr2)$statistic
                        ks.pval <- ks.test(distr1, distr2)$p.value
                        p1 <- ggplot() + 
                                geom_density(aes(x=distr1, y=..density..))+
                                geom_density(aes(x=distr2, y=..density..), colour='red')+
                                geom_vline(aes(xintercept=normresp))+
                                xlim(c(0,1))+
                                ggtitle(paste("response distribution", titleinf), subtitle = paste('KS stat', ks.stat, 'pvalue', ks.pval))
                        print(p1)
                })
        }
}

#----------------------------------------------------------
# Comparison ZD and IK mutation samplers simulated response distributions
#----------------------------------------------------------

joint <- testjoint
test.neighbors = neighbors.of.joint (testjoint)
chainlen <- 10
nchains <- 10000

testmsjoint <- mutsampler2(joint=testjoint, chainLen=chainlen, nChains=nchains, neighbors=test.neighbors)
respdistr.ZDMS <- genrespdistr(testmsjoint, betavar=1) 

res <- genchainsMS(len=chainlen, nchains=nchains, bias=0.5, joint=testjoint) 
testmsjoint2 <- genjoint(res)
respdistr.IKMS <- genrespdistr(testmsjoint2, betavar=1) 

normresps <- list()
normresps[[1]]<-condproballBay(testjoint, 0, 10) #chainlength (3rd arg) is irrelevant when betavar(2nd arg)=0
normresps <- preddistrBay(1, normresps, infnames=0) #need to make this in DF then apply preddistrBay to obtain same order as model predictions

compare.resp.distr.plots(respdistr.ZDMS, respdistr.IKMS, normresps)


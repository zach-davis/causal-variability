#### Create plots for the simulations


#----------------------------------------------------------
# set dirs and load required scripts
#----------------------------------------------------------
# # When sourcing this script
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)

# When running in console
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

plotdir <- file.path(dirname(dirname(getwd())), "plots")  #ugly, make standard directory setup

#----------------------------------------------------------
source ("../utilities/utilities.R")
source ("../utilities/joint.utilities.R")
source ("../utilities/cgm.R")
source ("../models/mutsampler.R")
source ("../models/param_variability.R")

source ("../models/BMS.R")
source ("../models/beta_variability_model.R")

source ("../testing/fncs.sim.plots.R")




#----------------------------------------------------------
# Generate joint and set variables as input for simulations and plots
#----------------------------------------------------------

var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['X1','Y'] = ms['Y', 'X2'] = .5
bs = c (.25, .25, .25)
testjoint1 = joint.cgm.generic3(ms, bs) 
test.neighbors1 = neighbors.of.joint (testjoint$joint)


var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['X1','Y'] = ms['Y', 'X2'] = .5
bs = c (.65, .65, .65)
testjoint2 = joint.cgm.generic3(ms, bs) 
test.neighbors2 = neighbors.of.joint (testjoint2$joint)


chainlen <- 13
meanlen <- 13
nchains <- 5000


#----------------------------------------------------------
# Run simulations
#----------------------------------------------------------

# simulations done: ZD vs IK sampler, no vs uniform prior chain ms=.5 bs=.25
# simulations to do: low bs vs high bs chain. low ms vs high ms chain.


# # Simulate responses ZD sampler
# set.seed(1234)
# testmsjoint <- mutsampler2(joint=testjoint, chainLen=chainlen, nChains=nchains, neighbors=test.neighbors)
# respdistr.ZDMS <- genrespdistr(testmsjoint, betavar=0) 

# Simulate responses IK sampler
set.seed(1234)
res <- genchainsMSpoislen(meanlen=chainlen, nchains=nchains, bias=0.5, joint=testjoint1) 
testmsjoint2 <- genjoint(res)
respdistr.IKMS <- genrespdistr(testmsjoint2, betavar=0)

set.seed(1234)
res <- genchainsMSpoislen(meanlen=chainlen, nchains=nchains, bias=0.5, joint=testjoint2)
testmsjoint2 <- genjoint(res)
respdistr.IKMS2 <- genrespdistr(testmsjoint2, betavar=1)

# simulate responses Beta var model

#################### NEED TO ADAPT INPUT PLOTS TO NOT try to plot meanjoint
respdistr.betavar <- betavariability.sim(joint=testjoint1, concentration=100, nSamps=5000) 

# simulate responses par var model
respdistr.parvar <- param_variability(ms=ms, 
                       bs=bs, 
                       ms_conc=.01,
                       bs_conc=.01,
                       nSamples=500)

# Correct responses 
normresps <- list()
normresps[[1]]<-condproballBay(testjoint, 0, 10) #chainlength (3rd arg) is irrelevant when betavar(2nd arg)=0
normresps <- preddistrBay(1, normresps, infnames=0) #need to make this in DF then apply preddistrBay to obtain same order as model predictions



#----------------------------------------------------------
# Generate plots comparison ZD and IK mutation samplers simulated response distributions
#----------------------------------------------------------


testplot <- compare.resp.distr.plots(respdistr.IKMS, respdistr.parvar) 

pdf(file.path(plotdir, "comparison_noprior_uniformprior_poislen.pdf"),onefile = TRUE)
testplot
dev.off()


#----------------------------------------------------------
# Generate plots to visualize markov violations and explaining away
#----------------------------------------------------------

testplot2 <- plot.mrkv(respdistr.IKMS, nms.chain.mrkv) # of a single simulated resp distribution


testplot3 <- compare.plot.mrkv(respdistr.IKMS, respdistr.parvar, nms.chain.mrkv) # comparing them from two different simulations
pdf(file.path(plotdir, "comparison_noprior_uniformprior_poislen_markov.pdf"),onefile = TRUE)
testplot3
dev.off()

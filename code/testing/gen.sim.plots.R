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

source ("../testing/fncs.sim.plots.R")




#----------------------------------------------------------
# Generate joint and set variables as input for simulations and plots
#----------------------------------------------------------

var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['X1','Y'] = ms['X2', 'Y'] = .5
bs = c (.5, .5, .5)
testjoint1 = joint.cgm.generic3(ms, bs) 
test.neighbors1 = neighbors.of.joint (testjoint1$joint)


var.names = c ('X1','Y','X2')
ms2 = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['X1','Y'] = ms['X2', 'Y'] = .5
bs2 = c (.15, .15, .15)
testjoint2 = joint.cgm.generic3(ms2, bs2) 
test.neighbors2 = neighbors.of.joint (testjoint2$joint)


chainlen <- 13
meanlen <- 13
nchains <- 5000


#----------------------------------------------------------
# Run simulations
#----------------------------------------------------------


# Simulate responses ZD sampler
set.seed(1234)
testmsjoint <- mutsampler2(joint=testjoint, chainLen=chainlen, nChains=nchains, neighbors=test.neighbors)
respdistr.ZDMS <- genrespdistr(testmsjoint, betavar=0)

# Simulate responses IK sampler
set.seed(1234)
res <- genchainsMSpoislen(meanlen=chainlen, nchains=nchains, bias=0.5, joint=testjoint1) 
testmsjoint1 <- genjoint(res)
respdistr.IKMS <- genrespdistr(testmsjoint1, betavar=1)

set.seed(1234)
res <- genchainsMSpoislen(meanlen=chainlen, nchains=nchains, bias=0.5, joint=testjoint2)
testmsjoint2 <- genjoint(res)
respdistr.IKMS2 <- genrespdistr(testmsjoint2, betavar=1)

# simulate responses Beta var model
respdistr.betavar <- betavariability.sim(joint=testjoint1, concentration=10, nSamps=5000) 

# simulate responses par var model
a <- param_variability(ms=ms, 
                       bs=bs, 
                       ms_conc=10,
                       bs_conc=10,
                       nSamples=5000)
respdistr.parvar <- genrespdistr(a, betavar=0)




#----------------------------------------------------------
# Generate plots comparison ZD and IK mutation samplers simulated response distributions
#----------------------------------------------------------


testplot <- compare.resp.distr.plots(respdistr.IKMS, respdistr.IKMS2)

pdf(file.path(plotdir, "testplot.pdf"),onefile = TRUE)
testplot
dev.off()


#----------------------------------------------------------
# Generate plots to visualize markov violations and explaining away
#----------------------------------------------------------

testplot2 <- plot.mrkv(respdistr.IKMS, nms.chain.mrkv) # of a single simulated resp distribution


testplot3 <- compare.plot.mrkv(respdistr.IKMS, respdistr.IKMS2, nms.comeff.expaw) # comparing them from two different simulations
pdf(file.path(plotdir, "testplot_markov.pdf"),onefile = TRUE)
testplot3
dev.off()

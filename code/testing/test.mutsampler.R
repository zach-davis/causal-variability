#----------------------------------------------------------
# Unit testing
#     Comparing sampling to normative Bayes nets
#----------------------------------------------------------

library (doMC)  # Multicore support 
library(foreach)
no.cpus = 4
registerDoMC(no.cpus)  # No. of cores to use

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source ("../utilities/utilities.R")
source ("../utilities/joint.utilities.R")
source ("../utilities/cgm.R")
source ("../models/mutsampler.R")

print.screening.off = function (j) {
  p.z.x1.y1 = v.p (j, 'z', gx=c (x=1, y=1))
  p.z.x0.y1 = v.p (j, 'z', gx=c (x=0, y=1))
  cat ('p(z=1|y=1,x=1) = ', p.z.x1.y1, '\n')
  cat ('p(z=1|y=1,x=0) = ', p.z.x0.y1, '\n')
  cat ('Markov violation: p(z=1|y=1,x=1) - p(z=1|y=1,x=0) = ', p.z.x1.y1 - p.z.x0.y1, '\n')
}

test.joints = function (chain.joint, sampled.joint) {
  cat ('In normative joint, Y screens Z off from X. In fact, the Markov condition is satisfied:\n')
  print.screening.off (chain.joint)
  cat ('In contrast, let\'s see what happens in the sampled joint:\n')
  print.screening.off (sampled.joint)
}

# CHAIN NETWORK TESTS

# Create joint distribution for chain network: X-> Y -> Z
#   Causal strengths (or "powers") = .75
#   Base rate of X = .50
#   Background causes of Y and Z = 1/3
var.names = c ('x','y','z')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['x','y'] = ms['y','z'] = .75
bs = c (.5, 1/3, 1/3)
chain.joint = joint.cgm.generic (ms, bs)
chain.neighbors = neighbors.of.joint (chain.joint)

cat ('*** TESTS OF A CHAIN NETWORK: X-> Y -> Z.\n')
cat ('*** Test of mutation sampling with actual sampling:\n')
chainLen = 10
sampled.joint = mutsampler (chain.joint, chainLen, nChains=1000, neighbors=chain.neighbors, sample.in.parallel=T)
test.joints (chain.joint, sampled.joint)

cat ('*** Test of mutation sampling with the expected joint (computed without sampling):\n')
exact.joint = mutsampler.exact (chain.joint, chainLen, chain.neighbors)
test.joints (chain.joint, exact.joint)

cat ('*** Test of mutation sampling with expected joint, where chain length is sampled from a Poisson distribution :\n')
poisson.joint = mutsampler.poisson (chain.joint, chainLen, chain.neighbors)
test.joints (exact.joint, poisson.joint)

cat ('*** Size of Markov violation decreases as chain length increases::\n')
for (chainLen in c (5, 15, 45)) {
  cat ('Chain length =', chainLen, '\n')
  exact.joint = mutsampler.exact (chain.joint, chainLen, chain.neighbors)
  print.screening.off (exact.joint)
}


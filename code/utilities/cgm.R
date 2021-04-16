this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source ("utilities.R")
source ("joint.utilities.R")
source ("condprob.utilities.R")


#-------------------------------------------------------------------------------------------
# Generic routine that return joint distributions for causal graphical models (cgms).
#
#   ms - square matrix of causal strengths (rows must be named, no cycles)
#   bs - background causes of each variable (appear in same order as in ms)
#
#-------------------------------------------------------------------------------------------
joint.cgm.generic = function (ms, bs) {
  vs = rownames (ms)
  j = joint.skeleton (vs)
  j$p = 1
  j = joint.cgm.generic.driver (j, 1:length (vs), ms, bs)
  return (j)
}

joint.cgm.generic2 = function (ms, bs) { #same as above, except reverses order of rows, so returns joint in IK convention (ie starting at 111)
  j = joint.cgm.generic(ms, bs)
  j = j[nrow(j):1,]
  return (j)
}

joint.cgm.generic3 = function (ms, bs) { #same as above (joint.cgm.generic2), except now includes ms,bs, and normative inferences as outputs
  j = joint.cgm.generic(ms, bs)
  j = j[nrow(j):1,]
  normresps <- list()
  normresps[[1]]<-condproballBay(j, 0, 10)
  normresps <- preddistrBay(1, normresps, infnames=0)
  res <- list(joint=j, ms=ms, bs=bs, normresps=normresps)
  return (res)
}


#-------------------------------------------------------------------------------------------
# Utility routine.
# Compute joint probabilities.
# j$p is the probability of the root nodes. child.vis are the child (non-root) nodes.
#-------------------------------------------------------------------------------------------
joint.cgm.generic.driver = function (j, child.vis, ms, bs) {
  # Probability that an effect is in state e given n independent causes.
  e.given.icn  = function (e, ms, b, ics) {
    e0.given.icn = function () { 
      if (length (ms) == 0)
        1 - b
      else 
        (1 - b) * prod (mapply (function (m, ic) {(1 - m)**ic}, ms, ics)) 
    }  
    p0 = e0.given.icn ()
    p = ifelse (e == 1, 1 - p0, p0)
    return (p)
  }
  ps.for.vi = function (vi) {
    pi = ms[,vi] > 0  # Parents of variable vi
    ms = ms[,vi][pi]  # Strength of the links from the parents
    b = bs[vi]        # Strength of variable vi's background causes.
    # Compute the probability of variable i in each network state.
    apply (j[,joint.vs (j)], 1, function (x)
      e.given.icn (x[vi], ms, b, x[pi])
    )      
  }
  stopifnot (nrow (ms) == length (bs), nrow (ms) == ncol (ms), ms >= 0, ms <= 1, all (diag (ms) == 0), bs >= 0, bs <= 1)
  # Compute ps for each child variable in each network state.
  child.ps = sapply (child.vis, ps.for.vi)
  # For each network state, multiply the probability of the root nodes and the probabilities of the children. 
  j$p = j$p * apply (child.ps, 1, prod)
  stopifnot (equal.within.tol (sum (j$p), 1))
  j = joint.w.rownames (joint.normalized (j))
  return (j)
}


#############################################################################
# Generic routine for arbitrary network of (independent, generative) relations.
#   ms - square matrix of causal strengths (rows must be named, no cycles)
#   bs - background causes of each variable (appear in same order as in ms)
#############################################################################
joint.cgm.generic = function (ms, bs) {
    vs = rownames (ms); j = joint.skeleton (vs); j$p = 1
    joint.cgm.generic.driver (j, 1:length (vs), ms, bs)
}

# Utility routine.
joint.cgm.generic.driver = function (j, child.vis, ms, bs) {
    ps.for.vi = function (i) {
        pi = ms[,i] > 0; ms = ms[,i][pi]; b = bs[i] # Parents of i, and the strengths of the links.
        apply (j[,joint.vs (j)], 1, function (x) e.given.icn (x[i], ms, b, x[pi]))      
    }
    stopifnot (nrow (ms) == length (bs), nrow (ms) == ncol (ms), ms >= 0, ms <= 1, all (diag (ms) == 0), bs >= 0, bs <= 1)
    # Compute ps for each child variable in each state.
    child.ps = laply (child.vis, ps.for.vi, .parallel = TRUE)
    # Multiply by the probability of the root nodes.
    j$p = j$p * apply (child.ps, 2, prod); stopifnot (equal.within.tol (sum (j$p), 1))
    joint.w.rownames (j)
}


# setting up sample call
sample_network1 <- t(array(c(0, 0, 0,
                             1, 0, 0,
                             1, 1, 0),
                           dim=c(3,3)))
vars <- c('x','y','z')
rownames(sample_network1) <- vars
colnames(sample_network1) <- vars

# sample call
joint.cgm.generic(ms=sample_network1, bs=rep(.1,3))

library(boot)
library(neldermead)
library(optimbase)
library(plyr)

cgm.data.likelihood = function (cgm.joint.fun, cgm.params, data, with.nas = FALSE) {
  j = do.call (cgm.joint.fun, cgm.params)
  sample.likelihood (j, data, with.nas)
}

cgm.data.logl = function (cgm.joint.fun, cgm.params, data, with.nas = FALSE) {
  j = do.call (cgm.joint.fun, cgm.params)
  sample.loglk (j, data, with.nas)
}

cgm.visible.data.likelihood = function (cgm.joint.fun, cgm.params, visible.vars, data) {
  j = do.call (cgm.joint.fun, cgm.params)
  j = integrate.over.others (j, visible.vars)
  sample.likelihood (j, data)
}

cgm.visible.data.logl = function (cgm.joint.fun, cgm.params, visible.vars, data) {
  log (cgm.visible.data.likelihood (cgm.joint.fun, cgm.params, visible.vars, data))
}

most.likely.cgm.params = function (cgm.joint.fun, data, with.nas = FALSE, logit.params = TRUE) { 
  transform.ps   = function (ps) mapply (function (p, logit.param) if (logit.param) logit     (p) else p, ps, logit.params)
  untransform.ps = function (ps) mapply (function (p, logit.param) if (logit.param) inv.logit (p) else p, ps, logit.params)
  outfun = function (x, optimValues, state) {
    if ((optimValues$iteration %% 20) == 0) { #cat("**** Iteration ", optimValues$iteration, " -- ", format (optimValues$fval, nsmall = 8), "\n", x, "\n") 
    }   
  }
  TolX = .0001; TolFun = .0001; #TolX = 100; TolFun = 100
  opt = optimset (OutputFcn = outfun, TolX = TolX, TolFun = TolFun, MaxIter = 10000000, MaxFunEvals = 10000000)
  eval.params = function (ps) { 
    ps = untransform.ps (ps)
    if (any (ps[logit.params] < min.cgm.p) | any (ps[logit.params] > max.cgm.p)) err = Inf
    else {
      err = - cgm.data.logl (cgm.joint.fun, as.list (ps), data, with.nas) 
      if (is.na (err)) { print (err); print (ps); stopifnot (!is.na (err)) }
    }
    err
  }
  
  # Set up parameter transformations. 
  num.p = length (formals (cgm.joint.fun))
  if (length (logit.params) == 1) logit.params = rep (logit.params, num.p) else stopifnot (length (logit.params) == num.p)
  
  x = transform.ps (rep (.5, num.p))
  # Do the search (repeated 5 times).
  iter = 0
  for (i in 1:5) { sol = fminsearch (eval.params, x, opt); iter = iter + sol$output$iterations; x = as.vector (t (neldermead.get(sol, 'xopt'))) }
  if (sol$exitflag != 1) { print (sol); stop ("Fminsearch failure....") }  
  #cat ("iterations=", iter,"err=", sol$fval,"x=", x,"\n"); #print.eval.ds (ds, x); 
  untransform.ps (x)
}
  
update.cgm.params = function (d, true.params, sample.size) {
  for (p in length (d)) { update.beta.pdist (d [[p]], true.params [p], sample.size) }
  d
} 

learn.link.in.cgm = function (cgm.joint.fun, true.params, link.to.learn, sample.size) {
  j = do.call (cgm.joint.fun, as.list (true.params))
  s = sample.joint (j, sample.size);               
  visible.s = occluded.sample (s, link.to.learn)
  colnames (visible.s) = c ("c", "e", "count")
  ps = most.likely.cgm.params (joint.cgm.c1e1, visible.s); colnames (ps) = c ("c", "m", "b")
  ps
}

learn.ic2e1.in.cgm = function (cgm.joint.fun, true.params, ic2e1.vars, sample.size) {
  j = do.call (cgm.joint.fun, as.list (true.params))
  s = sample.joint (j, sample.size);               
  visible.s = integrate.sample.over.others (s, ic2e1.vars)
  colnames (visible.s) = c ("ca", "cb", "e", "count")
  most.likely.cgm.params (joint.cgm.ic2e1.full, visible.s); colnames (ps) = c ("ca", "cb", "ma", "mb", "b")
}

learn.link.in.uncertain.cgm = function (d, cgm.joint.fun, link.to.learn, cgm.sample.size, data.sample.size) {
  process.sample = function (s) { learn.link.in.cgm (cgm.joint.fun, as.list (s), link.to.learn, data.sample.size) }
  mvsample = mv.pdist.sample (d, cgm.sample.size) # Sample parameters from distribution(s) in d.
  learned.ps = aaply (mvsample, 1, process.sample, .parallel = TRUE)
  mean (as.data.frame (learned.ps)); names (mean.ps) = c ("c", "m", "b")
} 

learn.ic2e1.in.uncertain.cgm = function (d, cgm.joint.fun, link.to.learn, cgm.sample.size, data.sample.size) {
  process.sample = function (s) { learn.ic2e1.in.cgm (cgm.joint.fun, as.list (s), link.to.learn, data.sample.size) }
  mvsample = mv.pdist.sample (d, cgm.sample.size) # Sample parameters from distribution(s) in d.
  learned.ps = aaply (mvsample, 1, process.sample, .parallel = TRUE)
  mean.ps = mean (as.data.frame (learned.ps)); names (mean.ps) = c ("ica", "icb", "ma", "mb", "b")
  mean.ps
} 


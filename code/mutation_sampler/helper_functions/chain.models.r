#############################################################################
# Routines that return joint distributions for chain graphs.
# 
# Example of naming conventions: ab.cd..wx.yz 
#   - Directed links between A and B, and C and D.
#   - Undirected links between W and X, and Y and Z.
#
# Routines for graphs that include no directed links just call the appropriate 
# routine for creating an MRF.
# 
# Routines implement the  chain graph Markov property as described in Lauritzen & 
# Richardson (2002), "Chain graph models and their causal interpretations."
#############################################################################
library (boot)

cgm.c1e1.stats = function (c, m, b) { 
  j = joint.cgm.c1e1 (c, m, b); e = v.p (j, 'e')
  list (marg = mean (c (c, e)), marg.e = e, assoc = lod (j, 'c', 'e'))
}

#############################################################################
# GRAPHS WITH TWO VARIABLES.
#############################################################################

#############################################################################
# X=Y
#############################################################################
joint.cm..xy = function (c, m, b) { 
  stats = cgm.c1e1.stats (c, m, b)
  joint.mrf.xy (stats$marg, stats$marg, stats$assoc)
}

#############################################################################
# GRAPHS WITH THREE VARIABLES.
#############################################################################

#############################################################################
# X=C->E
#############################################################################
j.cm.ce..xc = joint.skeleton (c ('x','c','e'))

joint.cm.ce..xc = function (c, m, b, fixed.marginals = TRUE) {  
  p = function (x) { x.p (f.ce, x['e'], gx = x['c']) * x.p (f.xc, x[c('x','c')]) }
  f.xc = joint.cm..xy (c, m, b); names (f.xc) = c('x','c','p')
  f.ce = joint.cgm.c1e1 (v.p (f.xc, 'c'), m, b)
  j = j.cm.ce..xc; j$p = apply (j, MARGIN = 1, p)
  j 
}

#############################################################################
# X->E<-Y; X=Y
# A common effect structure in which the two causes are associated. 
#############################################################################
j.cm.xe.ye..xy = joint.skeleton (c ('x','y','e'))

joint.cm.xe.ye..xy = function (c, m, b, assoc.xy, fixed.marginals = TRUE) {  
  p = function (x) { x.p (f.xye, x['e'], gx = x[c('x','y')]) * x.p (f.xy, x[c('x','y')]) }
  f.xy = joint.mrf.xy (c, c, assoc.xy)
  f.xye = joint.cgm.ic2e1 (v.p (f.xy, 'x'), m, b); names (f.xye) = c ('x','y','e', 'p')
  j = j.cm.xe.ye..xy; j$p = apply (j, MARGIN = 1, p)
  j 
}

joint.cm.xe.ye..xy.via.e = function (e, m, b, assoc.xy) {  
  j.fun = function (x) joint.cm.xe.ye..xy (inv.logit (x), m, b, assoc.xy, F)
  eval.fun = function (x) { j = j.fun (x); err = abs (v.p (j, 'e') - e); err }
  #browser ()
  joint.mrf.generic (eval.fun, j.fun, df = 1)
}

#############################################################################
# C->E=X
#############################################################################
j.cm.ce..xe = joint.skeleton (c ('c','e','x'))

joint.cm.ce..xe = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.cex, x[c('e','x')], gx=x['c']) * var.p (x['c'], c) 
  f.ce = joint.cgm.c1e1 (c, m, b); e = v.p (f.ce, 'e')
  f.ex = joint.mrf.xy (e, e, lod (f.ce)); names (f.ex) = c('e','x',"p")
  f.cex = merge.factors (list (f.ce, f.ex), fixed.marginals)
  j = j.cm.ce..xe; j$p = apply (j, MARGIN = 1, p)
  j 
}

#############################################################################
# X<-C->Y; X=Y
# A common cause structure in which the two causes are independently associated. 
#############################################################################
j.cm.cx.cy..xy = joint.skeleton (c ('c','x','y'))

joint.cm.cx.cy..xy = function (c, m, b, assoc.xy, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.cxy, x[c('x','y')], gx=x['c']) * var.p (x['c'], c) 
  f.cx = f.cy = joint.cgm.c1e1 (c, m, b)
  names (f.cx) = c('c','x','p'); names (f.cy) = c('c','y','p')
  e = v.p (f.cx, 'x'); f.xy = joint.mrf.xy (e, e, assoc.xy)
  f.cxy = merge.factors (list (f.cx, f.cy, f.xy), fixed.marginals)
  j = j.cm.cx.cy..xy; j$p = apply (j, MARGIN = 1, p)
  j 
}

joint.cm.cx.cy..xy.via.e = function (e, m, b, assoc.xy) {  
  j.fun = function (x) joint.cm.cx.cy..xy (inv.logit (x), m, b, assoc.xy, F)
  eval.fun = function (x) { j = j.fun (x); err = abs (v.p (j, 'x') - e) + abs (v.p (j, 'y') - e); err }
  joint.mrf.generic (eval.fun, j.fun, df = 1)
}

#############################################################################
# X=Y=Z (a chain).
#############################################################################
joint.cm..xy.yz = function (c, m, b, fixed.marginals = TRUE) { 
  stats = cgm.c1e1.stats (c, m, b)
  joint.mrf.xy.yz (stats$marg, stats$assoc, fixed.marginals)
}

#############################################################################
# X=Y, Y=Z, X=Z (a 3-cycle).
#############################################################################
joint.cm..xy.yz.xz = function (c, m, b, fixed.marginals = TRUE) { 
  stats = cgm.c1e1.stats (c, m, b)
  joint.mrf.xyz (stats$marg, stats$assoc, fixed.marginals)
}

#############################################################################
# GRAPHS WITH FOUR VARIABLES.
#############################################################################

#############################################################################
# W=X=Y=Z (a chain).
#############################################################################
joint.cm..wx.xy.yz = function (c, m, b, fixed.marginals = TRUE) { 
  stats = cgm.c1e1.stats (c, m, b)
  joint.mrf.wx.xy.yz (stats$marg, stats$assoc, fixed.marginals)
}

#############################################################################
# W=X, X=Y, Y=Z, Z=W (a 4-cycle)
#############################################################################
joint.cm..wx.xy.yz.zw = function (c, m, b, fixed.marginals = TRUE) { 
  stats = cgm.c1e1.stats (c, m, b)
  joint.mrf.wx.xy.yz.zw (stats$marg, stats$assoc, fixed.marginals)
}

#############################################################################
# E->X, X=Y, Y=Z, X=Z
#############################################################################
j.cm.ex..xy.yz.xz = joint.skeleton (c ('e','x','y','z'))

joint.cm.ex..xy.yz.xz = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.exyz, x[c ('x','y','z')], gx = x['e']) * var.p (x['e'], c)
  f.ex  = joint.cgm.c1e1 (c, m, b); names (f.ex) = c('e','x','p')
  f.xyz = joint.mrf.xyz (v.p (f.ex, 'x'), lod (f.ex), fixed.marginals)
  f.exyz = merge.factors (list (f.ex, f.xyz), fixed.marginals)
  j = j.cm.ex..xy.yz.xz; j$p = apply (j, MARGIN = 1, p)
  j 
}

joint.cm.ex..xy.yz.xz.old = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) { x.p (f.xyz, x[c ('y','z')], gx = x['x']) * x.p (f.ex, x[c('e','x')]) }
  f.ex  = joint.cgm.c1e1 (c, m, b); names (f.ex) = c('e','x','p'); x = v.p (f.ex, 'x')
  f.xyz = joint.mrf.xyz (x, lod (f.ex), fixed.marginals)
  j = j.cm.ex..xy.yz.xz; j$p = apply (j[,joint.vs (j)], MARGIN = 1, p)
  j 
}

joint.cm.ex..xy.yz.xz.old2 = function (c, m, b, fixed.marginals = T) { 
  p = function (x) { comp.exyz = joint.conditionalized (comp.exyz, x['e']); x.p (comp.exyz, x) * var.p (x['e'], c) }
  stats = cgm.c1e1.stats (c, m, b)
  comp.exyz = joint.mrf.ext.xyz (stats$marg, stats$assoc, fixed.marginals)
  j = j.cm.ex..xy.yz.xz; j$p = apply (j, MARGIN = 1, p)
  j
}

joint.cm.ex..xy.yz.xz.moralized = function (c, m, b) { 
  p = function (x) x.p (f.exyz, x[c ('x','y','z')], gx = x['e']) * var.p (x['e'], c)
  f.exz = setNames (joint.cgm.ic2e1 (c, m, b), c ('e','z','x','p'))
  f.xyz = joint.cm..xy.yz.xz (c, m, b)
  f.exyz = merge.factors (list (f.exz, f.xyz), T)
  j = j.cm.ex..xy.yz.xz; j$p = apply (j, MARGIN = 1, p)
  j
}

# Computes joint via sampling.
joint.cm.ex..xy.yz.xz.sampled = function (c, m, b) { 
  p = function (r) { 
    if (r[1] == 0) comp.xyz = comp.xyz.0
    else comp.xyz = comp.xyz.1
    x.p (comp.xyz, r[-1]) * var.p (r[1], c) 
  }
  f.ex = f.xy = f.yz = f.zx = joint.cgm.c1e1 (c, m, b)
  names (f.ex) = c('e','x','p'); names (f.xy) = c('x','y','p'); names (f.yz) = c('y','z','p'); names (f.zx) = c('z','x','p')
  comp.xyz.0 = joint.from.sampled.factors (list (f.xy, f.yz, f.zx), f.ex, c(e=0))
  comp.xyz.1 = joint.from.sampled.factors (list (f.xy, f.yz, f.zx), f.ex, c(e=1))
  j = joint.skeleton (c ('e','x','y',"z")); j$p = apply (j, MARGIN = 1, p)
  j
}

query.cm.ex..xy.yz.xz = function (c, m, b) {
  query.cgm.ext.xyz.joint (joint.cm.ex..xy.yz.xz (c, m, b))
} 

#############################################################################
# A->X, X=Y, Y<-B
#############################################################################
j.cm.ax.by..xy = joint.skeleton (c ('a','b','x','y'))

joint.cm.ax.by..xy = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.abxy, x[c('x','y')], gx = x[c('a','b')]) * var.p (x['a'], c) * var.p (x['b'], c)
  f.ax = f.by = joint.cgm.c1e1 (c, m, b)
  names (f.ax) = c('a','x','p'); names (f.by) = c('b','y','p')
  e = v.p (f.ax, 'x'); f.xy = joint.mrf.xy (e, e, lod (f.ax))
  f.abxy = merge.factors (list (f.ax, f.by, f.xy), fixed.marginals)
  j = j.cm.ax.by..xy; j$p = apply (j, MARGIN = 1, p)
  browser ()
  j
}

#############################################################################
# A->E<-B; E=X
#############################################################################
j.cm.ae.be..ex = joint.skeleton (c ('a','b','e','x'))

joint.cm.ae.be..ex = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.abex, x[c('e','x')], gx=x[c('a','b')]) * var.p (x['a'], c) * var.p (x['b'], c) 
  f.abe = joint.cgm.ic2e1 (c, m, b); names (f.abe) = c('a','b','e','p'); e = v.p (f.abe, 'e')
  stats = cgm.c1e1.stats (c, m, b)
  f.ex = joint.mrf.xy (e, e, stats$assoc); names (f.ex) = c('e','x',"p")
  f.abex = merge.factors (list (f.abe, f.ex), fixed.marginals)
  j = j.cm.ae.be..ex; j$p = apply (j, MARGIN = 1, p)
  j 
}

#############################################################################
# X=Y, Y=Z, X=Z, Z->E
#############################################################################
j.cm.ze..xy.yz.xz = joint.skeleton (c ('x','y','z','e'))

joint.cm.ze..xy.yz.xz = function (c, m, b, fixed.marginals = TRUE) {  
  p = function (x) { x.p (f.ze, x['e'], gx = x['z']) * x.p (f.xyz, x[c('x','y','z')]) }
  stats = cgm.c1e1.stats (c, m, b)
  f.xyz = joint.mrf.xyz (stats$marg, stats$assoc, fixed.marginals)
  f.ze  = joint.cgm.c1e1 (stats$marg, m, b); names (f.ze) = c('z','e','p')
  j = j.cm.ze..xy.yz.xz; j$p = apply (j[,joint.vs (j)], MARGIN = 1, p)
  j 
}

#############################################################################
# GRAPHS WITH FIVE VARIABLES.
#############################################################################

#############################################################################
# A->X, X=Y, Y<-B; X->E<-B
#############################################################################
j.cm.ax.by.xe.ye..xy = joint.skeleton (c ('a','b','x','y','e'))

joint.cm.ax.by.xe.ye..xy = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) {
    e.given.ic2 (x['e'], m, m, b, x['x'], x['y']) * 
      x.p (f.abxy, x[c('x','y')], gx = x[c('a','b')]) * var.p (x['a'], c) * var.p (x['b'], c)
  }
  f.ax = f.by = joint.cgm.c1e1 (c, m, b)
  names (f.ax) = c('a','x','p'); names (f.by) = c('b','y','p')
  xy = v.p (f.ax, 'x'); f.xy = joint.mrf.xy (xy, xy, lod (f.ax))
  f.abxy = merge.factors (list (f.ax, f.by, f.xy), fixed.marginals)
  j = j.cm.ax.by.xe.ye..xy; j$p = apply (j, MARGIN = 1, p)
  j
}

#############################################################################
# A->X, X=Y=Z, Z<-B
#############################################################################
j.cm.ax.by..xy.yz = joint.skeleton (c ('a','b','x','y','z'))

joint.cm.ax.by..xy.yz = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.abxyz, x[c('x','y','z')], gx = x[c('a','b')]) * var.p (x['a'], c) * var.p (x['b'], c)
  f.ax = f.bz = joint.cgm.c1e1 (c, m, b)
  names (f.ax) = c('a','x','p'); names (f.bz) = c('b','z','p')
  f.xyz = joint.mrf.xyz (v.p (f.ax, 'x'), lod (f.ax), fixed.marginals)
  f.abxyz = merge.factors (list (f.ax, f.bz, f.xyz), fixed.marginals)
  j = j.cm.ax.by..xy.yz; j$p = apply (j, MARGIN = 1, p)
  j
}

#############################################################################
# E->W, W=X, X=Y, Y=Z, X=W
#############################################################################
j.cm.ex..wx.xy.yz.xz = joint.skeleton (c ('e','w','x','y','z'))

joint.cm.ex..wx.xy.yz.xz = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) x.p (f.ewxyz, x[c ('w','x','y','z')], gx = x['e']) * var.p (x['e'], c)
  f.ew  = joint.cgm.c1e1 (c, m, b); names (f.ew) = c('e','w','p')
  f.wxyz = joint.mrf.wx.xy.yz.zw (v.p (f.ew, 'w'), lod (f.ew), fixed.marginals)
  f.ewxyz = merge.factors (list (f.ew, f.wxyz), fixed.marginals)
  j = j.cm.ex..wx.xy.yz.xz; j$p = apply (j, MARGIN = 1, p)
  j 
}

#############################################################################
# GRAPHS WITH SEVEN VARIABLES.
#############################################################################

#############################################################################
# E-> U = V; E->W, W=X, X=Y, Y=Z, X=W
#############################################################################
j.cm.eu.ew..uv.wx.xy.yz.xz = joint.skeleton (c ('e','u','v','w','x','y','z'))

joint.cm.eu.ew..uv.wx.xy.yz.xz = function (c, m, b, fixed.marginals = TRUE) { 
  p = function (x) { 
    x.p (f.ewxyz, x[c ('w','x','y','z')], gx = x['e']) * 
      x.p (f.euv, x[c ('u','v')], gx = x['e']) * var.p (x['e'], c)
  }
  f.euv   = joint.cm.ce..xe (c, m, b, fixed.marginals); names (f.euv) = c('e','u','v','p');
  f.ewxyz = joint.cm.ex..wx.xy.yz.xz (c, m, b, fixed.marginals)
  j = j.cm.eu.ew..uv.wx.xy.yz.xz; j$p = apply (j, MARGIN = 1, p)
  j
}

#############################################################################
# COMPOSITE MODELS (Independent combinations of the above.)
#############################################################################

#############################################################################
# Cyc2 + CE Composite Model
#############################################################################
joint.cm.cyc2.c1e1 = function (c, m, b) {
  merge.independent.joints (list (joint.cm..xy (c, m, b), joint.cgm.c1e1 (c, m, b) ))
}

query.cm.cyc2.c1e1 = function (c, m, b) query.cgm.cyc2.c1e1.joint (joint.cm.cyc2.c1e1 (c, m, b))

query.cm.double.cyc2.c1e1 = function (c.w, m.w, b.w, c.s, m.s, b.s) {
  q.w = query.cm.cyc2.c1e1 (c.w, m.w, b.w); q.w$str = "w"
  q.s = query.cm.cyc2.c1e1 (c.s, m.s, b.s); q.s$str = "s"
  q = rbind (q.w, q.s); rownames(q) = interaction (q$qtype, q$str, q$q)  
  q
}

#############################################################################
# Cyc2 + Chain3 Composite Model
#############################################################################
joint.cm.cyc2.ch3 = function (c, m, b) {
  merge.independent.joints (list (joint.cm..xy (c, m, b), joint.cgm.ch3 (c, m, b)))
}

query.cm.cyc2.ch3 = function (c, m, b) query.cgm.cyc2.ch3.joint (joint.cm.cyc2.ch3 (c, m, b))

query.cm.double.cyc2.ch3 = function (c.w, m.w, b.w, c.s, m.s, b.s) {
  q.w = query.cm.cyc2.ch3 (c.w, m.w, b.w); q.w$str = "w"
  q.s = query.cm.cyc2.ch3 (c.s, m.s, b.s); q.s$str = "s"
  q = rbind (q.w, q.s); rownames(q) = interaction (q$qtype, q$str, q$q)  
  q
}

#############################################################################
# Cyc2 + Cyc3 Composite Model
#############################################################################
j.cm.cyc2.cyc3 = joint.skeleton (c ("ya", "yb", "za", "zb", "zc"))

joint.cm.cyc2.cyc3 = function (c, m, b) {
  cyc2 = joint.cm..xy (c, m, b);       names (cyc2) = c("ya", "yb",'p')
  cyc3 = joint.cm..xy.yz.xz (c, m, b); names (cyc3) = c("za", "zb", "zc",'p')
  merge.independent.joints (list (cyc2, cyc3))
}

q.cm.cyc2.cyc3 = joint.query.skeleton (j.cm.cyc2.cyc3, max.conditioners = 0)
query.cm.cyc2.cyc3.joint      = function (j) joint.query.driver (j, q.cm.cyc2.cyc3)
query.cm.cyc2.cyc3            = function (c, m, b) query.cm.cyc2.cyc3.joint (joint.cm.cyc2.cyc3 (c, m, b))

#############################################################################
# Double Cyc2 + Cyc3 Composite Model. 
############################################################## ###############
query.cm.double.cyc2.cyc3 = function (c.w, m.w, b.w, c.s, m.s, b.s) {
  q.w = query.cm.cyc2.cyc3 (c.w, m.w, b.w); q.w$str = "w"
  q.s = query.cm.cyc2.cyc3 (c.s, m.s, b.s); q.s$str = "s"
  q = rbind (q.w, q.s); rownames(q) = interaction (q$qtype, q$str, q$q)  
  q
}

query.cm.double.cyc2.cyc3 = function (c.w, m.w, b.w, c.s, m.s, b.s) {
  q.w = query.cm.cyc2.cyc3 (c.w, m.w, b.w); q.w$str = "w"
  q.s = query.cm.cyc2.cyc3 (c.s, m.s, b.s); q.s$str = "s"
  q = rbind (q.w, q.s); rownames(q) = interaction (q$qtype, q$str, q$q)  
  q
}

#############################################################################
# Double E->X=Y=Z Composite Model. 
############################################################## ###############
# Separate parameters for each pair of networks ("w" and "s").
query.cm.double.ex..xy.yz.xz = function (c.w, m.w, y.w, c.s, m.s, y.s) {
  q.w = query.cm.ex..xy.yz.xz (c.w, m.w, y.w); q.w$str = "w"
  q.s = query.cm.ex..xy.yz.xz (c.s, m.s, y.s); q.s$str = "s"
  q = rbind (q.w, q.s); rownames(q) = interaction (q$qtype, q$str, q$q)  
  q
}

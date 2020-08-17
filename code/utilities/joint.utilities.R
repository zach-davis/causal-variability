#----------------------------------------------------------
# Miscellaneous routines to query and manipulate a joint distribution.
# Joints are represented as R data frames. Variables are binary by default.
#
# An example of a joint with three binary variables (x, y, and z): 
#
#       x y z          p
#   000 0 0 0 0.22222222
#   001 0 0 1 0.11111111
#   010 0 1 0 0.02777778
#   011 0 1 1 0.13888889
#   100 1 0 0 0.05555556
#   101 1 0 1 0.02777778
#   110 1 1 0 0.06944444
#   111 1 1 1 0.34722222
#
# The "p" column represent the probability of that state. Thus, "p" cannot be a variable name.
# By convention, rownames are the concatenation of the variable states.
# By convention, rows are ordered with "0" states first and "1" states last.
#
# Variable naming conventions:
#    j  -- a joint distibution
#    x  -- a system state (a setting for each system variables)
#    xs -- a vector of system states
#    v  -- a variable name
#    vs -- a vector of variable names
#    gx -- a to-be-conditioned on system state
#----------------------------------------------------------

binary.states = c (0, 1)

#----------------------------------------------------------
# Create the skeleton of a joint (columns for variables but no p column). 
#----------------------------------------------------------
joint.skeleton = function (vs, states=binary.states) {
  stopifnot (!'p' %in% vs) # Column 'p' represents the probability, so can't be a variable.
  j = joint.w.rownames (xs.of.vs (vs, states))
  return (j)
}   

#----------------------------------------------------------
# Set row names as interaction of all the vars.
# (Allows joint to be indexed like an array.)  
#----------------------------------------------------------
joint.w.rownames = function (j) {
  rownames (j) = interaction (j[joint.vs (j)], sep='')
  return (j)
}

#----------------------------------------------------------
# Returns a data frame in which all combinations of the values 
#   of the variables in vs have been instantiated.
# Variables are binary by default.
# If states is a list, then its elements represent the states 
#   for each variable in vs.
#----------------------------------------------------------
xs.of.vs = function (vs, states=binary.states) { 
  if (is.list (states))
    stopifnot (length (states) == length (vs))
  else 
    states = rep (list (states), length (vs))
  stopifnot (all (is.character (vs)), all (sapply (states, function (s) length (s) > 0)))
  # Create data frame representing joint (flip columns to put rows in canonical order).
  j = data.frame (do.call (expand.grid, rev (states))[,length(vs):1])
  colnames (j) = vs
  return (j)
}

#----------------------------------------------------------
# Return names of vars in j.
#----------------------------------------------------------
joint.vs = function (j) { 
  return (setdiff (names (j), 'p')) # All column names except "p"
}

#----------------------------------------------------------
# Marginal probability of v = 1 for each v in vs (suitable for binary variables).
# If vs is NULL, then marginals for all variables in j are returned.
# Joint is first conditionalized on a non-null gx.
#----------------------------------------------------------
v.p = function (j, vs = NULL, gx = NULL)  { 
  v.p.1 = function (v)
    sum (j$p[j[,v] == 1])
  
  if (is.null (vs))
    vs = joint.vs (j) else assert.names (j, vs)
  if (!is.null (gx))
    j = joint.conditionalized (j, gx)
  p.of.vs.given.gx = 
    if (length (vs) == 1) 
      v.p.1 (vs)
    else 
      sapply (vs, v.p.1)
  return (p.of.vs.given.gx)
}   

#----------------------------------------------------------
# Joint conditionalized on x.
#----------------------------------------------------------
joint.conditionalized = function (j, x) {
  j$p [!joint.idx (j, x)] = 0
  j = joint.normalized (j)
  return (j)
}

#----------------------------------------------------------
# Return an index vector indicating which rows of the joint correspond to state x.
# X can be a vector or list. If list, it can be a (one row) data frame.
#----------------------------------------------------------
joint.idx = function (j, x) {
  assert.names (j, names (x))
  if (is.list (x)) { 
    stopifnot (all (lengths (x) == 1))
    if (class (x) == 'data.frame')
      x = unlist (x) # Must unlist for recycling to work.
  } 
  else 
    stopifnot (is.vector (x))
  idx = if  (length (x) == 1) 
    j[, names(x)] == x 
  else if (length (x) == 0)
    rep (TRUE, nrow (j)) 
  else
    apply (sapply (names (x), function (v) j[, v] == x[v]), 1, all)
  return (idx)
}   

#----------------------------------------------------------
# Check that joint has the expected variables.
#----------------------------------------------------------
assert.names = function (j, vs) {
  if (is.null (vs))
    stop ('assert.names: Unnamed state.')
  if (!all (vs %in% joint.vs (j))) 
    stop ('assert.names: Names (', vs, '), not all in joint (', joint.vs (j), ').')
}   






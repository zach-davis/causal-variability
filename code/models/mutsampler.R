library ("plyr")
library ("reshape")
library ('gtools')

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source ("../utilities/joint.utilities.R")
source ("../utilities/cgm.R")

#-------------------------------------------------------------------------------------------
#   mutsampler                                                                                       
#                                                                                             
#   Parameters:                                                   
#     - joint:           the joint distribution to sample over, DF with column p indicating probability, and nrow being the nr of possible states                                  
#     - chainLen:        length of one chain                                                       
#     - nChains:         number of chains to run                                       
#     - neighbors:       indicates which states are "neighbors" of each other.   
#                        if specified, controls which states are included in the proposal distribution.   
#                        if NULL, all other states are eligible as the next state.   
#     - starting.states: start sampling at these states
#                        by default are the two prototypes of the distribution
#                        if NULL start sampling at random state
#     - bias:               bias toward the first starting state (by default, the 111... prototype).               
#     - method:          sampling method: mutate (default), gibbs, gibbs1.
#     - q.ps             an alternative proposal distribution (mutation sampling only)
#                                                                                             
#   mutsampler will start at a some point in the joint, and jump to another stimulus type based 
#   on the MH rule. The proposal structure looks as follows, where every line drawn between   
#   nodes means that they are 1 switch away from each other (and thus accessible by our      
#   proposal structure)                                                                      
#                                                                                            
#               001 - 011                                                                  
#             /     X     \                                                                   
#         000 - 010   101 - 111                                                              
#             \     X     /                                                                  
#               100 - 110                                                                     
#                                                                                            
#   Having a small "chainLen" parameter implies a large bias, since there is less time to     
#   escape from a local region. nChains does not change the bias, just reduces the variance. 
#
#   The standard mutation sampling proposal distribution only proposes states that are
#   neighbors of the current state. The neighbors are chosen with equal probability. 
#   When the resulting joint is going to be used for a conditional probability query, 
#   this can result in the sampling of states that are not useful for the query. Parameter
#   q.ps can be used to bias the sampling towards only useful states. See conditional.proposal.distribution
#   below and Appendix F of Davis & Rehder, 2020.
#-------------------------------------------------------------------------------------------
mutsampler <- function (joint,                               # joint distribution
                        chainLen,                        # number of steps in a chain
                        nChains,                         # number of chains to run
                        neighbors=NULL,                  # use neighbors proposal distribution or not?
                        starting.states=c(nrow(joint), 1), # possible states to start at, prototypes by default
                        bias=1 / length(starting.states), 
                        method='mutate', 
                        q.ps=NULL,                       # An alternative proposal distribution.
                        sample.in.parallel=F) {
  # function for running a chain
  chain.runner = function (i) {
    mcmcWinner <- rep(1e-10, no.states) # to store data, 1e-10 is to avoid na's, shouldn't matter
    # sample a starting state
    current <- sample(starting.states, 1, prob=starting.state.dist)
    # run one chain
    for (t in 1:chainLen) {
      mcmcWinner[current] <- mcmcWinner[current] + 1  # counts the starting point
      current <- sample(1:no.states, 1, prob=transition.ps[current,])
    }
    # return normalized run
    return (mcmcWinner / sum(mcmcWinner))
  }
  
  stopifnot (nChains > 0 & chainLen > 0)
  no.states <- nrow(joint)
  
  # defining the starting state distribution
  if (is.null (starting.states)) {
    # Randomly starting state.
    starting.state.dist = 1:no.states
    starting.state.dist <- rep(1/no.states, no.states)
  }
  else {
    # Use caller specified starting states.
    starting.state.dist <- c(bias, 1-bias)
  }
  
  # Precompute transition probability of any state to any other state
  transition.ps = transition.ps.fun (joint, neighbors, method, q.ps)
  
  # Run chains, potentially in parallel
  mcmcWinner <- ldply (1:nChains, chain.runner, .parallel=sample.in.parallel)
  
  # Average over chains
  mcmcWinner <- apply (mcmcWinner, 2, mean)
  stopifnot (equal.within.tol (sum (mcmcWinner), 1))
  joint$p = mcmcWinner / sum (mcmcWinner)
  return(joint)
}
#-------------------------------------------------------------------------------------------
#   same as above mutsampler, but also outputs the joint distr of each chain and chainlengths in addition to mean joint of all chains
#   !!! Only works for set chainlength, not for varying chainlength based on poisson distr !!!
mutsampler2 <- function (joint,                               # joint distribution
                        chainLen,                        # number of steps in a chain
                        nChains,                         # number of chains to run
                        neighbors=NULL,                  # use neighbors proposal distribution or not?
                        starting.states=c(nrow(joint), 1), # possible states to start at, prototypes by default
                        bias=1 / length(starting.states), 
                        method='mutate', 
                        q.ps=NULL,                       # An alternative proposal distribution.
                        sample.in.parallel=F) {
  # function for running a chain
  chain.runner = function (i) {
    mcmcWinner <- rep(1e-10, no.states) # to store data, 1e-10 is to avoid na's, shouldn't matter
    # sample a starting state
    current <- sample(starting.states, 1, prob=starting.state.dist)
    # run one chain
    for (t in 1:chainLen) {
      mcmcWinner[current] <- mcmcWinner[current] + 1  # counts the starting point
      current <- sample(1:no.states, 1, prob=transition.ps[current,])
    }
    # return normalized run
    return (mcmcWinner / sum(mcmcWinner))
  }
  
  stopifnot (nChains > 0 & chainLen > 0)
  no.states <- nrow(joint)
  
  # defining the starting state distribution
  if (is.null (starting.states)) {
    # Randomly starting state.
    starting.state.dist = 1:no.states
    starting.state.dist <- rep(1/no.states, no.states)
  }
  else {
    # Use caller specified starting states.
    starting.state.dist <- c(bias, 1-bias)
  }
  
  # Precompute transition probability of any state to any other state
  transition.ps = transition.ps.fun (joint, neighbors, method, q.ps)
  
  # Run chains, potentially in parallel
  chainjoints <- ldply (1:nChains, chain.runner, .parallel=sample.in.parallel) #joint distr per chain (previously mcmcWinner)
  
  # Average over chains
  meanjointdistr <- apply (chainjoints, 2, mean)
  stopifnot (equal.within.tol (sum (meanjointdistr), 1))
  joint$p = meanjointdistr / sum (meanjointdistr)
  
  chainlens <- rep(chainLen, nChains)
  
  res <- list(meanjoint=joint, chainjoints=chainjoints, chainlens=chainlens)
  return(res)
}

#--------------------------------------------------------------
# Returns expected result of sampling without actually sampling. 
mutsampler.exact = function (joint, 
                             chainLen, 
                             neighbors=NULL, 
                             biased_prototype=T, 
                             prototypes=c (nrow (joint), 1), 
                             bias=1 / length (prototypes), 
                             method='mutate', 
                             q.ps=NULL) {
  
  chain.fun = function (ps, transition.ps) {
    ps = next.ps = matrix (ps, nrow=1)
    if (chainLen > 1) {
      # Run the chain.
      for (i in 2:ceiling (chainLen)) {
        next.ps = next.ps %*% transition.ps # Expected probabilities at next iteration.
        prev.ps = ps
        ps = ps + next.ps # Running sum of expected no. of visits to each state.
      }
      # Compute weighted average of chainLen and (chainLen - 1).
      w = chainLen%%1; if (w > 0) ps = w * ps + (1 - w) * prev.ps 
      ps = ps / sum (ps) # Normalize
    }
    ps
  }
  stopifnot (chainLen >= 1)
  if (chainLen == Inf) {
    return (joint)
  } else {
    joint = mutsampler.generic (joint, chainLen, neighbors, biased_prototype, prototypes, bias, method, q.ps, chain.fun=chain.fun)
    return(joint)
  }
}

#-------------------------------------------------------------------------------------------
#   Returns expected result of sampling without actually sampling.                                                
#                                                                                             
#   Lambda is parameter of a Poisson distribution over chain lengths. 
#   Since chain lengths are at least 2, so too must be lambda
#   (i.e., actual Poisson distribution parameter is lambda - 2).
#-------------------------------------------------------------------------------------------
mutsampler.poisson = function (joint, 
                               lambda, 
                               neighbors=NULL, 
                               biased_prototype=T, 
                               prototypes=c (nrow (joint), 1), 
                               bias=1 / length (prototypes), 
                               method='mutate', 
                               q.ps=NULL) {
  chain.fun = function (ps, transition.ps) {
    ps = matrix (ps, nrow=1)
    stopifnot (lambda >= 2)
    lambda = lambda - 2 # Distribution over possible chain lengths (2, 3, ..) 
    max.steps = qpois (.999, lambda) # Almost all of the probability mass.
    p.per.chain.length = dpois (0:max.steps, lambda)
    p.per.chain.length = p.per.chain.length / sum (p.per.chain.length) # Normalize
    ps.weighted = 0; expected.no.visits = ps
    for (i in 1:length (p.per.chain.length)) { 
      ps = ps %*% transition.ps
      # Contribution at chain length i + 1
      expected.no.visits = expected.no.visits + ps
      joint.ps.at.length = expected.no.visits / sum (expected.no.visits)
      ps.weighted = ps.weighted + joint.ps.at.length * p.per.chain.length[i] 
    }
    stopifnot (equal.within.tol (sum (ps.weighted, na.rm=T), 1))
    ps.weighted / sum (ps.weighted, na.rm=T) # Normalize, just in case.
  }
  
  mutsampler.generic (joint, lambda, neighbors, biased_prototype, prototypes, bias, method, q.ps, chain.fun=chain.fun)
}

#-------------------------------------------------------------------------------------------
#   Helper function that returns matrix of state transition probabilities.        
#   
#   params
#       q.ps = transition probabilities, if NULL equal probability of transitioning
#-------------------------------------------------------------------------------------------
transition.ps.fun = function (joint, 
                              neighbors=NULL, 
                              method='mutate', 
                              q.ps=NULL, 
                              transitions.in.parallel=F) {
  
  # Transition probabilities for the standard "neighbors" proposal distribution. 
  mutate.transitions = function (xi) {
    nis = neighbors [[xi]]
    # Each neighbor's Hastings ratio, times the probability that we'll sample that neighbor.
    ps.of.xi = rep (0, nrow (joint))
    ps.of.xi[nis] = sapply (nis, function (ni) { 
      min (c (joint$p[ni] * q.ps[ni, xi] / (joint$p[xi] * q.ps[xi, ni]), 1)) * q.ps[xi, ni]
    })
    stopifnot (!any (is.na (ps.of.xi)))
    # Rest of probability mass goes to current state.
    ps.of.xi[xi] = 1 - sum (ps.of.xi) 
    ps.of.xi
  }
  
  # Transition probabilities based on Gibbs sampling
  gibbs.transitions = function (xi, depth=length (joint.vs (joint))) {
    # Recurse down the tree of variables/values in the order specified by vis.
    ps1 = function (vis) {
      recurse = function (xi, vii, ps) {
        stopifnot (rownames (joint)[xi] == names (neighbors)[xi])
        vi = vis[vii] # Variable to apply Gibbs update.
        ni = neighbors[[xi]][vi] # New state if Gibbs variable changes.
        new.val.p = joint$p[ni] / (joint$p[ni] + joint$p[xi]) # Probability that it will change.
        if (vii < length (vis))
          ps = ps +
          recurse (ni, vii + 1, ps) * new.val.p +
          recurse (xi, vii + 1, ps) * (1 - new.val.p)
        else { # Recursion has bottomed-out.
          ps[ni] = new.val.p
          ps[xi] = 1 - new.val.p
        }
        ps
      }
      ps.of.xi = recurse (xi, 1, rep (0, nrow (joint)))
      ps.of.xi / sum (ps.of.xi) 
    }
    # All possible variable orders (but only go to depth specified by depth).
    vis = permutations (length (joint.vs (joint)), depth) 
    # Transition probabilities for each order.
    ps.of.xi = aaply (vis, 1, ps1)
    # Average over orders.
    apply (ps.of.xi, 2, mean)
  }
  
  gibbs1.transitions = function (xi) {
    gibbs.transitions (xi, depth=1)
  }
  
  nr = nrow (joint)
  # Each row xi has transition probabilities for state j[xi,] 
  joint = joint.normalized (joint)
  if (is.null (neighbors)) {
    # By default, all nodes (states) can be proposed from any state.
    neighbors = setNames (lapply (1:nr, function (xi) seq (1, nr)[-xi]), rownames (joint))
  }
  if (is.null (q.ps)) {
    # By default, proposal distribution is uniform over neighbors.
    q.ps = sapply (neighbors, 
                   function (nis) {
                     sapply (1:nr, 
                             function (ni) {
                               if (ni %in% nis) {
                                 1 / length (nis)
                               } else {
                                 0
                               }
                             })
                   })
  }
  else {
    stopifnot (method == 'mutate') # q.ps only applies for mutation sampler.
  }
  
  transition.funs = list (mutate=mutate.transitions, gibbs=gibbs.transitions, gibbs1=gibbs1.transitions)
  stopifnot (!is.null (method) & method %in% names (transition.funs))
  transition.ps = laply (1:nr, transition.funs[[method]], .parallel=transitions.in.parallel)
  return(transition.ps)
}

#-------------------------------------------------------------------------------------------
#  Helper function that, for each entry in joint 'joint', compute its neighbors (entries that differ by     
#  one binary variable).                                                        
neighbors.of.joint = function (joint) {
  neigbors.of.xi = function (xi) {
    flip1 = function (vi) { 
      substr (xi.rowname, vi, vi) = if (substr (xi.rowname, vi, vi) == '1') '0' else '1'; 
      xi.rowname 
    }
    xi.rowname = rownames (joint)[xi]
    # Get rownames of neighbors.
    ns.of.xi = sapply (1:length (joint.vs (joint)), flip1)
    # Convert to row numbers.
    sapply (ns.of.xi, function (n) seq (1, nrow (joint))[which (n == rownames (joint))])
  }
  joint = joint.w.rownames (joint)
  setNames (lapply (1:nrow (joint), neigbors.of.xi), rownames (joint))
}

#--------------------------------------------------------------
#  Used for both exact and poisson sampling. 
mutsampler.generic = function (joint, 
                               chainLen, 
                               neighbors=NULL, 
                               biased_prototype=TRUE, 
                               prototypes=c (nrow (joint), 1), 
                               bias=1 / length (prototypes), 
                               method='mutate', 
                               q.ps=NULL, 
                               chain.fun) {
  nr = nrow (joint)
  # Each row xi has transition probabilities for state j[xi,] 
  joint = joint.normalized (joint)
  transition.ps = transition.ps.fun (joint, neighbors, method, q.ps)
  
  # if we're dealing with a biased prototype...
  if (biased_prototype & length (prototypes) > 0) {
    # Start with probability mass on prototypes.
    if (is.character (prototypes[1])) {
      # Row names to row numbers.
      prototypes = sapply (prototypes, function (rn) which (rn == rownames (joint)))
    }
    
    ps = rep (NA, nr)
    ps[-prototypes] = 1e-10
    p.prototypes = 1 - sum (ps[-prototypes])
    ps[prototypes[1]] = bias * p.prototypes # By default, mass to all variables = 1. 
    p.prototypes = p.prototypes - ps[prototypes[1]]
    ps[prototypes[2:length (prototypes)]] = p.prototypes / (length (prototypes) - 1) # By default, mass to all variables = 0.
  }
  else {
    ps = rep (1 / nr, nr)
  }
  
  # Run the chain
  joint$p = as.vector (do.call (chain.fun, list (ps, transition.ps)))
  return(joint)
}

#--------------------------------------------------------------
# Returns a proposal distribution that biases sampling away from those system 
# states that will not be used in a subsequent conditional probabiity judgment. 
# The amount of the bias is determined condition.true.weight. States that will be
# used in a subsequent conditional probabiity judgment condition.true.weight times more
# than the states won't be used.
# (See Appendix F of Davis & Rehder, 2020.)
conditional.proposal.distribution = function (joint, condition, neighbors, condition.true.weight=10) {
  ps1 = function (nis) {
    ps = rep (0, nrow (joint))
    ps[nis] = sapply (nis, function (ni) if (joint[ni,]$p == 0) 1 else condition.true.weight)
    ps = ps / sum (ps) # Normalize
  }
  if (!is.null (condition))
    joint = joint.conditionalized (joint, condition)
  laply (neighbors, ps1)
}
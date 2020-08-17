#---------------------------------
# parameter variability model
#   what variability would we expect if people used bayes nets but had uncertainty over parameters
#---------------------------------
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source ("../utilities/utilities.R")
source ("../utilities/cgm.R")

# input: 
#   causal strength matrix (ms) and base rates for each cause (bs)
#   sds for parameters
#   number of samples to run
# output: distribution of joint distributions
param_variability <- function(ms,
                              bs,
                              ms_sd = .1,
                              bs_sd = .1,
                              nSamples = 100) {
    # creating storage for joints
    joint_store <- joint.skeleton(vs=c('x','y','z'))
    joint_store <- cbind(joint_store, 
                         data.frame(
                             matrix(ncol=nSamples,
                                    nrow=2**length(bs))))
    
    # randomly generating parameters, for now all the same ms/bs values
    stopifnot('can only handle one causal strength value right now' =
                  (ms[ms!=0] %>% unique %>% length) == 1)
    
    # sampling causal strengths
    ms_dist <- rnorm(n    = nSamples, 
                     mean = ms[ms!=0] %>% unique, # picking out the causal strength
                     sd   = ms_sd)
    ms_dist <- ms_dist %>% truncater(lower=0, upper=1)
    # sampling base rates
    bs_dist <- rnorm(n    = nSamples, 
                     mean = bs[bs!=0] %>% unique, # picking out base rates
                     sd   = bs_sd)
    bs_dist <- ms_dist %>% truncater(lower=0, upper=1)
    
    ms_indices <- ms != 0
    for (e in 1:nSamples) {
        # applying new causal strengths
        tmp_ms    <- ms
        tmp_ms[ms_indices] <- ms_dist[e]
        # applying new base rates
        tmp_bs <- rep(bs_dist[e], length(bs))
        
        # getting temporary joint
        tmp_joint <- joint.cgm.generic (ms=tmp_ms, 
                                        bs=tmp_bs)
        joint_store[,e+length(bs)] <- tmp_joint$p
    }
    
    return(joint_store)
}

#--------------------
# helper truncation function
truncater <- function(vec, lower, upper) {
    vec[vec < lower] <- lower
    vec[vec > upper] <- upper
    return(vec)
}

# Create joint distribution for chain network: X-> Y -> Z
#   Causal strengths (or "powers") = .75
#   Base rate of X = .50
#   Background causes of Y and Z = 1/3
var.names = c ('x','y','z')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['x','y'] = ms['y','z'] = .75
bs = c (.5, 1/3, 1/3)

a <- param_variability(ms=ms, 
                       bs=bs, 
                       ms_sd=.01,
                       bs_sd=.01,
                       nSamples=6)

a[1,4:ncol(a)] %>% as.numeric %>% hist



















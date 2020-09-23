#---------------------------------
# parameter variability model
#   what variability would we expect if people used bayes nets but had uncertainty over parameters
#---------------------------------
library(tidyverse)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source ("../utilities/utilities.R")
source ("../utilities/joint.utilities.R")
source ("../utilities/cgm.R")

# input: 
#   causal strength matrix (ms) and base rates for each cause (bs)
#   concentration for distribution around parameters
#   number of samples to run
# output: 
#   meanjoint: mean of sampled probability of each state
#   chainjoints: distribution of joint distributions
#   chainlens: chainlengths (not relevant for this model)
param_variability <- function(ms,
                              bs,
                              ms_conc = 1,
                              bs_conc = 1,
                              nSamples = 100) {
    # creating storage for joints
    joint_store <- joint.skeleton(vs=c('X1','Y','X2'))
    joint_store <- cbind(joint_store, 
                         data.frame(
                             matrix(ncol=nSamples,
                                    nrow=2**length(bs))))
    normjoint <- joint.cgm.generic3(ms, bs)
    normresps <- normjoint$normresps
    normjoint <- normjoint$joint
    

    # sampling causal strengths --------------------------------
    # array of ms (third dimension is each new sampled ms)
    ms_dist <- array(NaN, dim=c(dim(ms)[1], dim(ms)[2], nSamples))
    ms_indices <- which(ms != 0, arr.ind=TRUE)
    for (ms_idx in 1:nrow(ms_indices)) {
        ms_row <- ms_indices[ms_idx,'row']
        ms_col <- ms_indices[ms_idx,'col']
        # beta noise around parameters
        shape1 <- (ms[ms_row, ms_col] * (ms_conc - 2)) + 1
        shape2 <- ms_conc - shape1
        ms_dist[ms_row, ms_col, ] <- 
            rbeta(nSamples, shape1, shape2)
    }
    ms_dist[is.nan(ms_dist)] <- 0
    
    # sampling base rates --------------------------------
    bs_dist <- array(NaN, dim=c(length(bs), nSamples))
    for (bs_idx in 1:length(bs)) {
        # beta noise around parameters
        shape1 <- (bs[bs_idx] * (bs_conc - 2)) + 1
        shape2 <- bs_conc - shape1
        bs_dist[bs_idx, ] <- 
            rbeta(nSamples, shape1, shape2)
    }
    
    for (e in 1:nSamples) {
        # applying new causal strengths
        tmp_ms    <- ms_dist[,,e]
        rownames(tmp_ms) <- colnames(tmp_ms) <- c('x','y','z')
        # applying new base rates
        tmp_bs <- bs_dist[,e]
        
        # getting temporary joint
        tmp_joint <- joint.cgm.generic (ms=tmp_ms, 
                                        bs=tmp_bs)
        joint_store[,e+length(bs)] <- tmp_joint$p
    }
    
    #output needs to be list with meanjoint chainjoints and chainlens
    meanjoint <- joint_store[,1:3]
    meanjoint$p <- apply(joint_store[,4:ncol(joint_store)], 1, mean) #mean joint distribution
    meanjoint <- meanjoint[nrow(meanjoint):1,] #reverse order so it follows IK order of states convention 111-000
    
    chainjoints <- data.frame(matrix(ncol=2**length(bs), nrow=nSamples))
    chainjoints <- t(joint_store[,4:ncol(joint_store)])
    chainjoints <- chainjoints[,ncol(chainjoints):1] #reverse order so it follows IK convention
    
    chainlens <- rep(1, nSamples) #vector of fakechainlengths, required by fnc genpreddistr, can only use this fnc with betavar=0 for this model.
    
    res <- list(meanjoint=meanjoint, chainjoints=chainjoints, 
                chainlens=chainlens, normjoint=normjoint, ms=ms, bs=bs, ms_conc=ms_conc, bs_conc=bs_conc, nSamps=nSamples, normresps=normresps,
                model="parameter variability model")
    
    return(res)
}


# Create joint distribution for chain network: X-> Y -> Z
#   Causal strengths (or "powers") = .75
#   Base rate of X = .50
#   Background causes of Y and Z = 1/3

# var.names = c ('x','y','z')
# ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
# ms['x','y'] = ms['y','z'] = .75
# bs = c (.75, .75, .75)
# 
# a <- param_variability(ms = ms, 
#                        bs = bs, 
#                        ms_conc = .01,
#                        bs_conc = .01,
#                        nSamples = 6)
# 
# a$chainjoints[1,] %>% as.numeric %>% hist



















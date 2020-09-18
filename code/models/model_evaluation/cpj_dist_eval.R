#------------------------------------
# evaluating the distributions of conditional probability judgments
#------------------------------------
library(tidyverse)
library(LaplacesDemon)   # gives kl-divergence through function KLD
library(multimode)       # for evaluating the modality of a distribution
library(mixtools)        # for identifying the parameters of mixed gaussians

# KL-Divergence -----------------
# basic KL-Divergence that can be used to evaluate how similar two networks are
px <- dnorm(runif(100),0,1)
py <- dnorm(runif(100),0.1,0.9)
KLD(px,py)

# qualitative summary -------------------
qualitative_eval <- function(judgments, bw=.1, threshold=.05) {
    # number of modes
    n_modes <- nmodes(judgments, bw=bw, full.result = FALSE)
    # location of modes
    loc_modes <- normalmixEM(pz, k = n_modes, fast=TRUE)
    # had to get a little clever but threshold is distance from 0 or 1
    # e.g. threshold of .05 finds all <.05 and >.95
    extremes <- mean(abs(judgments - .5) > (.5 - threshold))
    # variance
    judge_var <- var(judgments)
    
    return(list(n_modes = n_modes,
                loc_modes = loc_modes,
                extremes = extremes,
                judge_var = judge_var))
}









#------------------------------------
# evaluating the distributions of conditional probability judgments
#------------------------------------
library(LaplacesDemon)   # gives kl-divergence through function KLD
library(multimode)       # for evaluating the modality of a distribution
library(mixtools)        # for identifying the parameters of mixed gaussians

# basic KL-Divergence that can be used to evaluate how similar two networks are
px <- dnorm(runif(100),0,1)
py <- dnorm(runif(100),0.1,0.9)
KLD(px,py)

# generating some mixture data
p1 <- rnorm(200, mean=0, sd=1)
p2 <- rnorm(200, mean=4, sd=1)
p3 <- rnorm(200, mean=8, sd=1)
pz <- c(p1,p2,p3)

# testing modality
nmodes(c(p1,p2,p3), bw=1, full.result=TRUE) %>% print

# gaussian mixture modeling to identify modes
normalmixEM(pz, fast=TRUE)

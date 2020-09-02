#------------------------------------
# evaluating the distributions of conditional probability judgments
#------------------------------------
library(LaplacesDemon)   # gives kl-divergence through function KLD
library(multimode)

px <- dnorm(runif(100),0,1)
py <- dnorm(runif(100),0.1,0.9)
KLD(px,py)

print(hist(px))
nmodes(px, .01) %>% print

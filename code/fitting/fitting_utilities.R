#----------------------------------------------------------
# Model fitting utilities, PDA method 
#----------------------------------------------------------


#----------------------------------------------------------
# Smoothing kernel function
#----------------------------------------------------------  

#inputs: predicted response distribution for a single inference(type) based on samples (so not for betavar analytic model)
#output: smoothed distribution over each percentage point

smoothkernel <- function(pred.pdf){
  #simulated PDF, ie density of predictions at each percentage point, according to Turner&Sederberg PDA paper
  dens <- density(pred.pdf, from=0, to=1, n=101, bw="nrd0", kernel="epanechnikov") 
  return(dens$y)
}



#----------------------------------------------------------
# PDA function
#----------------------------------------------------------  
# inputs: SPDF (synthetic PDF) (output from smoothkernel fnc above or guessfnc) and response data, for single inference(type)
# output: SLL (approximated summed log likelihood.)

PDA.SLL <- function(spdf, responses){
  SLL <-  rep(NA,length(responses))
  SLL <- pmax(spdf[(responses*100)+1], 1e-10) # take the density at the points (x-value of dens) of the responses, if value is zero take 1e-10 to prevent log(0)
  SLL <- sum(log(SLL)) #take the logarithm and them sum them to get SLL
  return(SLL)
}


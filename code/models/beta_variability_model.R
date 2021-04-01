#----------------------------------------------------------
# Betavar model
#----------------------------------------------------------

# Predicts responses are distributed according to beta distribution with mode centered on correct/CBN answer
# free parameter: concentration, ie how much variability.

# Beta distr has 2 shape parameters. The mode = (shape1-1)/(shape1+shape2-2)  see https://en.wikipedia.org/wiki/Beta_distribution#Mode
# so when optimizing parameters shape1 and shape 2 to fit to data, we need to constrain as follows:
# shape1 = normresp*(shape1 + shape2 -2) + 1, form which we get: 
# shape1 = (normresp*shape2 - 2*normresp + 1)/(1 -  normresp) #the constraint on fitting shape1 and shape 2
# shape1 + shape2 is the socalled 'concentration' parameter, this is what we will have to be optimized under above constraint when fitting.


inferencenames <- c( 'X1|X2==0', #used here by betavar function to name output columns, should put this in utilities script
                     'X1|X2==1', 
                     'X1|Y==0',
                     'X1|Y==0 & X2==0',
                     'X1|Y==0 & X2==1',
                     'X1|Y==1',
                     'X1|Y==1 & X2==0',
                     'X1|Y==1 & X2==1',
                     'Y|X2==0',
                     'Y|X2==1', 
                     'Y|X1==0',
                     'Y|X1==0 & X2==0',
                     'Y|X1==0 & X2==1',
                     'Y|X1==1',
                     'Y|X1==1 & X2==0',
                     'Y|X1==1 & X2==1',
                     'X2|Y==0',
                     'X2|Y==1',
                     'X2|X1==0',
                     'X2|Y==0 & X1==0',
                     'X2|Y==1 & X1==0',
                     'X2|X1==1', 
                     'X2|Y==0 & X1==1',
                     'X2|Y==1 & X1==1',
                     'X1',
                     'Y',
                     'X2'
)


#----------------------------------------------------------
# Betavar model functions (04/2021)
#----------------------------------------------------------  
#inputs: causal strengths, baserates, concentration (free parameter Betavar model)


betavar.pdf.samps <- function(ms, bs, concentration, nSamps=10000){
  #output: predicted distributions for each inference
  joint = joint.cgm.generic3(ms, bs) #this computes multiple things while only normresp is needed here
  respdistr <- data.frame(matrix(NA, nSamps, 27))
  colnames(respdistr) <- inferencenames
  
  for (i in 1:27){
    shape1 <- as.numeric(joint$normresps[i]*(concentration-2) + 1) #solve eq for mode for shape1, plugging in norm resp and shape1=concentration-shape2
    shape2 <- concentration - shape1
    respdistr[,i] <- rbeta(nSamps, shape1, shape2)
  }
  res <- respdistr
  return(res)
}


betavar.pdf.analytic <- function(ms, bs, concentration){
  # output: density values at each percentage point, from 0 to 100, ie at 101 points
  # (density at x% = index-1, e.g. density[1] is density at response 0%, density[33] is density at response 32%)
  joint = joint.cgm.generic3(ms, bs) #this is a waste, only need normresp here
  respdistr <- data.frame(matrix(NA, 101, 27))
  colnames(respdistr) <- inferencenames
  xvals <- seq(0, 1, .01)
  
  for (i in 1:27){
    shape1 <- as.numeric(joint$normresps[i]*(concentration-2) + 1) #solve eq for mode for shape1, plugging in norm resp and shape1=concentration-shape2
    shape2 <- concentration - shape1
    respdistr[,i] <- dbeta(xvals, shape1, shape2)
  }
  res <- respdistr
  return(res) 
}


#----------------------------------------------------------
# Original Betavar model simulation function (used for exploratory simulations in 2020)
#----------------------------------------------------------  

betavariability.sim <- function(joint, concentration, nSamps){
  # inputs: joint distribution (list with joint, ms, bs, normresps), concentration parameter, number of samples
  # output: distribution of inferences, and a lot of other information
  normresp <- joint$normresps
  ms <- joint$ms
  bs <- joint$bs
  joint <- joint$joint
  
  respdistr <- data.frame(matrix(NA, nSamps, 27))
  colnames(respdistr) <- inferencenames
  
  for (i in 1:27){
    shape1 <- as.numeric(normresp[i]*(concentration-2) + 1) #solve eq for mode for shape1, plugging in normresp and shape1=concentration-shape2
    shape2 <- concentration - shape1
    respdistr[,i] <- rbeta(nSamps, shape1, shape2)
  }
  
  res <- list(model="Beta distr normative inference", respdistr=respdistr, normjoint=joint, concentration=concentration, nSamps=nSamps, normresps=normresp, ms=ms, bs=bs)
  return(res)
}



#----------------------------------------------------------
# Helper fnc
#----------------------------------------------------------

# not used now
betashape2to1 <- function (shape2, normresp){
  shape1 <- (normresp*shape2 - 2*normresp + 1)/(1 -  normresp)
  return(shape1)
}

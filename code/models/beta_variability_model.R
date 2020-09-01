# Beta variability model, functions as a sort of null hypothesis
# Predicts responses are distributed according to beta distribution with mode centered on correct/normative CBN answer
# free parameter: concentration, ie how much variability.


# Beta distr has 2 shape parameters. The mode = (shape1-1)/(shape1+shape2-2)  see https://en.wikipedia.org/wiki/Beta_distribution#Mode
# so when optimizing parameters shape1 and shape 2 to fit to data, we need to constrain as follows:
# shape1 = normresp*(shape1 + shape2 -2) + 1, form which we get: 
# shape1 = (normresp*shape2 - 2*normresp + 1)/(1 -  normresp) #the constraint on fitting shape1 and shape 2
# shape1 + shape2 is the socalled 'concentration' parameter, this is what we will have to be optimized under above constraint when fitting..

inferencenames <- c( 'X1|X2==0',
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
# Simulate beta variability model
#----------------------------------------------------------


betavariability.sim <- function(joint, concentration, nSamps){
  # inputs: joint distribution, concentration parameter, number of samples
  # output: distribution of inferences
  normresp <- list()
  normresp[[1]]<-condproballBay(joint, 0, 10) #nested list for compatibility with fnc preddistrBay
  normresp <- preddistrBay(1, normresp) #DF with each normative response
  
  res <- data.frame(matrix(NA, nSamps, 27))
  colnames(res) <- inferencenames
  
  for (i in 1:27){
    shape1 <- as.numeric(normresp[i]*(concentration-2) + 1) #solve eq for mode for shape1, plugging in norm resp and shape1=concentration-shape2
    shape2 <- concentration - shape1
    res[,i] <- rbeta(nSamps, shape1, shape2)
  }
  return(res)
}



#----------------------------------------------------------
# Fit beta variability model: objective function
#----------------------------------------------------------

# !!!!! input needs to be checked. optimize for shape parameter1 or for concentration parameter? !!!!!!!!!

#Objective function to optimize for fitting
Beta_obj1 <- function(par, response, normresp) {
  # inputs: shape parameter 1(to be optimized), vector of responses, vector of correct answers
  # output: Negative summed log-likelihood (NegSLL) for all inferences, fitted concentration parameter
  
  shape2 <- par #shape parameter2
  SLLinfs <- rep(NA, length(response))
  
  #vectorize loop using sapply?
  for (i in length(response)){ #loop through each inference #could also loop through all unique values of normresp. Or, do sapply
    resp <- round(response[i], digits=2)
    shape1 <- (normresp[i]*shape2 - 2*normresp[i] + 1)/(1 -  normresp[i]) #constrain betadistr shape1 par so that mode = normresp
    SLLinfs[i] <-  log(max(dbeta(resp,shape1,shape2), 1e-10)) #save density at value of response as syntethic SLL??
  }
  negSLL <- -sum(SLLinfs)
  return(negSLL)
}

#----------------------------------------------------------
# Helper fnc
#----------------------------------------------------------

# not used now
betashape2to1 <- function (shape2, normresp){
  shape1 <- (normresp*shape2 - 2*normresp + 1)/(1 -  normresp)
  return(shape1)
}

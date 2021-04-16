#----------------------------------------------------------
# Functions to calculate conditional and marginal probabilities of events based on joint distributions
# includes option to incorporate symmetric beta prior, if no prior wanted set betavar=0
#----------------------------------------------------------



condprobBay <- function (x, betavar, len, event = NULL, given = NULL) {
  # calculates conditional or marginal probability based on combination of prior and relative freq in samples
  # X is dataframe with column names as variables X1, Y, X2, and p  (probability)
  # betavar refers to prior, prior is beta(betavar,betavar)
  # len = amount of samples / chain length.
  # event is the variable and value of which the probability is calculated, given is the condition statement
  # event and given need to evaluate to logical expreesions (e.g. X1==1 or (X1=x=1&X2==0))
  # returns NaN if event and conditioning state is never observed
  if (!is.numeric(betavar)||!is.numeric(len)) {
    stop("betavar and or len is missing")
  }
  if (missing(event)) {
    r <- TRUE
  }
  else {
    e <- substitute(event)
    r <- eval(e, x, parent.frame())
    r <- r & !is.na(r)
    if (!isTRUE(all.equal(sum(x$p), 1))) 
      warning("'space' does not have probability 1.")
  }
  A <- x[r, ] #part of x where event holds
  if (missing(given)) {
    p <- sum(A$p)
  }
  else {
    f <- substitute(given)
    g <- eval(f, x, enclos = parent.frame())
    if (!is.logical(g)) {
      B <- given
    }
    else {
      g <- g & !is.na(g)
      B <- x[g, ] #part of x where given holds
    }
    p <- ((sum(intersect(A, B)$p)*len)+betavar)/(sum((B$p)*len) + 2*betavar) #the estimated probability, with the beta prior incorporate as a 'pseudo count'
  }
  return(p)
}

condproballBay <- function(dagdata, betavar, len){
  #only for discrete data, values need to be 0 and 1
  
  cndprbs <- list(
    condprobBay(dagdata, betavar, len, X1==1),
    condprobBay(dagdata, betavar, len, X1==1, (Y==0)),
    condprobBay(dagdata, betavar, len, X1==1, (Y==1)),
    condprobBay(dagdata, betavar, len, X1==1, (X2==0)),
    condprobBay(dagdata, betavar, len, X1==1, (X2==1)),
    condprobBay(dagdata, betavar, len, X1==1, (Y==0 & X2==0)),
    condprobBay(dagdata, betavar, len, X1==1, (Y==1 & X2==0)),
    condprobBay(dagdata, betavar, len, X1==1, (Y==0 & X2==1)),
    condprobBay(dagdata, betavar, len, X1==1, (Y==1 & X2==1)),
    condprobBay(dagdata, betavar, len, Y==1),
    condprobBay(dagdata, betavar, len, Y==1, (X1==0)),
    condprobBay(dagdata, betavar, len, Y==1, (X1==1)),
    condprobBay(dagdata, betavar, len, Y==1, (X2==0)),
    condprobBay(dagdata, betavar, len, Y==1, (X2==1)),
    condprobBay(dagdata, betavar, len, Y==1, (X1==0 & X2==0)),
    condprobBay(dagdata, betavar, len, Y==1, (X1==1 & X2==0)),
    condprobBay(dagdata, betavar, len, Y==1, (X1==0 & X2==1)),
    condprobBay(dagdata, betavar, len, Y==1, (X1==1 & X2==1)),
    condprobBay(dagdata, betavar, len, X2==1),
    condprobBay(dagdata, betavar, len, X2==1, (Y==0)),
    condprobBay(dagdata, betavar, len, X2==1, (Y==1)),
    condprobBay(dagdata, betavar, len, X2==1, (X1==0)),
    condprobBay(dagdata, betavar, len, X2==1, (X1==1)),
    condprobBay(dagdata, betavar, len, X2==1, (Y==0 & X1==0)),
    condprobBay(dagdata, betavar, len, X2==1, (Y==1 & X1==0)),
    condprobBay(dagdata, betavar, len, X2==1, (Y==0 & X1==1)),
    condprobBay(dagdata, betavar, len, X2==1, (Y==1 & X1==1))
  )
  
  
  names(cndprbs) <- c(
    'X1',
    'X1|Y==0',
    'X1|Y==1',
    'X1|X2==0',
    'X1|X2==1', 
    'X1|Y==0 & X2==0',
    'X1|Y==1 & X2==0',
    'X1|Y==0 & X2==1',
    'X1|Y==1 & X2==1',
    'Y',
    'Y|X1==0',
    'Y|X1==1',
    'Y|X2==0',
    'Y|X2==1', 
    'Y|X1==0 & X2==0',
    'Y|X1==1 & X2==0',
    'Y|X1==0 & X2==1',
    'Y|X1==1 & X2==1',
    'X2',
    'X2|Y==0',
    'X2|Y==1',
    'X2|X1==0',
    'X2|X1==1', 
    'X2|Y==0 & X1==0',
    'X2|Y==1 & X1==0',
    'X2|Y==0 & X1==1',
    'X2|Y==1 & X1==1'
  )
  return(cndprbs)
}


condprobMSBay <- function(betavar, nchains, varsamples, varjoint){ 
  # Function to do conditional or marginal prob query on samples from runchainsMS
  chainprobs <- lapply(1:nchains, function(i) {
    condproballBay(varjoint[[i]], betavar, nrow(varsamples[[i]])) 
  })
  return(chainprobs)
}

condprobMSBay2 <- function(betavar, nchains, chainlens, varjoint){ 
  # Function to do conditional or marginal prob query on samples from runchainsMS
  chainprobs <- lapply(1:nchains, function(i) {
    condproballBay(varjoint[[i]], betavar, chainlens[i]) 
  })
  return(chainprobs)
}

guessMS <- function(x, nchains){
  # function to change NaNs to .5 in the prob inferences, according to MS model this happens when states are not visited
  # input should be output of condprobMS function.
  for (i in 1:nchains){ #loop through chains #LvM09062020: Here, I did not replace the loop with lapply, because lapply removes the names of the list
    for (ii in 1:27){   #loop through each prob query
      if (is.nan(x[[i]][[ii]])){
        x[[i]][ii] <- .5
      }
    }
  }
  return(x)
}


#function to get DF with distribution of all inferences
preddistrBay <- function(nchains, chainMSprobsBay, infnames=0){
  # puts marginal and conditional probability queries into tidy dataframe, can do multiple chains (and thus joints) or single joint
  # argument chainMSprobsBay is output from condproballBay (either for chain joints or normative joint)
  preddistr <- data.frame(matrix(ncol = 27, nrow = nchains))
  
  if (infnames==0) {
    colnames(preddistr) <- c( 'X1|X2==0',
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
  } else {
    colnames(preddistr) <- paste0("infnr", as.character(seq(1:27)))
  }
  
  for (i in 1:nchains){ 
    preddistr[i,1] <- chainMSprobsBay[[i]]$`X1|X2==0` #infnr 1
    preddistr[i,2] <- chainMSprobsBay[[i]]$`X1|X2==1` #infnr 2
    preddistr[i,3] <- chainMSprobsBay[[i]]$`X1|Y==0` #infnr 3
    preddistr[i,4] <- chainMSprobsBay[[i]]$`X1|Y==0 & X2==0` #infnr 4
    preddistr[i,5] <- chainMSprobsBay[[i]]$`X1|Y==0 & X2==1` #infnr 5
    preddistr[i,6] <- chainMSprobsBay[[i]]$`X1|Y==1` #infnr 6
    preddistr[i,7] <- chainMSprobsBay[[i]]$`X1|Y==1 & X2==0` #infnr 7
    preddistr[i,8] <- chainMSprobsBay[[i]]$`X1|Y==1 & X2==1` #infnr 8
    preddistr[i,9] <- chainMSprobsBay[[i]]$`Y|X2==0` #infnr 9
    preddistr[i,10] <- chainMSprobsBay[[i]]$`Y|X2==1` #infnr 10
    preddistr[i,11] <- chainMSprobsBay[[i]]$`Y|X1==0` #infnr 11
    preddistr[i,12] <- chainMSprobsBay[[i]]$`Y|X1==0 & X2==0` #infnr 12
    preddistr[i,13] <- chainMSprobsBay[[i]]$`Y|X1==0 & X2==1` #infnr 13
    preddistr[i,14] <- chainMSprobsBay[[i]]$`Y|X1==1` #infnr 14
    preddistr[i,15] <- chainMSprobsBay[[i]]$`Y|X1==1 & X2==0` #infnr 15
    preddistr[i,16] <- chainMSprobsBay[[i]]$`Y|X1==1 & X2==1` #infnr 16
    preddistr[i,17] <- chainMSprobsBay[[i]]$`X2|Y==0` #infnr 17
    preddistr[i,18] <- chainMSprobsBay[[i]]$`X2|Y==1` #infnr 18
    preddistr[i,19] <- chainMSprobsBay[[i]]$`X2|X1==0` #infnr 19
    preddistr[i,20] <- chainMSprobsBay[[i]]$`X2|Y==0 & X1==0` #infnr 20
    preddistr[i,21] <- chainMSprobsBay[[i]]$`X2|Y==1 & X1==0` #infnr 21
    preddistr[i,22] <- chainMSprobsBay[[i]]$`X2|X1==1` #infnr 22
    preddistr[i,23] <- chainMSprobsBay[[i]]$`X2|Y==0 & X1==1` #infnr 23
    preddistr[i,24] <- chainMSprobsBay[[i]]$`X2|Y==1 & X1==1` #infnr 24
    preddistr[i,25] <- chainMSprobsBay[[i]]$`X1` #infnr 25
    preddistr[i,26] <- chainMSprobsBay[[i]]$`Y` #infnr 26
    preddistr[i,27] <- chainMSprobsBay[[i]]$`X2` #infnr 27
  }
  return(preddistr)
}

normativeresps <- function(ms, bs){
  #outputs all conditional and marginal probabilities (ie normative responses) based on causal strengths and baserates
  j = joint.cgm.generic2(ms, bs)
  normresps <- list()
  normresps[[1]]<-condproballBay(j, 0, 10)
  normresps <- preddistrBay(1, normresps, infnames=0)
  return (normresps)
}

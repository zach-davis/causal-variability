this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library ("plyr")
library ("reshape")
library ('gtools')

if(1) {
    source ("helper_functions/lib.r")
    source ("helper_functions/print.lib.r")
    source ("helper_functions/hw.r")
    source ("helper_functions/plib.r")
    source ("helper_functions/cgm.r")
    source ("helper_functions/chain.models.r")
    source ("helper_functions/sample.cgm.r")
    source ("helper_functions/learn.cgm.r")
}

#######################################################################################################
#####                                                                                             #####
#####   Returns matrix of state transition probabilities.                                         #####
#####                                                                                             #####
#######################################################################################################
jl.transitions.in.parallel = FALSE
transition.ps = function (j, neighbors=NULL, method='mutate', q.ps=NULL, gx=NULL) {
    
    # Transition probabilities for the standard "neighbors" proposal distribution. 
    mutate.transitions = function (xi) {
        nis = neighbors [[xi]]
        # Each neighbor's Hastings ratio, times the probability that we'll sample that neighbor.
        ps.of.xi = rep (0, nrow (j))
        ps.of.xi[nis] = sapply (nis, function (ni) { 
            min (c (j$p[ni] * q.ps[ni, xi] / (j$p[xi] * q.ps[xi, ni]), 1)) * q.ps[xi, ni]
        })
        stopifnot (!any (is.na (ps.of.xi)))
        # Rest of probability mass goes to current state.
        ps.of.xi[xi] = 1 - sum (ps.of.xi) 
        ps.of.xi
    }
    
    # Transition probabilities based on Gibbs sampling
    gibbs.transitions = function (xi, depth=length (joint.vs (j))) {
        # Recurse down the tree of variables/values in the order specified by vis.
        ps1 = function (vis) {
            recurse = function (xi, vii, ps) {
                stopifnot (rownames (j)[xi] == names (neighbors)[xi])
                vi = vis[vii] # Variable to apply Gibbs update.
                ni = neighbors[[xi]][vi] # New state if Gibbs variable changes.
                new.val.p = j$p[ni] / (j$p[ni] + j$p[xi]) # Probability that it will change.
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
            ps.of.xi = recurse (xi, 1, rep (0, nrow (j)))
            ps.of.xi / sum (ps.of.xi) 
        }
        # All possible variable orders (but only go to depth specified by depth).
        vis = permutations (length (joint.vs (j)), depth) 
        # Transition probabilities for each order.
        ps.of.xi = aaply (vis, 1, ps1)
        # Average over orders.
        apply (ps.of.xi, 2, mean)
    }
    
    gibbs1.transitions = function (xi) 
        gibbs.transitions (xi, 1)
    
    nr = nrow (j)
    # Each row xi has transition probabilities for state j[xi,] 
    j = joint.normalized (j)
    if (is.null (neighbors))
        # By default, all nodes (states) can be proposed.
        neighbors = setNames (lapply (1:nr, function (xi) seq (1, nr)[-xi]), rownames (j))
    if (is.null (q.ps)) 
        # By default, proposal distribution is uniform over neighbors.
        q.ps = sapply (neighbors, function (nis) {
            sapply (1:nr, function (ni) if (ni %in% nis) 1 / length (nis) else 0)
        })
    else
        stopifnot (method == 'mutate') # q.ps only applies for mutation sampler.
    transition.funs = list (mutate=mutate.transitions, gibbs=gibbs.transitions, gibbs1=gibbs1.transitions)
    stopifnot (!is.null (method) & method %in% names (transition.funs))
    transition.ps = laply (1:nr, transition.funs[[method]], .parallel=jl.transitions.in.parallel)
    if (!is.null (gx)) {
        jc = joint.conditionalized (j, gx)
        transition.ps = sapply (1:nr, function (r) if (jc[r,]$p == 0) rep (0, nr) else transition.ps[,r])
        transition.ps = t (apply (transition.ps, 1, function (ps) ps / sum (ps)))
    }
    transition.ps
}

#######################################################################################################
#####                                                                                             #####
#####   Master function, takes in 5 parameters:                                                   #####
#####     - dist:       the joint distribution to sample over                                     #####
#####     - chainLen:      length of 1 chain                                                      #####
#####     - nChains:      number of runs (number of chains)                                       #####
#####     - bias:       start at outside nodes, or start randomly anywhere in the space?          #####
#####     - propStruct: use proposal structure, or sample from any point at random?               #####
#####                                                                                             #####
#####   jlMaster will start at a some point in the joint, and jump to another stimulus type based #####
#####   on the MH rule. The proposal structure looks as follows, where every line drawn between   #####
#####   nodes means that they are 1 switch away from each other (and thus accessible by our       #####
#####   proposal structure)                                                                       #####
#####                                                                                             #####
#####               001 - 011                                                                     #####
#####             /     X     \                                                                   #####
#####         000 - 010   101 - 111                                                               #####
#####             \     X     /                                                                   #####
#####               100 - 110                                                                     #####
#####                                                                                             #####
#####   Having a small "chainLen" parameter implies a large bias, since there is less time to escape #####
#####   from a local region. "nChains" does not change the bias, just reduces the variance.         #####
#####                                                                                             #####
#######################################################################################################
jl.sample.in.parallel = FALSE

jlMaster <- function (j, chainLen, nChains, bias, propStruct, neighbors=NULL, method='mutate', q.ps=NULL, gx=NULL) {
    sample.one.cpu = function (nChains) {
        # precomputing random numbers to access next rows
        mcmcWinner <- rep(1e-10,no.states) # to store data, 1e-10 is to avoid na's, shouldn't matter
        # number of times to simulate the experiment
        for (simNum in 1:nChains) {
            # if biased is requested, then only start at 000 or 111, else sample randomly
            if (bias)
                current <- sample(c(1,no.states),1)
            else
                current <- sample(1:no.states,1,prob=rep(1/no.states,no.states))
            # Run one chain
            for (t in 1:chainLen) {
                mcmcWinner[current] <- mcmcWinner[current] + 1  # counts the starting point
                current <- sample((1:nrow(j)), 1, prob=transition.ps[current,])
            }
        }
        return (mcmcWinner / sum (mcmcWinner))
    }
    stopifnot (nChains > 0 & chainLen > 0)
    if (propStruct)
        stopifnot (!is.null (neighbors))
    else
        neighbors = NULL
    no.states <- nrow(j)
    transition.ps = transition.ps (j, neighbors, method, q.ps, gx)
    # allocate chains to cpus
    Chains.per.cpu = jobs.per.cpu (nChains)
    mcmcWinner <- ldply (Chains.per.cpu, sample.one.cpu, .parallel=jl.sample.in.parallel, .inform=F)
    # average over cpus
    mcmcWinner <- apply (mcmcWinner, 2, mean)
    stopifnot (equal.within.tol (sum (mcmcWinner), 1))
    mcmcWinner = mcmcWinner / sum (mcmcWinner)
    
    #displaying
    a <- rbind(t(j[,"p"]),mcmcWinner) # concatenating the output of the normative values and simulation values
    a <- rbind(a, a[2,]-a[1,])           # computing the difference between the two
    rownames(a) <- c("true","mcmc","diff")
    colnames(a) <- rownames (j)
    return(a)
}

jlMaster.original <- function (dist, chainLen, nChains, bias, propStruct, neighbors=NULL, q.ps=NULL) {
    
    sample.one.cpu = function (nChains) {
        # precomputing random numbers to access next rows
        ndimensions <- ncol(dist)-1
        if (propStruct) {
            rowVec <- sample(ndimensions, chainLen*nChains, replace=TRUE)
            dim(rowVec) <- c (nChains, chainLen)
            stopifnot (!is.null (neighbors))
        }
        mcmcWinner <- rep(1e-10,nrow(dist)) # to store data, 1e-10 is to avoid na's, shouldn't matter
        hastVec <- runif(chainLen*nChains) # precomputing random numbers for the hastings ratio
        dim(hastVec) <- c(nChains, chainLen)
        # number of times to simulate the experiment
        for (simNum in 1:nChains) {
            # if biased is requested, then only start at 000 or 111, else sample randomly
            if (bias)
                current <- sample(c(1,nrow(dist)),1)
            else
                current <- sample(1:nrow(dist),1,prob=rep(1/nrow(dist),nrow(dist)))
            # Run one chain
            for (t in 1:chainLen) {
                mcmcWinner[current] <- mcmcWinner[current] + 1  # counts the starting point
                # if there is a proposal structure, choose a proposal node that differs by only
                # 1 value. Otherwise, choose randomly
                if (propStruct)
                    propRow <- neighbors[[current]][rowVec[simNum, t]]
                else
                    propRow <- sample((1:nrow(dist))[(1:nrow(dist))!=current], 1)
                # Hastings ratio which is either 1 (if the proposal is more likely than the current)
                # or the ratio between the proposal and current (if proposal is less likely)
                hastRatio <- min(c(dist[propRow,"p"]/dist[current,"p"], 1))
                if (hastRatio > hastVec[simNum, t]) 
                    current <- propRow
            }
        }
        return (mcmcWinner / sum (mcmcWinner))
    }
    transition.ps = transition.ps (j, neighbors, method, q.ps)
    # allocate chains to cpus
    stopifnot (nChains > 0 & chainLen > 0)
    Chains.per.cpu = jobs.per.cpu (nChains)
    mcmcWinner <- ldply (Chains.per.cpu, sample.one.cpu, .parallel = jl.sample.in.parallel, .inform = F)
    # average over cpus
    mcmcWinner <- apply (mcmcWinner, 2, mean)
    stopifnot (equal.within.tol (sum (mcmcWinner), 1))
    mcmcWinner = mcmcWinner / sum (mcmcWinner)
    
    #displaying
    a <- rbind(t(dist[,"p"]),mcmcWinner) # concatenating the output of the normative values and simulation values
    a <- rbind(a, a[2,]-a[1,])           # computing the difference between the two
    rownames(a) <- c("true","mcmc","diff")
    colnames(a) <- for (i in 1:(nrow(dist))) (paste(dist[i,1:ncol(dist)-1],sep="",collapse=""))
    return(a)
}


numRuns <- 100       # number of times to run the simulation
chainLen <- 100      # length of chain
nChains <- 10        # number of chains
var1 <- 'za'         # variables to test correlation
var2 <- 'zb'         # variables to test correlation
propStruct=TRUE      # should we use the proposal structure?
bias=FALSE           # should we start at the ends?

########################################################################################
###########          Gives a distribution of participants                    ###########
###########                                                                  ###########
########### - Introduces a new parameter "numRuns" to designate how many     ###########
###########   participants to simulate                                       ###########
########################################################################################
jlDist <- function(dist, chainLen, nChains, bias, propStruct, numRuns) {
    c <- data.frame()
    for (t in 1:numRuns) {
        b <- jlMaster(dist=dist, chainLen=chainLen, nChains=nChains, bias=bias, propStruct=propStruct)
        #c[t,] <- b[2,]
        c <- rbind(c, b[2,])
        if (t%%100==0) print(t)
    }
    
    d <- melt(c)
    plot(d[,1],jitter(d[,2]), ylim=c(0,1), xaxt='n')
    
    lines(1:length(dist$p),dist$p,col="red")
    
    title(main=paste0(if(bias==FALSE)"un","biased distribution with", if(propStruct==FALSE)"out"," proposal structure\n",
                      chainLen," chain length, ",nChains," chains"))
    
    axislabels <- rep(0,nrow(dist))
    for (i in 1:(nrow(dist))) axislabels[i] <- (paste(dist[i,1:ncol(dist)-1],sep="",collapse=""))
    axis(1,at=1:length(dist$p),labels= axislabels)
    return(d)
}
#dist <- rep(1/8,8)
#jlDist(dist,chainLen,nChains,bias,propStruct,numRuns)

########################################################################################
###########              Returns independence violations                     ###########
###########                                                                  ########### 
########### var1 and var2 denote the two variables to return the lods for    ###########
########################################################################################
bq_jlFit <- function(dist, chainLen, nChains, bias, propStruct, numRuns, var1, var2) {
    lodStore <- c()
    for (i in 1:numRuns) {
        a <- jlMaster(dist,chainLen,nChains,bias,propStruct)
        sampleOut <- dist
        sampleOut$p <- a[2,]
        logOdds <- lod(sampleOut, var1,var2)
        
        lodStore <- c(lodStore, logOdds)
        if (i%%100==0) print(i)
    }
    return(lodStore)
}


# Provides a histogram of Markov violations
if(0) {
    lodStore <- bq_jlFit(dist,chainLen,nChains,bias,propStruct,numRuns,var1,var2)
    lod_noInf <- lodStore[!lodStore>100]
    lod_noInf <- lod_noInf[!lod_noInf<(-100)]
    lodInfo <- c(mean(lod_noInf, na.rm=T), sd(lod_noInf,na.rm=T), sum(lodStore>100, na.rm=T), sum(lodStore<(-100), na.rm=T))
    names(lodInfo) <- c("mean","sd","Inf","-Inf")
    print(lodInfo)
    
    hist(lodStore,col="grey", breaks=20, xlab=paste0("log odds ratio for ",chainLen," chain length, ",nChains," chains\nmean = ", 
                                                     round(mean(lod_noInf,na.rm=T),digits=2),"; sd = ",round(sd(lod_noInf,na.rm=T),digits=2)), 
         main=paste0(if(bias==FALSE)"un","biased distribution with", if(propStruct==FALSE)"out"," proposal structure\n", var1,":",var2))
}

###########################################################################################
####### For each entry in joint j, compute its neighbors (entries that differ by     ######
####### one binary variable).                                                        ######
###########################################################################################
neighbors.of.j = function (j) {
    neigbors.of.xi = function (xi) {
        flip1 = function (vi) { 
            substr (xi.rowname, vi, vi) = if (substr (xi.rowname, vi, vi) == '1') '0' else '1'; 
            xi.rowname 
        }
        xi.rowname = rownames (j)[xi]
        # Get rownames of neighbors.
        ns.of.xi = sapply (1:length (joint.vs (j)), flip1)
        # Convert to row numbers.
        sapply (ns.of.xi, function (n) seq (1, nrow (j))[which (n == rownames (j))])
    }
    j = joint.w.rownames (j)
    setNames (lapply (1:nrow (j), neigbors.of.xi), rownames (j))
}

######################################################################################
#######     Gets joints and queries for specific causal networks.              ######
######################################################################################

# One cause, one effect.
joint.c1e1.neighbors = neighbors.of.j (j.cgm.c1e1)
joint.sampler.c1e1 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.c1e1(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.c1e1.neighbors)[2,]
    return(dist)
}

# Common cause with two effects.
joint.c1e2.neighbors = neighbors.of.j (j.cgm.c1e2)
joint.sampler.c1e2 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.c1e2(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.c1e2.neighbors)[2,]
    return(dist)
}
samplerQuery.c1e2 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.c1e2.joint (joint.sampler.c1e2 (c,m,b,chainLen,nChains,bias,propStruct))

# Common cause with two effects (with names = x, ya, and yb).
joint.c1e2.xy.neighbors = neighbors.of.j (j.cgm.c1e2.xy)
joint.sampler.c1e2.xy <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.c1e2.xy(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.c1e2.xy.neighbors)[2,]
    return(dist)
}

# Common effect.
joint.ic2e1.neighbors = neighbors.of.j (j.cgm.c2e1)
joint.sampler.ic2e1 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ic2e1(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.ic2e1.neighbors)[2,]
    return(dist)
}
samplerQuery.ic2e1 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.c2e1.joint (joint.sampler.ic2e1 (c,m,b,chainLen,nChains,bias,propStruct))

# Common effect (with names = x, ya, and yb).
joint.c2e1.xy.neighbors = neighbors.of.j (j.cgm.c2e1.xy)
joint.sampler.ic2e1.xy <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ic2e1.xy(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruc,joint.c2e1.xy.neighborst)[2,]
    return(dist)
}

joint.sampler.ic2e1.full <- function (ca,cb,ma,mb,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ic2e1.full(ca,cb,ma,mb,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct)[2,]
    return(dist)
}

# Three element chains.
joint.ch3.neighbors = neighbors.of.j (j.cgm.ch3)
joint.sampler.ch3 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ch3(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.ch3.neighbors)[2,]
    return(dist)
}
samplerQuery.ch3 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.ch3.joint (joint.sampler.ch3 (c,m,b,chainLen,nChains,bias,propStruct))

# Common cause with three effects.
joint.c1e3.neighbors = neighbors.of.j (j.cgm.c1e3)
joint.sampler.c1e3 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.c1e3(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.c1e3.neighbors)[2,]
    return(dist)
}
samplerQuery.c1e3 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.c1e3.joint (joint.sampler.c1e3 (c,m,b,chainLen,nChains,bias,propStruct))

# Common effect with three causes.
joint.ic3e1.neighbors = neighbors.of.j (j.cgm.ic3e1)
joint.sampler.ic3e1 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ic3e1(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.ic3e1.neighbors)[2,]
    return(dist)
}
samplerQuery.ic3e1 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.ic3e1.joint (joint.sampler.ic3e1 (c,m,b,chainLen,nChains,bias,propStruct))

# Four element chains.
joint.ch4.neighbors = neighbors.of.j (j.cgm.ch4)
joint.sampler.ch4 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ch4(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.ch4.neighbors)[2,]
    return(dist)
}
samplerQuery.ch4 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.ch4.joint (joint.sampler.ch4 (c,m,b,chainLen,nChains,bias,propStruct))

# A 122 structure
joint.122.neighbors = neighbors.of.j (j.cgm.122)
joint.sampler.122 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.122(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.122.neighbors)[2,]
    return(dist)
}
samplerQuery.122 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.122.joint (joint.sampler.122 (c,m,b,chainLen,nChains,bias,propStruct))

# A 221 structure.
joint.221.neighbors = neighbors.of.j (j.cgm.221)
joint.sampler.221 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.221(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.221.neighbors)[2,]
    return(dist)
}
samplerQuery.221 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.221.joint (joint.cgm.221 (c,m,b,chainLen,nChains,bias,propStruct))

# A 212 structure.
joint.212.neighbors = neighbors.of.j (j.cgm.212)
joint.sampler.212 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.212(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.212.neighbors)[2,]
    return(dist)
}
samplerQuery.212 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.212.joint (joint.cgm.212 (c,m,b,chainLen,nChains,bias,propStruct))

# A 311 structure.
joint.311.neighbors = neighbors.of.j (j.cgm.311)
joint.sampler.311 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.311(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.311.neighbors)[2,]
    return(dist)
}
samplerQuery.311 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.311.joint (joint.cgm.311 (c,m,b,chainLen,nChains,bias,propStruct))

# A 113 structure.
joint.113.neighbors = neighbors.of.j (j.cgm.113)
joint.sampler.113 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.113(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.113.neighbors)[2,]
    return(dist)
}
samplerQuery.113 <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.113.joint (joint.cgm.113 (c,m,b,chainLen,nChains,bias,propStruct))

# Student network.
joint.student.neighbors = neighbors.of.j (j.cgm.student)
joint.sampler.student <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.student(c,m,b)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.student.neighbors)[2,]
    return(dist)
}
samplerQuery.student <- function (c,m,b,chainLen,nChains,bias=T,propStruct=T) 
    query.cgm.student.joint (joint.cgm.student (c,m,b,chainLen,nChains,bias,propStruct))

# A ic2e1 + cc2e1 network. 
joint.ic2e1.cc2e1.neighbors = neighbors.of.j (j.cgm.c2e1.c2e1)
joint.sampler.ic2e1.cc2e1 <- function (c1,m1,b1,c2,m2,b2,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ic2e1.cc2e1(c1, m1, b1, c2, m2, b2)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.ic2e1.cc2e1.neighbors)[2,]
    return(dist)
}

# Double ic2e1 network. 
joint.ic2e1.ic2e1.neighbors = neighbors.of.j (j.cgm.c2e1.c2e1)
joint.sampler.ic2e1.ic2e1 <- function (c1,m1,b1,c2,m2,b2,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.ic2e1.ic2e1(c1, m1, b1, c2, m2, b2)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.ic2e1.ic2e1.neighbors)[2,]
    return(dist)
}

# Double cc2e1 network. 
joint.cc2e1.cc2e1.neighbors = neighbors.of.j (j.cgm.c2e1.c2e1)
joint.sampler.cc2e1.cc2e1 <- function (c1,m1,b1,c2,m2,b2,chainLen,nChains,bias=T,propStruct=T) {
    dist   <- joint.cgm.cc2e1.cc2e1(c1, m1, b1, c2, m2, b2)
    dist$p <- jlMaster(dist,chainLen,nChains,bias,propStruct,joint.cc2e1.cc2e1.neighbors)[2,]
    return(dist)
}

########################################################################################
#######            Returns expected joints without sampling                       ######
########################################################################################
jl.transitions.in.parallel = FALSE
joint.sampler.generic = function (j, chainLen, neighbors=NULL, d=.50, bias=T, method='mutate', q.ps=NULL, chain.fun) {
    nr = nrow (j)
    # Each row xi has transition probabilities for state j[xi,] 
    j = joint.normalized (j)
    transition.ps = transition.ps (j, neighbors, method, q.ps)
    if (bias) {
        # Start with probability mass on prototypes.
        p.min = 1e-10; ps = rep (p.min, nr)
        p.prototype = (1 - (nr - 2) * p.min) # Mass to distribute between 111 and 000
        ps[nr] = d * p.prototype             # Mass to 111
        ps[1]  = (1 - d) * p.prototype       # Mass to 000
    }
    else
        ps = rep (1 / nr, nr)
    j$p = as.vector (do.call (chain.fun, list (ps, transition.ps)))
    j
}

joint.sampler.exact = function (j, chainLen, neighbors=NULL, d=.50, bias=T, method='mutate', q.ps=NULL) {
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
    if (chainLen == Inf) return (j)
    j = joint.sampler.generic (j, chainLen, neighbors, d, bias, method, q.ps, chain.fun=chain.fun)
    j
}

# Lambda is parameter of a Poisson distribution over chain lengths. 
# Since chain lengths are at least 2, so too must be lambda
# (i.e., actual Poisson distributin parameter is lambda - 2).
joint.sampler.poisson = function (j, lambda, neighbors=NULL, d=.50, bias=T, method='mutate', q.ps=NULL) {
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
        ps.weighted / sum (ps.weighted, na.rm=T) # Just in case.
    }
    joint.sampler.generic (j, lambda, neighbors, d, bias, method, q.ps, chain.fun=chain.fun)
}

conditional.proposal.distribution = function (j, gx, neighbors, condition.true.weight=10) {
    ps1 = function (nis) {
        ps = rep (0, nrow (j))
        ps[nis] = sapply (nis, function (ni) if (j[ni,]$p == 0) 1 else condition.true.weight)
        ps = ps / sum (ps)
    }
    if (!is.null (gx))
        j = joint.conditionalized (j, gx)
    laply (neighbors, ps1)
}

# One cause, one effect.
joint.sampler.c1e1.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e1.neighbors else NULL
    joint.sampler.exact (joint.cgm.c1e1 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.c1e1.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e1.neighbors else NULL
    joint.sampler.poisson (joint.cgm.c1e1 (c, m, b), lambda, neighbors, d, bias, method)
}

# Common cause with two effects.
joint.sampler.c1e2.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.neighbors else NULL
    joint.sampler.exact (joint.cgm.c1e2 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.c1e2.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.neighbors else NULL
    joint.sampler.poisson (joint.cgm.c1e2 (c, m, b), lambda, neighbors, d, bias, method)
}

# Common cause with two effects (with names = x, ya, and yb).
joint.sampler.c1e2.xy.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.xy.neighbors else NULL
    joint.sampler.exact (joint.cgm.c1e2.xy (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.c1e2.xy.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.xy.neighbors else NULL
    joint.sampler.poisson (joint.cgm.c1e2.xy (c, m, b), lambda, neighbors, d, bias, method)
}

# Common effect with independent causes.
joint.sampler.ic2e1.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic2e1.neighbors else NULL
    joint.sampler.exact (joint.cgm.ic2e1 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.ic2e1.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.xy.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ic2e1 (c, m, b), lambda, neighbors, d, bias, method)
}

# Common effect with independent causes (with names = x, ya, and yb).
joint.sampler.ic2e1.xy.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c2e1.xy.neighbors else NULL
    joint.sampler.exact (joint.cgm.ic2e1.xy (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.ic2e1.xy.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.xy.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ic2e1.xy (c, m, b), lambda, neighbors, d, bias, method)
}

# Common effect with conjunctive causes (with names = x, ya, and yb).
joint.sampler.cc2e1.xy.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c2e1.xy.neighbors else NULL
    joint.sampler.exact (joint.cgm.cc2e1.xy (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.cc2e1.xy.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e2.xy.neighbors else NULL
    joint.sampler.poisson (joint.cgm.cc2e1.xy (c, m, b), lambda, neighbors, d, bias, method)
}

# Three element chain.
joint.sampler.ch3.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ch3.neighbors else NULL
    joint.sampler.exact (joint.cgm.ch3 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.ch3.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ch3.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ch3 (c, m, b), lambda, neighbors, d, bias, method)
}

# Common cause with three effects.
joint.sampler.c1e3.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e3.neighbors else NULL
    joint.sampler.exact (joint.cgm.c1e3 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.c1e3.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.c1e3.neighbors else NULL
    joint.sampler.poisson (joint.cgm.c1e3 (c, m, b), lambda, neighbors, d, bias, method)
}

# Common effect with three causes
joint.sampler.ic3e1.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic3e1.neighbors else NULL
    joint.sampler.exact (joint.cgm.ic3e1 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.ic3e1.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic3e1.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ic3e1 (c, m, b), lambda, neighbors, d, bias, method)
}

# Four element chain.
joint.sampler.ch4.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ch4.neighbors else NULL
    joint.sampler.exact (joint.cgm.ch4 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.ch4.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ch4.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ch4 (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.122.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.122.neighbors else NULL
    joint.sampler.exact (joint.cgm.122 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.122.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.122.neighbors else NULL
    joint.sampler.poisson (joint.cgm.122 (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.221.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.221.neighbors else NULL
    joint.sampler.exact (joint.cgm.221 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.221.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.221.neighbors else NULL
    joint.sampler.poisson (joint.cgm.221 (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.212.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.212.neighbors else NULL
    joint.sampler.exact (joint.cgm.212 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.212.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.212.neighbors else NULL
    joint.sampler.poisson (joint.cgm.212 (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.311.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.311.neighbors else NULL
    joint.sampler.exact (joint.cgm.311 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.311.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.311.neighbors else NULL
    joint.sampler.poisson (joint.cgm.311 (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.113.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.113.neighbors else NULL
    joint.sampler.exact (joint.cgm.113 (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.113.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.113.neighbors else NULL
    joint.sampler.poisson (joint.cgm.113 (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.student.exact = function (c, m, b, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.student.neighbors else NULL
    joint.sampler.exact (joint.cgm.student (c, m, b), chainLen, neighbors, d, bias, method)
}
joint.sampler.student.poisson = function (c, m, b, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.student.neighbors else NULL
    joint.sampler.poisson (joint.cgm.student (c, m, b), lambda, neighbors, d, bias, method)
}

joint.sampler.ic2e1.cc2e1.exact = function (c1, m1, b1, c2, m2, b2, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic2e1.cc2e1.neighbors else NULL
    joint.sampler.exact (joint.cgm.ic2e1.cc2e1 (c1, m1, b1, c2, m2, b2), chainLen, neighbors, d, bias, method)
}
query.sampler.ic2e1.cc2e1.exact = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') 
    query.cgm.c2e1.c2e1.subjoints (
        joint.sampler.ic2e1.xy.exact (c1, m1, b1, lambda, d=d, propStruct=propStruct, bias=bias, method=method), 
        joint.sampler.cc2e1.xy.exact (c2, m2, b2, lambda, d=d, propStruct=propStruct, bias=bias, method=method)
    )
joint.sampler.ic2e1.cc2e1.poisson = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic2e1.cc2e1.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ic2e1.cc2e1 (c1, m1, b1, c2, m2, b2), lambda, neighbors, d, bias, method)
}
query.sampler.ic2e1.cc2e1.poisson = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') 
    query.cgm.c2e1.c2e1.subjoints (
        joint.sampler.ic2e1.xy.poisson (c1, m1, b1, lambda, d=d, propStruct=propStruct, bias=bias, method=method), 
        joint.sampler.cc2e1.xy.poisson (c2, m2, b2, lambda, d=d, propStruct=propStruct, bias=bias, method=method)
    )

joint.sampler.ic2e1.ic2e1.exact = function (c1, m1, b1, c2, m2, b2, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic2e1.ic2e1.neighbors else NULL
    joint.sampler.exact (joint.cgm.ic2e1.ic2e1 (c1, m1, b1, c2, m2, b2), chainLen, neighbors, d, bias, method)
}
query.sampler.ic2e1.ic2e1.exact = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') 
    query.cgm.c2e1.c2e1.subjoints (
        joint.sampler.ic2e1.xy.exact (c1, m1, b1, lambda, d=d, propStruct=propStruct, bias=bias, method=method), 
        joint.sampler.ic2e1.xy.exact (c2, m2, b2, lambda, d=d, propStruct=propStruct, bias=bias, method=method)
    )
joint.sampler.ic2e1.ic2e1.poisson = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.ic2e1.ic2e1.neighbors else NULL
    joint.sampler.poisson (joint.cgm.ic2e1.ic2e1 (c1, m1, b1, c2, m2, b2), lambda, neighbors, d, bias, method)
}
query.sampler.ic2e1.ic2e1.poisson = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') 
    query.cgm.c2e1.c2e1.subjoints (
        joint.sampler.ic2e1.xy.poisson (c1, m1, b1, lambda, d=d, propStruct=propStruct, bias=bias, method=method), 
        joint.sampler.ic2e1.xy.poisson (c2, m2, b2, lambda, d=d, propStruct=propStruct, bias=bias, method=method)
    )

joint.sampler.cc2e1.cc2e1.exact = function (c1, m1, b1, c2, m2, b2, chainLen, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.cc2e1.cc2e1.neighbors else NULL
    joint.sampler.exact (joint.cgm.cc2e1.cc2e1 (c1, m1, b1, c2, m2, b2), chainLen, neighbors, d, bias, method)
}
query.sampler.cc2e1.cc2e1.exact = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') 
    query.cgm.c2e1.c2e1.subjoints (
        joint.sampler.cc2e1.xy.exact (c1, m1, b1, lambda, d=d, propStruct=propStruct, bias=bias, method=method), 
        joint.sampler.cc2e1.xy.exact (c2, m2, b2, lambda, d=d, propStruct=propStruct, bias=bias, method=method)
    )
joint.sampler.cc2e1.cc2e1.poisson = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') {
    neighbors = if (propStruct) joint.cc2e1.cc2e1.neighbors else NULL
    joint.sampler.poisson (joint.cgm.cc2e1.cc2e1 (c1, m1, b1, c2, m2, b2), lambda, neighbors, d, bias, method)
}
query.sampler.cc2e1.cc2e1.poisson = function (c1, m1, b1, c2, m2, b2, lambda, d=.50, propStruct=T, bias=T, method='mutate') 
    query.cgm.c2e1.c2e1.subjoints (
        joint.sampler.cc2e1.xy.poisson (c1, m1, b1, lambda, d=d, propStruct=propStruct, bias=bias, method=method), 
        joint.sampler.cc2e1.xy.poisson (c2, m2, b2, lambda, d=d, propStruct=propStruct, bias=bias, method=method)
    )















#-------------------
# takes in causal strengths (ms), base rates (bs), 2 models, parameter values
# outputs evaluation metrics for all causal variables 
# ZD::: something wrong with genchainsMSpoislen?
#-------------------

# model1 and model2 need to be the function that takes in ms and bs and outputs IK conventions
testing_framework <- function (ms, bs, BMS_pars, betavar_pars, parvar_pars) {
    # joint.cgm.generic2 follows IK convention, starting at 111
    testjoint = joint.cgm.generic2(ms, bs) 
    
    # generate chains with lengths according poisson distr
    # calculate mean joint, and joint per chain, and vector chainlengths
    # function  going from chainjoints + chainlens to response distributions.
    BMS_out <- genchainsMSpoislen(meanlen=BMS_pars[['meanlen']], 
                                  nchains=BMS_pars[['nchains']], 
                                  bias=BMS_pars[['bias']], 
                                  joint=testjoint) %>% 
        genjoint %>% 
        genrespdistr(betavar=0) %>% 
        .$respdistr
    
    # calculates normative responses to cpjs and distorts them by some amount
    betavar_out <- betavariability.sim(joint=testjoint, 
                                       concentration=betavar_pars[['concentration']], 
                                       nSamps=betavar_pars[['nSamps']])
    
    # assumes optimal responding, but uncertainty about the parameters
    parvar_out <- param_variability(ms=ms, 
                                    bs=bs, 
                                    ms_conc=parvar_pars[['ms_conc']],
                                    bs_conc=parvar_pars[['bs_conc']],
                                    nSamples=parvar_pars[['nSamples']]) %>% 
        genrespdistr(betavar=0) %>% 
        .$respdistr
    
    # evaluating the similarity --------------------------
    browser()
}

# setting up causal parameters -----------------------
var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
ms['X1','Y'] = ms['Y','X2'] = .5
bs = c (.5, .25, .25)

# running the comparison
nSamples <- 10
testing_framework(ms=ms, 
                  bs=bs, 
                  BMS_pars=c(meanlen=10, nchains=nSamples, bias=.5), 
                  betavar_pars=c(concentration=100, nSamps=nSamples), 
                  parvar_pars=c(ms_conc=100, bs_conc=100, nSamples=nSamples))





















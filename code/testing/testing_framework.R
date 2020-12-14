#-------------------
# takes in causal strengths (ms), base rates (bs), 2 models, parameter values
# outputs evaluation metrics for all causal variables 
#-------------------
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
library(tidyverse)
library(ggplot2)

# spits out difference between distributions for each of the 3 models ----------------
testing_framework <- function (ms, bs, BMS_pars, betavar_pars, parvar_pars) {
    # joint.cgm.generic2 follows IK convention, starting at 111
    testjoint = joint.cgm.generic3(ms, bs) 
    
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
                                       nSamps=betavar_pars[['nSamps']]) %>% 
        .$respdistr
    
    # assumes optimal responding, but uncertainty about the parameters
    parvar_out <- param_variability(ms=ms, 
                                    bs=bs, 
                                    ms_conc=parvar_pars[['ms_conc']],
                                    bs_conc=parvar_pars[['bs_conc']],
                                    nSamples=parvar_pars[['nSamples']]) %>% 
        genrespdistr(betavar=0) %>% 
        .$respdistr
    
    # evaluating the similarity --------------------------
    similarity_mat <- array(NaN, dim=c(3, ncol(BMS_out)))
    rownames(similarity_mat) <- c('BMS:parvar', 'BMS:betavar', 'parvar:betavar')
    colnames(similarity_mat) <- colnames(BMS_out)
    for (judge_idx in 1:ncol(BMS_out)) {
        similarity_mat['BMS:parvar',judge_idx]     <- KLD(BMS_out[,judge_idx],    parvar_out[,judge_idx])$intrinsic.discrepancy
        similarity_mat['BMS:betavar',judge_idx]    <- KLD(BMS_out[,judge_idx],    betavar_out[,judge_idx])$intrinsic.discrepancy
        similarity_mat['parvar:betavar',judge_idx] <- KLD(parvar_out[,judge_idx], betavar_out[,judge_idx])$intrinsic.discrepancy
    }
    return(similarity_mat)
}



# setting up causal parameters -----------------------
var.names = c ('X1','Y','X2')
ms = matrix (0, nrow=3, ncol=3, dimnames = list (var.names, var.names))
plot_type = 'ce'
if (plot_type == 'chain') {
    ms['X1','Y'] = ms['Y','X2'] = .25 #chain
} else if (plot_type == 'cc') {
    ms['Y','X1'] = ms['Y','X2'] = .25 # common cause
} else if (plot_type == 'ce') {
    ms['X1','Y'] = ms['X2','Y'] = .25 # common effect
}
bs = c (.5, .5, .5)

# running the comparison
nSamples <- 1000
judgment_compare <- testing_framework(ms=ms, 
                                      bs=bs, 
                                      BMS_pars=c(meanlen=10, nchains=nSamples, bias=.5), 
                                      betavar_pars=c(concentration=10, nSamps=nSamples), 
                                      parvar_pars=c(ms_conc=10, bs_conc=10, nSamples=nSamples))

# plotting heatmap of differences -----------------------
dd <- judgment_compare %>% 
    as_tibble %>% 
    mutate(model_comp = c('BMS:parvar', 'BMS:betavar', 'parvar:betavar')) %>% 
    pivot_longer(-model_comp, names_to = 'judgment', values_to = 'intrinsic_discrepancy')

ddplot <- dd %>% 
    ggplot() +
    geom_tile(aes(x=judgment, y=model_comp, fill=intrinsic_discrepancy)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90))

ggsave(filename = paste0('judge_comp_', plot_type, '.pdf'),
       plot = ddplot,
       path = '../../figures/judgment_compare/',
       width = 10,
       height = 4,
       units = 'in')
ggsave(filename = paste0('judge_comp_', plot_type, '.png'),
       plot = ddplot,
       path = '../../figures/judgment_compare/',
       width = 10,
       height = 4,
       units = 'in')
















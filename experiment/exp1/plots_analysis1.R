## plots and analysis exp 1

require(ggplot2)
require(tidyr)
require(ggridges)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("datexp1_clean.RData") # data from exp double graph format and inference incentivization.
#datlong is Df with cleaned data

################# Descriptive plots


ggplot(datlong)+ 
  geom_histogram(aes(x=response), binwidth=1)+
  facet_wrap(facets=vars(participant), nrow=5)

#Plot RTs
ggplot(datlong)+ 
  geom_histogram(aes(x=rt), binwidth=1)+
  facet_wrap(facets=vars(participant))


# response distributions per pp and inftype
ggplot(datlong, aes(x=inftype, y=response))+
  geom_violin(trim=F, scale="area", bw=5)+
  geom_dotplot(aes(color=normalqueried),
               binaxis='y', stackdir='center',
               stackratio=0.7, dotsize=0.8, alpha=.7)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red")+
  ylim(c(0,100))+
  facet_wrap(facets=vars(participant), nrow=3)


# ridge plots
# can put group=interaction(participant,normalqueried) to see effect normalquery
ggplot(datlong, aes(x=respflip, y=participant, group=participant, fill=participant))+
  geom_density_ridges(rel_min_height = 0.01, bandwidth=4, alpha=.9,
                      jittered_points = TRUE,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x))+
  theme_minimal()+
  facet_wrap(facets=vars(inftype), nrow=1)

ggplot(datlong, aes(x=respflip, y=participant, group=participant))+
  geom_density_ridges(rel_min_height = 0.01, bandwidth=4, alpha=.9,
                      aes(scale=3))+
  theme_minimal()+
  facet_wrap(facets=vars(inftype), nrow=1)


# response distributions per pp, blockindex, and inftype
ggplot(datlong, aes(x=inftype, y=response))+
  geom_violin(trim=F, scale="area", bw=5)+
  geom_dotplot(aes(color=normalqueried),
               binaxis='y', stackdir='center',
               stackratio=0.7, dotsize=2, alpha=.7)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red")+
  ylim(c(0,100))+
  facet_wrap(facets=vars(participant, blockindex), nrow=5)

# plot group level judgements
ggplot(datlong, aes(x=inftype, y=respflip))+
  geom_violin(trim=F, scale="area", bw=5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red")+
  ylim(c(0,100))



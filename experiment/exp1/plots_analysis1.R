## plots and analysis exp 1

require(ggplot2)
require(tidyr)
require(ggridges)
require(BayesFactor)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("datexp1_clean.RData") # data from exp double graph format and inference incentivization.
#datlong is Df with cleaned data
datlong$participant <- droplevels(datlong$participant)

#group data

datgroupnormal <-  group_by(datlong, participant, inftype,normalqueried) %>% 
  summarise(
    count = n(), 
    totalduration = mean(as.numeric(as.character(Q_TotalDuration))/60),
    meanrt1 = mean(Last.Click, na.rm = TRUE),
    medianrt1 = median(Last.Click, na.rm = TRUE),
    sdrt1 = sd(Last.Click, na.rm = TRUE),
    meanrt2 = mean(Page.Submit, na.rm = TRUE),
    sdrt2 = sd(Page.Submit, na.rm = TRUE),
    maxrt2 = max(Page.Submit, na.rm = TRUE),
    meanacc = mean(acc, na.rm = TRUE),
    medianacc = median(acc, na.rm = TRUE),
    sdacc = sd(acc, na.rm = TRUE),
    sdresp = sd(response, na.rm = TRUE),
    sdrespflip = sd(respflip, na.rm = TRUE),
  )

res1 <- anovaBF(sdrespflip ~ inftype*normalqueried + participant, whichRandom="participant", data=datgroupnormal,
        whichModels = 'top')
res1 # BF for including interaction inftype*normalqueried is 1/0.215=4.65, hence some evidence that normal queries significantly change SD per inftype. Dont know if relevant for now.

datgroup <-  group_by(datlong, participant, inftype) %>% 
  summarise(
    count = n(), 
    totalduration = mean(as.numeric(as.character(Q_TotalDuration))/60),
    meanrt1 = mean(Last.Click, na.rm = TRUE),
    medianrt1 = median(Last.Click, na.rm = TRUE),
    sdrt1 = sd(Last.Click, na.rm = TRUE),
    meanrt2 = mean(Page.Submit, na.rm = TRUE),
    sdrt2 = sd(Page.Submit, na.rm = TRUE),
    maxrt2 = max(Page.Submit, na.rm = TRUE),
    meanacc = mean(acc, na.rm = TRUE),
    medianacc = median(acc, na.rm = TRUE),
    meanrespflip = mean(respflip, na.rm=T),
    sdacc = sd(acc, na.rm = TRUE),
    sdresp = sd(response, na.rm = TRUE),
    sdrespflip = sd(respflip, na.rm = TRUE),
  )


res2 <- anovaBF(sdrespflip ~ inftype + participant, whichRandom="participant", data=datgroup)
res2 #BF for including inftype = 154335364, so very significant
samples1 <- posterior(res2, iterations=10000, columnFilter="^$")
plot(samples1)
samplesinftype <- samples1[,grep("^inftype", colnames(samples1), value=T)]

samplesinftype <- pivot_longer(data.frame(samplesinftype), cols=colnames(data.frame(samplesinftype)))

############ plot posteriors for effects of inftype
ggplot(samplesinftype, aes(x=value, y=name))+ 
  geom_density_ridges(rel_min_height = 0.01)+
  geom_vline(aes(xintercept=0), linetype=2)+
  ggtitle("Posterior distributions of coefficients")+
  ylab("Inference type")+
  theme_minimal()



############ plot with SDs per inftype, and with mean responses as well?
plot(sdrespflip ~ inftype, data=datgroup)
#this plot based on datgroup, hence improper averaging, but is able to show standard error of sd per pp
ggplot(datgroup, aes(x=inftype))+
  geom_bar(stat = "summary", fun.y = "mean", aes(y=sdrespflip))+
  stat_summary(aes(y=sdrespflip), fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="red")+
  geom_point(aes(y=meanrespflip/2), fun.y = "mean", stat = "summary",
             shape='-', size=18)+
  scale_y_continuous(
    name = "Standard Deviation",
    sec.axis = sec_axis(~.*2, name="Mean Response")
  )  

#this plot based on means per inftype, thus proper averaging
datinftype <-  group_by(datlong, inftype) %>% 
  summarise(
    count = n(), 
    meanrespflip = mean(respflip, na.rm=T),
    sdresp = sd(response, na.rm = TRUE),
    sdrespflip = sd(respflip, na.rm = TRUE),
  )

ggplot(datinftype, aes(x=inftype))+
  geom_bar(stat='identity',aes(y=sdrespflip))+
  geom_point(aes(y=meanrespflip/2), shape='-', size=18)+
  scale_y_continuous(
    name = "Standard Deviation",
    sec.axis = sec_axis(~.*2, name="Mean Response")
  )


############ scatter plot mean resp and sd resp per inftype per pp
ggplot(datgroup, aes(x=sdrespflip, y=meanrespflip))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  facet_wrap(facets=vars(inftype), nrow=1)
# this seems to look as if there is a negative corr for nondiagnostic/markov inferences, while theres a positive relationship for diagnostic infs


# maybe driven by people who consistently responded exactly 50% for diag resps, or had other signatures of bad data: try without maybe badpps (see secriptives_cleaning1.R)
maybebadpps <- c("2", "15", "11", "16", "22", "31", "32", "33" ,"34")
ggplot(datgroup[-which(datgroup$participant %in% maybebadpps),], aes(x=sdrespflip, y=meanrespflip))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  facet_wrap(facets=vars(inftype), nrow=1)
# shows kinda the same 




########### Vincentized plot, showing variability abstracting from the participant level per inference type (see Ratcliff, 1979)
### probably not do this at all.
vincentize <- function(responses, nbins, whichbinmean=seq(1,nbins)){
  #takes vector of responses, bins it into nbins equal sized bins, outputs the means of these bins (or of a single bin if whichbinmean)
  responses <- sort(responses)
  v1<- floor(rank(responses, ties.method="random") * nbins / (length(test1) + 1)) #4 here indicate nr of bins
  means <- as.vector(by(responses, v1, mean))
  means <- means[whichbinmean]
  return(means)
}





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



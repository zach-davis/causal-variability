# Descriptives and plots of pilot variability experiment

require(ggplot2)
require(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load("pilotdat1.RData") #pilot 1 is with single graph format and no incentivization
load("pilotdat2.RData") #pilot 2 is with double graph format and inference incentivization.
#load("pilotdat3.RData") #pilot 3 is with double graph format and no incentivization.



######## REMOVE many RT outliers first!?#########################


################# Descriptive stats

#completion times
completiontimesminutes <- as.numeric(as.character(ppdat$Q_TotalDuration))/60

#attention questions
rownames(aqdat)<- seq(1,dim(aqdat)[1])
aqsubmittedvars <- grep("d*aqt[234]" , colnames(aqdat), value=T)
aqdat[,aqsubmittedvars] <- as.numeric(as.character(unlist(aqdat[,aqsubmittedvars])))
wrongaqs <- 10-apply(aqdat[,aqsubmittedvars], 1, function(x) sum(is.na(x))) #number of attention check sets wrong, not including the third try. (ie counts views of 2 and 3 page with aqs), hence 10 is maximum wrong

completiontimesminutes
wrongaqs
#10 is max wrong


#summary stats per pp and inftype
sumstats1 <-  group_by(datlong, participant) %>% 
  summarise(
    count = n(), 
    meanrt1 = mean(Last.Click, na.rm = TRUE),
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
sumstats1

sumstats2 <-  group_by(datlong, inftype) %>% 
  summarise(
    count = n(), 
    meanrt1 = mean(Last.Click, na.rm = TRUE),
    sdrt1 = sd(Last.Click, na.rm = TRUE),
    meanrt2 = mean(Page.Submit, na.rm = TRUE),
    sdrt2 = sd(Page.Submit, na.rm = TRUE),
    meanacc = mean(acc, na.rm = TRUE),
    medianacc = median(acc, na.rm = TRUE),
    sdacc = sd(acc, na.rm = TRUE),
    sdresp = sd(response, na.rm = TRUE),
    sdrespflip = sd(respflip, na.rm = TRUE),
  )
sumstats2


################# Descriptive plots


ggplot(datlong)+ 
  geom_histogram(aes(x=response), binwidth=1)+
  facet_wrap(facets=vars(participant), nrow=1)

#Plot RTs
ggplot(datlong)+ 
  geom_histogram(aes(x=Last.Click), binwidth=1)+
  xlim(c(0,15))+
  facet_wrap(facets=vars(participant))

ggplot(datlong)+ 
  geom_histogram(aes(x=Page.Submit), binwidth=1)+
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
  facet_wrap(facets=vars(participant), nrow=2)


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



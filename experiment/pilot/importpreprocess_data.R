#import and preprocess data pilot 1
require(ggplot2)
require(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read.csv('datapilot1.csv')
dat <- dat[6:8,] #hardcoded for initial 3 particpants, rows 1-5 is other stuff and previews!!!!!!!!!!!!!!!!!!!!!!!!!!!!
oldcolnames <- names(dat)

#gather data on presentation order of inferences and blocks
ppvars <- c("ProlificID", "Q_TotalDuration", "FL_12_DO", "Domain1economicsinferences_DO", "Domain2biologyinferences_DO", "Domain3astronomyinferences_DO", "Domain4meteorologyinferences_DO","Domain5sociologyinferences_DO")
ppdat <- dat[, ppvars]

#gather data on participant level, ie duration, attention checks, ID
aqvars <- grep("d[1-5]aq", oldcolnames, value=T)
aqvars <- aqvars[!grepl("Click", aqvars)]
aqvars
aqdat <- dat[,c("ProlificID", "Q_TotalDuration",aqvars)]
View(aqdat) #third pp failed many aq's

#gather data on inference level
trialvarnames <- grep("_1$", oldcolnames, value=T) #names of columns with responses
trialdat <- dat[,c("ProlificID", trialvarnames)]
colnames(trialdat)<-sub("_1", "", colnames(trialdat))
trialdatlong <- pivot_longer(trialdat, cols=sub("_1", "", trialvarnames))
trialdatlong$domain <- sub("q..|q.","", trialdatlong$name)
trialdatlong$inf <- sub("d.","",trialdatlong$name)
trialdatlong$response <- as.numeric(as.character(trialdatlong$value))

trialtimernames <- grep("^d.t", oldcolnames, value=T) #names of columns with the timers for the inferenes
trialtimernames <- c(trialtimernames, grep("^d124", oldcolnames, value=T)) #timer d1t24 was wrongly name d124 in qualtrics, so added here
trialtimerdat <- dat[,c("ProlificID", trialtimernames)]
trialtimerdat[,-1] <- as.data.frame(sapply(trialtimerdat[,-1], function(x) as.numeric(as.character(x)))) #change vars into numeric

#loop to add timer data to judgement data
trialdatlong$First.Click <- NA
trialdatlong$Last.Click <- NA
trialdatlong$Page.Submit <- NA
trialdatlong$Click.Count <- NA
trialdatlong$inforder <- NA #var for when the inference was presented in block
timervarstoadd <- c("First.Click", "Last.Click", "Page.Submit", "Click.Count")
for (pp in ppdat$ProlificID){
  for (dom in grep("^Domain", colnames(ppdat), value=T)){
    orderinfs <- unlist(strsplit(toString(ppdat[ppdat$ProlificID==pp,dom]), "[|]"))
    for (i in seq(1,47,2)){
      timer <- orderinfs[i]
      timer <- paste0(timer, "_")
      inference <- orderinfs[i+1]
      
      trialdatlong[trialdatlong$ProlificID==pp&trialdatlong$name==inference,timervarstoadd] <- trialtimerdat[trialtimerdat$ProlificID==pp,grep(timer, colnames(trialtimerdat), value=T)]
      trialdatlong[trialdatlong$ProlificID==pp&trialdatlong$name==inference,'inforder'] <- (i+1)/2
      
    }
    
  }
}

# we have 6 inferencetypes
# each type occurs 4 times in each block, repitition due to asking for normal value and by switching X1 with X2
# non diagnostic:
#   1:P(Xi|Y=1,Xj=1) infs: 1,4,7,10
#   2:P(Xi|Y=1,Xj=?) infs: 2,5,8,11
#   3:(Xi|Y=1,Xj=0) infs: 3,6,9,12
# diagnostic:
#   4:P(Y|Xi=1,Xj=1) infs: 13,16,19,22
#   5:P(Y|Xi=1,Xj=?) infs: 14,17,20,23
#   6:P(Y|Xi=1,Xj=0) infs: 15,18,21,24

trialdatlong$inftype <- NA
for (i in dim(trialdatlong)[1]){ # not finished
  if (trialdatlong$inf[i] %in% c("q1","q4","q7","q10")) {trialdatlong$inftype[i]<-1}  
}



#trialdatlong is the dataframe we can use for analysis


# this shows that basically only the second participant properly understood and performed on the task, they were also the person taking half an hour for the exp
# and strong bias for responding at upper half of scale, and on increments of 25
ggplot(trialdatlong)+ 
  geom_histogram(aes(x=response), binwidth=1)+
  facet_wrap(facets=vars(ProlificID))

#first and third participant also seem to have responded after waiting quite a period
ggplot(trialdatlong)+ 
  geom_histogram(aes(x=Last.Click), binwidth=1)+
  facet_wrap(facets=vars(ProlificID))

ggplot(trialdatlong)+ 
  geom_histogram(aes(x=Page.Submit), binwidth=1)+
  facet_wrap(facets=vars(ProlificID))

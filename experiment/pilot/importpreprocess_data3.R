#import and preprocess data pilot3 

rm(list=ls())



require(ggplot2)
require(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read.csv('datapilot3.csv')
dat <- dat[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),] #for pilot2. rows 1-5 is other stuff and previews, row 6-10 is participants from pilot 1, row 11-16 is from pilot 2
dat$participant <- seq(1, dim(dat)[1])
oldcolnames <- names(dat)

#gather data on presentation order of inferences and blocks
ppvars <- c("ProlificID","participant", "Q_TotalDuration", "EndDate","FL_12_DO", "Domain1economicsinferences_DO", "Domain2biologyinferences_DO", "Domain3astronomyinferences_DO", "Domain4meteorologyinferences_DO","Domain5sociologyinferences_DO")
ppdat <- dat[, ppvars]

#gather data on participant level, ie duration, attention checks, ID
aqvars <- grep("d[1-5]aq", oldcolnames, value=T)
aqvars <- aqvars[!grepl("Click", aqvars)]
aqvars
aqdat <- dat[,c("ProlificID","participant", "Q_TotalDuration", "EndDate",aqvars)]
View(aqdat) #pp3 had most aqs wrong

#gather data on inference level
trialvarnames <- grep("_1$", oldcolnames, value=T) #names of columns with responses
trialdat <- dat[,c("ProlificID", "participant", trialvarnames)]
colnames(trialdat)<-sub("_1", "", colnames(trialdat))
datlong <- pivot_longer(trialdat, cols=sub("_1", "", trialvarnames))
datlong$domain <- sub("q..|q.","", datlong$name)
datlong$inf <- sub("d.","",datlong$name)
datlong$response <- as.numeric(as.character(datlong$value))

trialtimernames <- grep("^d.t", oldcolnames, value=T) #names of columns with the timers for the inferenes
trialtimernames <- c(trialtimernames, grep("^d124", oldcolnames, value=T)) #timer d1t24 was wrongly name d124 in qualtrics, so added here
trialtimerdat <- dat[,c("ProlificID", trialtimernames)]
trialtimerdat[,-1] <- as.data.frame(sapply(trialtimerdat[,-1], function(x) as.numeric(as.character(x)))) #change vars into numeric


############### check if this stuff still works with the double network per inf, as those image q's probably show up in the randomization order variable!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### easy solution: remove each third entry in orderinfs
#loop to add timer data to judgement data and to add index variables
datlong$First.Click <- NA
datlong$Last.Click <- NA
datlong$Page.Submit <- NA
datlong$Click.Count <- NA
datlong$inforder <- NA #var for when the inference was presented in block
timervarstoadd <- c("First.Click", "Last.Click", "Page.Submit", "Click.Count")
for (pp in ppdat$ProlificID){
  for (dom in grep("^Domain", colnames(ppdat), value=T)){
    orderinfs <- unlist(strsplit(toString(ppdat[ppdat$ProlificID==pp,dom]), "[|]"))
    orderinfs <- orderinfs[-seq(3,72,3)]
    for (i in seq(1,47,2)){
      timer <- orderinfs[i]
      timer <- paste0(timer, "_")
      inference <- orderinfs[i+1]
      datlong[datlong$ProlificID==pp&datlong$name==inference,timervarstoadd] <- trialtimerdat[trialtimerdat$ProlificID==pp,grep(timer, colnames(trialtimerdat), value=T)]
      datlong[datlong$ProlificID==pp&datlong$name==inference,'inforder'] <- (i+1)/2
    }
  }
}


datlong$blockindex <- NA #var indexing order of domains
dommatch <- matrix(c("d1","d2","d3","d4","d5","FL_21", "FL_43", "FL_44", "FL_27", "FL_28"), nrow=2, byrow=T)
for (pp in ppdat$ProlificID){
  orderblocks <- unlist(strsplit(toString(ppdat[ppdat$ProlificID==pp,"FL_12_DO"]), "[|]"))
  for (blockind in  seq(1:5)){
    curblock <- orderblocks[blockind]
    curblock <- dommatch[1,which(dommatch[2,]==curblock)]
    datlong[datlong$ProlificID==pp&datlong$domain==curblock, "blockindex"] <- blockind
  }
}
datlong$trialindex <- datlong$inforder + (datlong$blockindex-1)*24  #var indexing order of all trials
# we have 6 inferencetypes
# each type occurs 4 times in each block, repitition due to asking for normal value and by switching X1 with X2
# normative answers calculated using joint.cgm.generic3 with ms['Y', 'X1'] = ms['Y','X2'] = .75 and bs = c (.2, .5, .2). (which aligns with instructions)
# non diagnostic:
#   1:P(Xi|Y=1,Xj=1) infs: 1,4,7,10, norm ans: 0.8
#   2:P(Xi|Y=1,Xj=?) infs: 2,5,8,11, norm ans: 0.8
#   3:(Xi|Y=1,Xj=0) infs: 3,6,9,12, norm ans: 0.8
# diagnostic:
#   4:P(Y|Xi=1,Xj=1) infs: 13,16,19,22, norm ans: 0.9411765 
#   5:P(Y|Xi=1,Xj=?) infs: 14,17,20,23, norm ans: 0.8
#   6:P(Y|Xi=1,Xj=0) infs: 15,18,21,24, norm ans: 0.5
datlong$inftype <- NA
datlong$inftype[which(datlong$inf %in% c("q1","q4","q7","q10"))] <- 1
datlong$inftype[which(datlong$inf %in% c("q2","q5","q8","q11"))] <- 2
datlong$inftype[which(datlong$inf %in% c("q3","q6","q9","q12"))] <- 3
datlong$inftype[which(datlong$inf %in% c("q13","q16","q19","q22"))] <- 4
datlong$inftype[which(datlong$inf %in% c("q14","q17","q20","q23"))] <- 5
datlong$inftype[which(datlong$inf %in% c("q15","q18","q21","q24"))] <- 6
datlong$inftype <- factor(datlong$inftype,
                                  levels = c(1,2,3,4,5,6),
                                  labels = c("P(Xi|Y=1,Xj=1)", "P(Xi|Y=1)", "P(Xi|Y=1,Xj=0)", "P(Y|Xi=1,Xj=1)", "P(Y|Xi=1)", "P(Y|Xi=1,Xj=0)"))

#variables for which the normal value is queried: 4-6, 10-12, 16-18, 22-24
datlong$normalqueried <- 0
datlong$normalqueried[which(datlong$inf %in% c("q4","q5","q6","q10","q11","q12","q16","q17","q18","q22","q23","q24"))] <- 1
datlong$normalqueried <- as.factor(datlong$normalqueried)

#response variable where normalqueried inferences are flipped around midpoint 50% (ie using P(x=1)=1-P(x=0))
datlong$respflip <- NA
datlong$respflip <- datlong$response
datlong$respflip[which(datlong$normalqueried==1)]<- 100-datlong$response[which(datlong$normalqueried==1)]

#variable indicating flipped normative answer (answers for normal queried inferences are flipped around midpoint 50%)
datlong$normflip <- NA
datlong$normflip[which(datlong$inf %in% c("q1","q4","q7","q10","q2","q5","q8","q11","q3","q6","q9","q12"))] <- 80
datlong$normflip[which(datlong$inf %in% c("q13","q16","q19","q22"))] <- 94
datlong$normflip[which(datlong$inf %in% c("q14","q17","q20","q23"))] <- 80
datlong$normflip[which(datlong$inf %in% c("q15","q18","q21","q24"))] <- 50

#variable indicating normative answer
datlong$normresp <- NA
datlong$normresp <- datlong$normflip
datlong$normresp[which(datlong$normalqueried==1)]<- 100-datlong$normflip[which(datlong$normalqueried==1)]

#variable indicating accuracy (response distance from normative answer)
datlong$acc <- NA
datlong$acc <- datlong$normresp - datlong$response
datlong$acc <- abs(datlong$acc)
  
# save the dataframes

#save(datlong, ppdat, aqdat, dat, file="pilotdat3.RData") #this was first 5 participants without incentivization and only 1 graph


# Descriptives and plots of pilot variability experiment

rm(list=ls())


require(ggplot2)
require(tidyr)
require(ggridges)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("datexp1.RData") #is with double graph format and inference incentivization.




######## REMOVE many RT outliers first!?#########################


################# Descriptive stats

#attention questions
rownames(aqdat)<- seq(1,dim(aqdat)[1])
aqsubmittedvars <- grep("d*aqt[234]" , colnames(aqdat), value=T)
aqdat[,aqsubmittedvars] <- as.numeric(as.character(unlist(aqdat[,aqsubmittedvars])))
wrongaqs <- 10-apply(aqdat[,aqsubmittedvars], 1, function(x) sum(is.na(x))) #number of attention check sets wrong, not including the third try. (ie counts views of 2 and 3 page with aqs), hence 10 is maximum wrong


wrongaqs
#10 is max wrong
#pps 1, 3, 6(+1end) and 7 have 2+ wrong. 2 and 15 have 2 wrong.
#pps 21, 32, 33, 37 have 4+ wrong, 23,26,28,29 have 3+ wrong, 19 has 2 wrong


#summary stats per pp and inftype
sumstats1 <-  group_by(datlong, participant) %>% 
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
sumstats1

#Fastest meanrt: pp3 meanrt=2, pp7 meanrt=4, pp2 meanrt=5, pp1 meanrt=6
#fastest mean rt: 28 and 21 <2, pp 19 has 2.9 median, pp 37 has 3.4 median.

#pp15 has meanrt=29, pp18 meanrt=23, pp14 meanrt=19
#Fastes totduration: pp1 26m, pp2 28m, pp12 31m, pp3 32m.

ggplot(datlong[datlong$participant=="15",])+
  geom_density(aes(x=Last.Click), bw=1)

ggplot(datlong[datlong$participant=="14",])+
  geom_density(aes(x=Last.Click), bw=1)


max(datlong$Last.Click) #=944, by pp14


#sdresflip>sdresp is sign of not distinguish normal and nonnormal values.

sumstats2 <-  group_by(datlong, inftype) %>% 
  summarise(
    count = n(), 
    meanrt1 = mean(Last.Click, na.rm = TRUE),
    medianrt1 = median(Last.Click, na.rm = TRUE),
    sdrt1 = sd(Last.Click, na.rm = TRUE),
    meanrt2 = mean(Page.Submit, na.rm = TRUE),
    sdrt2 = sd(Page.Submit, na.rm = TRUE),
    meanacc = mean(acc, na.rm = TRUE),
    medianacc = median(acc, na.rm = TRUE),
    sdacc = sd(acc, na.rm = TRUE),
    meanrespflip = mean(respflip, na.rm = TRUE),
    sdresp = sd(response, na.rm = TRUE),
    sdrespflip = sd(respflip, na.rm = TRUE),
  )
sumstats2

ppdat$Q382


# response distributions per pp to check if sensitive to nonnormalquery
ggplot(datlong, aes(x=participant, y=respflip))+
  geom_violin(trim=F, scale="area", bw=5)+
  geom_dotplot(aes(color=normalqueried),
               binaxis='y', stackdir='center',
               stackratio=0.7, dotsize=0.4, alpha=.7)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red")+
  ylim(c(0,100))

###############################################################
# clean data, remove pps, remove outliers
###############################################################

#rejected on prolific, pps: 3, 7, 1, (wrong aqs, unreasonably short rts, not distinguishing between normal/nonnormal)
#rejected on prolific, pps 19,21 completiontime <17min, 28 median rt <2, 37 <4 median rt, not sensitive to non/normal, and many aqs wrong, 
badpps <- c("1","3","6","7", "19", "21", "28", "37")
datlong <- datlong[-which(datlong$participant %in% badpps),]

# maybe remove these?
# pp2: 2 aqs wrong. pp2 has meanrt=5, totaldur=28, responded to 4 inftypes with just 50%, and says sometimes she confues normal-nonnormal,
# pp15: 2 aqs wrong, has meanrt=29, totaldur=119
# pp11: did not respond differently to nonnormal vs normal q's for conflict type infs(where conditioning on 1,0)!!!!  weird response profile per inftype, 
# pp16&34 (&25 with a bit of noise): responded 50% to all diagnostic questions, and 25/75% to all non-diagnostic questions !!
# pp22: responded 100% for nonnormal queries, 50% for almost all normal queries
# pp31: did kinda 25/75% for all infs except one, where he did 50%
# pp32&33: 4+ aqs wrong, 23,26,29 3+ aqs wrong
maybebadpps <- c("2", "15", "11", "16", "22", "31", "32", "33" ,"34")



# remove rt outliers

#remove long rts (+2sd?) (keep response though, just fill in NA for rt)
outliersNA <- function(x){
  x[scale(x)>2]<-NA
  return(x)
}
datlong$rt <- NA
datlong$rt <- with(datlong, unlist(tapply(Last.Click, participant, outliersNA)))
#remove too fast rts (remove observation completely, cannot trust response)
datlong <- subset(datlong, rt>=3)

#scale RT within participants
datlong$rt.S <- NA
datlong$rt.S <- with(datlong, unlist(tapply(rt, participant, scale)))


save(datlong, ppdat, aqdat, dat, file="datexp1_clean.RData") #dat exp with 2 graphs and incentivization.

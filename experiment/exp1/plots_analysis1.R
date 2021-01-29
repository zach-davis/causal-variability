## plots and analysis exp 1

require(ggplot2)
require(tidyr)
require(ggridges)
require(ggpubr)
require(BayesFactor)
require(lme4)
require(lmerTest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("datexp1_clean.RData") # data from exp double graph format and inference incentivization.
#datlong is Df with cleaned data

### two new variables
datlong$participant <- droplevels(datlong$participant)

datlong$diagnostic <- NA #variable indicating diagnostic inference or not
datlong$diagnostic[which(datlong$inftype %in% c("P(Xi|Y=1,Xj=1)", "P(Xi|Y=1)", "P(Xi|Y=1,Xj=0)"))] <- 0
datlong$diagnostic[which(datlong$inftype %in% c("P(Y|Xi=1,Xj=1)", "P(Y|Xi=1)", "P(Y|Xi=1,Xj=0)"))] <- 1
datlong$diagnostic <- factor(datlong$diagnostic,
                             levels = c(0,1),
                             labels = c("Non-diagnostic", "Diagnostic"))
datlong$information <- NA #variable indicating type of information provided in inference, consistent, inconsistent, incomplete
datlong$information[which(datlong$inftype %in% c("P(Xi|Y=1,Xj=1)", "P(Y|Xi=1,Xj=1)"))] <- 1
datlong$information[which(datlong$inftype %in% c("P(Xi|Y=1)", "P(Y|Xi=1)"))] <- 2
datlong$information[which(datlong$inftype %in% c("P(Xi|Y=1,Xj=0)", "P(Y|Xi=1,Xj=0)"))] <- 3
datlong$information <- factor(datlong$information,
                             levels = c(1,2,3),
                             labels = c("Consistent", "Incomplete", "Inconsistent"))



### figure 1: group level distributions
p1 <- ggplot(datlong, aes(x=respflip, y=information))+ #can also change to y=inftype
  geom_density_ridges(rel_min_height = 0.01, bandwidth=4, alpha=.9,
                      aes(scale=1.5))+
  xlab("Response (in %)")+
  ylab("Inference type")+
  facet_grid(rows=vars(diagnostic), scales = "free_x", space = "free_x",
             switch="both")+
  theme_minimal(base_size = 25)+
  theme(panel.spacing = unit(-2.5, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text=element_text(hjust=.25),
        axis.title.y = element_blank())
p1
ggsave("figures/figure1.png", p1,  width = 30, height = 20, units = "cm")

#or: 
ggplot(datlong, aes(x=interaction(information, diagnostic), y=respflip))+
  geom_violin(rel_min_height = 0.01, bandwidth=4, alpha=.9,
                      aes(scale=3))+
  theme_minimal()



#group data

datgroup <-  group_by(datlong, participant, inftype, diagnostic, information) %>% 
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

################## Analysis with factors Diagnostic and Information

res1 <- anovaBF(sdrespflip ~ diagnostic*information + participant, whichRandom="participant", data=datgroup, whichModels = 'withmain')
res1 #so there is evidence for main effects and interaction
res1[4]/res1[3] #BF for interaction

res1.f <- lmer(sdrespflip ~ diagnostic*information + (1|participant), data=datgroup) #frequentist
anova(res1.f)



############ plot posteriors for effects of inftype. maybe use this plot?
samples1 <- posterior(res1, 1, iterations=10000, columnFilter="^$")
plot(samples1)
samples1inftype <- samples1[,grep("^diagn|^informat", colnames(samples1), value=T)]
samples1inftype <- pivot_longer(data.frame(samples1inftype), cols=colnames(data.frame(samples1inftype)))
ggplot(samples1inftype, aes(x=value, y=name))+ 
  geom_density_ridges(rel_min_height = 0.01)+
  geom_vline(aes(xintercept=0), linetype=2)+
  ggtitle("Posterior distributions of coefficients")+
  ylab("Inference type factors")+
  theme_minimal(base_size = 18)


############ plot with SDs per inftype, and with mean responses as well, plot 2

#this plot based on datgroup, hence improper averaging, but is able to show standard error of sd per pp
p2 <- ggplot(datgroup, aes(x=information))+
  geom_bar(stat = "summary", fun.y = "mean", aes(y=sdrespflip))+
  stat_summary(aes(y=sdrespflip), fun.data=mean_se, fun.args = list(mult=1), 
               geom="errorbar",color="black", width=.2, size=1)+
  geom_point(aes(y=meanrespflip/2), fun.y = "mean", stat = "summary",
             shape='-', size=30, stroke=20, color="grey40")+
  stat_summary(aes(y=meanrespflip/2), fun.data=mean_se, fun.args = list(mult=1), 
               geom="errorbar", color="black", width=.2, size=1)+
  scale_y_continuous(
    name = "Mean Standard Deviation",
    sec.axis = sec_axis(~.*2, name="Mean Response (in %)"))+
  facet_grid(. ~ diagnostic, switch = "x", scales = "free_x", space = "free_x")+
  theme_minimal(base_size = 23)+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank())
p2
ggsave("figures/figure2.png", p2,  width = 30, height = 20, units = "cm")


############ test and plot accuracy (mean distance from correct response) and sd resp per inftype per pp
cors <- by(datgroup, interaction(datgroup$diagnostic,datgroup$information), 
           FUN = function(X) cor.test(X$sdrespflip, X$meanacc, method = "pearson"))

# Default bayesian hypothesis test for correlation, employing JZS prior
# from: Wetzels & Wagenmakers (2012) A default Bayesian hypothesis test for correlations and partial correlations
jzs_corbf <- function(r ,n){
  int=function(r,n,g){
    (1+g)^((n-2)/2)*(1+(1-r^2)*g)^(-(n-1)/2)*
      g^(-3/2)*exp(-n/(2*g))};
  bf10=sqrt((n/2))/gamma(1/2)*integrate(int,lower=0,upper=Inf,r=r,n=n)$value;
  return(bf10)
  
}
corBFs<-data.frame(inftype=names(unlist(lapply(cors, "[[", c("estimate")))),
                   cor=as.vector(unlist(lapply(cors, "[[", c("estimate")))),
                   n=rep(29,6))
corBFs$BF10 <- as.character(format(mapply(jzs_corbf, corBFs$cor, corBFs$n),digits=3, scientific=T))
corBFs$information <- c("Consistent", "Consistent", "Incomplete", "Incomplete", "Inconsistent", "Inconsistent")
corBFs$diagnostic <- c("Non-diagnostic", "Diagnostic", "Non-diagnostic", "Diagnostic", "Non-diagnostic", "Diagnostic")
corBFs$BF10 <- paste("BF","=", corBFs$BF10)

p3 <- ggplot(datgroup, aes(x=sdrespflip, y=meanacc))+
  geom_point()+
  geom_smooth(method='lm', se=F, colour="grey40")+
  stat_cor(method = "pearson", label.x = 20, label.y = 42,
           aes(label = ..r.label..), size=5)+
  #facet_wrap(facets=vars(diagnostic, information), nrow=2)+
  facet_grid(rows=vars(information), cols=vars(diagnostic), switch="both")+
  geom_text(x=26.5, y=37, aes(label=BF10, fontface=2), size=4, data=corBFs)+
  theme_minimal(base_size = 23)+
  theme(panel.spacing = unit(-.5, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Mean accuracy")+
  xlab("Standard deviation")
p3
ggsave("figures/figure3.png", p3,  width = 20, height = 20, units = "cm")


### plot not separated by inftype

ggplot(datgroup, aes(x=sdrespflip, y=meanacc))+
  geom_point()+
  geom_smooth(method='lm', se=F, colour="grey40")+
  stat_cor(method = "pearson", label.x = 20, label.y = 45,
           label.sep="\n", size=5)+
  theme_minimal(base_size = 23)+
  theme(panel.spacing = unit(1, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Mean accuracy")+
  xlab("Standard deviation")

datoverallcorr <-  group_by(datlong, participant) %>% 
  summarise(
    meanacc = mean(acc, na.rm = TRUE),
    sdrespflip = sd(respflip, na.rm = TRUE)
  )

overallcor <- cor.test(datoverallcorr$sdrespflip, datoverallcorr$meanacc, method = "pearson")
BFoverallcor <- jzs_corbf(overallcor$estimate, 29)
BFoverallcor

### plot separated by information and diagnostic factors

p31 <- ggplot(datgroup, aes(x=sdrespflip, y=meanacc, group=diagnostic, colour=diagnostic))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  stat_cor(method = "pearson", label.x = 20,
           label.sep=",", size=5)+
  theme_minimal(base_size = 23)+
  theme(panel.spacing = unit(1, "lines"))+
  ylab("Mean accuracy")+
  xlab("Standard deviation")

p32 <- ggplot(datgroup, aes(x=sdrespflip, y=meanacc, group=information, colour=information))+
  geom_point()+
  geom_smooth(method='lm', se=F)+
  stat_cor(method = "pearson", label.x = 20,
           label.sep=",", size=5)+
  theme_minimal(base_size = 23)+
  theme(panel.spacing = unit(1, "lines"))+
  ylab("Mean accuracy")+
  xlab("Standard deviation")

ggarrange(p31, p32, nrow=2)


### plot markov violations
#Plotting all participants at once is a mess, need to group them over their standard deviations.

#create grouping variable
datpp <-  group_by(datlong[datlong$diagnostic=="Non-diagnostic",], participant) %>% 
  summarise(
    sdrespflipmrkv = sd(respflip, na.rm = TRUE),
  )  
datpp$sdgroup <- Hmisc::cut2(datpp$sdrespflipmrkv, g=3)
datpp$sdgroup <- factor(datpp$sdgroup,
                        levels=levels(datpp$sdgroup),
                        labels=c("Low variability", "Medium variability", "High variability"))
datgroup <- merge(datgroup, datpp, by="participant")

datgroupmean<-  group_by(datgroup, sdgroup, inftype, diagnostic, information) %>% 
  summarise(
    meangroupmean = mean(meanrespflip, na.rm = TRUE),
  )  
datgroup <- merge(datgroup, datgroupmean)

p4 <- ggplot(datgroup[datgroup$diagnostic=="Non-diagnostic",])+
  facet_wrap(facets=vars(sdgroup))+
  geom_line(aes(y=meanrespflip, x=information, group=participant),
            color='grey40')+
  #geom_point(aes(y=meanrespflip, x=information, group=participant))+
  geom_line(aes(y=meangroupmean, x=information, group=participant), 
           size =3)+
  stat_summary(aes(y=meanrespflip, x=information), fun.data=mean_se, fun.args = list(mult=1), 
               geom="errorbar", color="black", width=.2, size=1)+
  geom_segment(aes(y=80, yend=80, x=1, xend=3),
               color='grey40', linetype="dashed", size=2)+
  ggtitle("Markov violations (non-diagnostic inferences)")+
  ylab("Mean response (in %)")+
  theme_minimal(base_size = 24)+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 40, hjust=1))


p4
ggsave("figures/figure4.png", p4,  width = 30, height = 20, units = "cm")

## statistics influence variability on markov violations
## (this could be done better using mixed model, long data, and meansd as continuous covariate)

datlong <- merge(datlong, datpp, by="participant")## add sdgroup var to long data
res3 <- anovaBF(respflip ~ information*sdgroup + participant, whichRandom="participant", 
                data=datlong[datlong$diagnostic=="Non-diagnostic",], whichModels = 'withmain')
res3 #so there is evidence for main information and interaction
res3[4]/res3[3] #BF for interaction.

res1.f <- lmer(respflip ~ information*sdgroup + (1|participant), data=datlong[datlong$diagnostic=="Non-diagnostic",]) #frequentist
anova(res1.f)
summary(res1.f)




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



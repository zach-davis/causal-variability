# Plot relevant aspects of simulated responses


#----------------------------------------------------------
# set dirs and load required scripts
#----------------------------------------------------------
# # When sourcing this script
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)

# When running in console
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

plotdir <- file.path(dirname(dirname(getwd())), "plots")  #ugly, make standard directory setup

#----------------------------------------------------------
source ("../utilities/utilities.R")
source ("../utilities/joint.utilities.R")
source ("../utilities/cgm.R")
source ("../models/mutsampler.R")
source ("../models/param_variability.R")

source ("../models/BMS.R")
source ("../models/beta_variability_model.R")


library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
library(gtable)


#----------------------------------------------------------
# Generate plots for each inference based on two sets of response distributions to compare
#----------------------------------------------------------


# fnc to plot tables in ggplots with titles
myTableGrob <- function(data_dt, title_v, fontsize_v = 14){
  # Create custom table grob with title
  # Data.table that the grob will be made out of
  # Title for display
  # Fontsize for title. Default is 14
  
  ## Table
  table_grob <- tableGrob(data_dt, rows = rep('', nrow(data_dt)))
  ## Title
  title_grob <- textGrob(title_v, gp = gpar(fontsize = fontsize_v))
  ## Add title
  table_grob <- gtable_add_rows(table_grob, heights = grobHeight(title_grob) + unit(5,'mm'), pos = 0)
  table_grob <- gtable_add_grob(table_grob, title_grob, 1, 1, 1, ncol(table_grob), clip = "off")
}


#fnc to generate plot of input parameters and joint, so that it can be added to the other plots
plot.input.pars <- function(simresp, nrplot=1){
  # first input is simulated responses, ie output genrespdistr
  # second nrplot input to declare that plotting will involve two sets of inputs, then input 2 for second set will change title to 'red lines'
    nrplot <- nrplot
    if (nrplot==1){plotcolor<-"(black lines)"}
    if (nrplot==2){plotcolor<-"(red lines)"}
    nms.pos.inputs <- c("meanChainlen", "betavar", "nChains", "bias", "model", "concentration", "nSamps") #possible input parameters to plot
    actual.inputs <- match(nms.pos.inputs, names(simresp)) #returns NA if not present, otherwise index.
    actual.inputs <- actual.inputs[!is.na(actual.inputs)] #remove NAs
    testjoint <- testjoint
    bs <- simresp$bs
    ms <- simresp$ms
    plot.inputs <- simresp[actual.inputs]
    normjoint1 <- simresp$normjoint
    simjoint1 <- simresp$meanjoint
    dfinputplot <- data.frame(posind=seq(1, length(plot.inputs)),inputnames=names(plot.inputs), plottext=as.character(plot.inputs)) #forcce list of inputs into DF
    dfinputplot$plottext <- sapply(dfinputplot$plottext, function(x){gsub('),', '), \n', x)}) # input linebreaks at ), for the joint inputs, NOT NEEDED ANYMORE
    
    p1 <- ggplot(dfinputplot)+
      geom_blank()+
      geom_label(aes(x=-1, y=posind, label=inputnames))+
      geom_text(aes(x=1, y=posind, label=plottext))+
      theme_void()+
      xlim(c(-1.5, 8))+
      annotation_custom(myTableGrob(normjoint1, "input joint distr"), xmin=4, xmax=4, ymin=4, ymax=5)+
      annotation_custom(myTableGrob(ms, "causal strengths"), xmin=7, xmax=7, ymin=4, ymax=5)+
      annotation_custom(textGrob(paste("baserates", bs[1], bs[2], bs[3])), xmin=7, xmax=7, ymin=3.5, ymax=4)+
      annotation_custom(myTableGrob(simjoint1, "mean sim joint distr"), xmin=4, xmax=4, ymin=2, ymax=3)+
      ggtitle(paste("Input parameters", plotcolor))
    
    return(p1)
  }


# fnc generates 27 density plots to compare simulated response distributions.
compare.resp.distr.plots <- function(respdistr1, respdistr2){
  #function to input two response inference distributions and plot all densities to compare
  plots.compare.allinfs <- list()
  for (inf in 1:29){ #27 plots of inferences, 2 of input parameters of simulations
    plots.compare.allinfs[[inf]] <- local({
      inf <- inf
      
      if (inf==1){ #plot with input parameters 1
        respdistr <- respdistr1
        p1 <- plot.input.pars(respdistr)
      }
      if (inf==2){ #plot with input parameters 2
        respdistr <- respdistr2
        p1 <- plot.input.pars(respdistr, 2)
      }
      
      if (inf>2) {
          distr1 <- respdistr1$respdistr[,inf-2]
          distr2 <- respdistr2$respdistr[,inf-2]
          
          normresp1 <- respdistr1$normresps[[inf-2]]
          normresp2 <- respdistr2$normresps[[inf-2]]
          
          titleinf <- names(respdistr1$respdistr)[inf-2]
          p1 <- ggplot() + 
            geom_density(aes(x=distr1, y=..density..))+
            geom_density(aes(x=distr2, y=..density..), colour='red')+
            geom_vline(aes(xintercept=normresp1))+
            geom_vline(aes(xintercept=normresp2), colour='red')+
            xlim(c(0,1))+
            ggtitle(paste("response distribution", titleinf))
      }
      
      
      print(p1)
    })
  }
  return(plots.compare.allinfs)
}



#----------------------------------------------------------
# Fncs for plots to visualize markov violations and explaining away
#----------------------------------------------------------


#Markov trials for chain & common cause
nms.chain.mrkv = list()
nms.chain.mrkv[[1]] = c("X1|Y==1 & X2==1","X1|Y==1", "X1|Y==1 & X2==0")
nms.chain.mrkv[[2]] = c("X1|Y==0 & X2==1","X1|Y==0", "X1|Y==0 & X2==0")
nms.chain.mrkv[[3]] = c("X2|Y==1 & X1==1","X2|Y==1", "X2|Y==1 & X1==0")
nms.chain.mrkv[[4]] = c("X2|Y==0 & X1==1","X2|Y==0", "X2|Y==0 & X1==0")

#Markov trials for Common effect
nms.comeff.mrkv = list()
nms.comeff.mrkv[[1]] = c("X1|X2==1","X1", "X1|X2==0")
nms.comeff.mrkv[[2]] = c("X2|X1==1","X2", "X2|X1==0")

#Explaining away trials for Common effect
nms.comeff.expaw = list()
nms.comeff.expaw[[1]] = c("X1|Y==1 & X2==0","X1|Y==1", "X1|Y==1 & X2==1") #inequalities, first should be the largest
nms.comeff.expaw[[2]] = c("X2|Y==1 & X1==0","X2|Y==1", "X2|Y==1 & X1==1")



plot.mrkv <- function(respdistr, inf.names){
  nplots <- length(inf.names)
  plots.mrkv <- list()
  for (i in 1:(nplots+1)){
    plots.mrkv[[i]] <- local({
      i <- i
      if (i==1){
        respdistr <- respdistr
        p1 <- plot.input.pars(respdistr)
      }
      if (i>1){
        resps <- respdistr$respdistr[,inf.names[[(i-1)]]]
        resps <- stack(resps)
        p1 <- ggplot(resps, aes(x=ind, y=values)) + 
          geom_violin() +
          stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                       geom="errorbar", color="red", width=0.2) +
          stat_summary(fun.y=mean, geom="point", color="red", size=4) +
          ggtitle("Markov violations / Failure to explain away")
      }
      print(p1)
    })
  }
  return(plots.mrkv)
}


compare.plot.mrkv <- function(respdistr1, respdistr2, inf.names){
  plots1 <- list()
  plots2 <- list()
  nPlots <- length(inf.names)
  plots1 <- plot.mrkv(respdistr1, inf.names)
  plots2 <- plot.mrkv(respdistr2, inf.names)
  allplots <- list()
  allplots[1] <- plots1[1]
  allplots[2] <- plots2[1]
  for (i in 1:(nPlots)){
    allplots[[i+2]] <- ggarrange(plots1[[i+1]], plots2[[i+1]])
  }
  return(allplots)
}

compare.plot.mrkv2 <- function(respdistr1, respdistr2, inf.names){
  plots1 <- list()
  plots2 <- list()
  plots1 <- plot.mrkv(respdistr1, inf.names)
  plots2 <- plot.mrkv(respdistr2, inf.names)
  plots1 <- c(plots1, plots2)
  return(plots1)
}



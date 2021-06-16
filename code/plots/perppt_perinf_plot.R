#-------------------
# plotting per-participant, per-inference type variability
library(ggplot2)
library(ggridges)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
load('../../experiment/exp1/datexp1_clean.RData')

p = datlong %>% 
    filter(domain!='d1') %>% #, inf %in% c('q16','q13')) %>% 
    mutate(inftype = recode(inftype, 
                            `P(Xi|Y=1,Xj=1)`= 'Predictive\nConsistent',
                            `P(Xi|Y=1)`= 'Predictive\nIncomplete',
                            `P(Xi|Y=1,Xj=0)`= 'Predictive\nInconsistent',
                            `P(Y|Xi=1,Xj=1)`= 'Diagnostic\nConsistent',
                            `P(Y|Xi=1)`= 'Diagnostic\nIncomplete',
                            `P(Y|Xi=1,Xj=0)`= 'Diagnostic\nInconsistent')) %>% 
    ggplot(aes(x = respflip, y = participant)) +
    stat_density_ridges() +
    xlim(c(0,100)) +
    xlab('conditional probability judgment') +
    facet_wrap(~inftype, ncol=6) +
    theme_ridges() + 
    theme(
        panel.grid.major.y = element_line(),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x  = element_text(size = 15),
        #axis.text.y  = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
    ) +
    NULL
ggsave(filename='../../figures/variability_ridgeplot.png', 
       plot=p,
       width=18,
       height=7,
       units='in')

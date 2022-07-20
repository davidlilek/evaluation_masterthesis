library(ggplot2)
library(reshape)
library(chron)
library(RColorBrewer)

###########################
#
#
# load color palette
#
#
###########################

mypalette_red <-brewer.pal(9,"PuRd")
#reverse color to get the more intense for MBR
mypalette_red <- rev(mypalette_red[c(7:9)])
mypalette_blue <-brewer.pal(9,"PuBu")
mypalette_blue <- rev(mypalette_blue[c(7:9)])
color_manual <- rep(c(mypalette_red,mypalette_blue),3)


timing <- read.csv2("timing_experiment.csv")
timing_4plot <- melt(timing[-c(9:10),c(1,3,5,7)],
                     id=c("Evaluation.method"))
timing_4plot$threads <- rep(c(1,10,30),each=8)
timing_4plot$method <- rep(rep(c("MBR","noMBR"),each=4),3)
timing_4plot$value <- timing_4plot$value/10
timing_4plot$Evaluation.method <- as.factor(timing_4plot$Evaluation.method)
       
timing_4plot$Evaluation.method <- factor(timing_4plot$Evaluation.method, 
       levels = c("MBR",
                  "noMBR",
                  "MBR_pooled",
                  "noMBR_pooled",
                  "MBR_pooled_oldversion",
                  "noMBR_pooled_oldversion",
                  "MBR_pooled_secondpeptides",
                  "noMBR_pooled_secondpeptides"
                  ))


ggplot(data=timing_4plot, aes(x=threads, y=value, group = Evaluation.method, color = method)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("Run-time for 15 samples [min]") + 
  xlab("Used threads") +
  scale_color_manual(values = c(mypalette_red[1],mypalette_blue[1])) +
  theme(legend.position="bottom",
        legend.box="vertical",
        legend.margin=margin()) +
  facet_wrap(~ Evaluation.method, ncol = 2)  

  pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/timing_experiment.png"

ggsave(pth,
       width = 7,
       height = 7) 

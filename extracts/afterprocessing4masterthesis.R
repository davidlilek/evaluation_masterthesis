##########################################
#
# extracts
#
##########################################

library(ggpubr)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(corrplot)

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

############################
#
# plots for base evaluation
#
############################
plots_no_twopeptides <- readRDS("plots_twopeptides.rds")
plots_no_onepeptide <- readRDS("plots_onepeptide.rds")
sample_names <- c("_A","_B","_C","_D","_E","_F")




for (sample_name in sample_names){
  ggarrange(plots_onepeptide[[sample_name]],plots_twopeptides[[sample_name]],
            ncol = 1, nrow = 2,
            common.legend = TRUE, legend="bottom")
  pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun",sample_name,".png",sep="")
  ggsave(pth,
         width = 7,
         height = 7)  
}

### combine all

sample_name1 <- "_A"
sample_name2 <- "_B"
sample_name3 <- "_C"
sample_name4 <- "_D"
sample_name5 <- "_E"
sample_name6 <- "_F"

p <- ggarrange(plots_no_onepeptide[[sample_name1]]+rremove("ylab")+rremove("xlab"),plots_no_onepeptide[[sample_name2]]+rremove("ylab")+rremove("xlab"),plots_no_onepeptide[[sample_name3]]+rremove("ylab")+rremove("xlab"),
               plots_no_onepeptide[[sample_name4]]+rremove("ylab")+rremove("xlab"),plots_no_onepeptide[[sample_name5]]+rremove("ylab")+rremove("xlab"),plots_no_onepeptide[[sample_name6]]+rremove("ylab")+rremove("xlab"),
               plots_no_twopeptides[[sample_name1]]+rremove("ylab")+rremove("xlab"),plots_no_twopeptides[[sample_name2]]+rremove("ylab")+rremove("xlab"),plots_no_twopeptides[[sample_name3]]+rremove("ylab")+rremove("xlab"),
               plots_no_twopeptides[[sample_name4]]+rremove("ylab")+rremove("xlab"),plots_no_twopeptides[[sample_name5]]+rremove("ylab")+rremove("xlab"),plots_no_twopeptides[[sample_name6]]+rremove("ylab")+rremove("xlab"),
               labels = c("- Extract A", "- Extract B", "- Extract C",
                          "- Extract D", "- Extract E", "- Extract F",
                          "- Extract A", "- Extract B", "- Extract C",
                          "- Extract D", "- Extract E", "- Extract F"),
               ncol = 3, nrow = 4,
               common.legend = TRUE, legend="bottom",
               hjust = c(-1.6,-1.7,-1.6,-1.7),
               align = "hv",
               font.label = list(size = 14, color = "black", face = "plain"))
annotate_figure(p, left = textGrob("No of protein identifications", rot = 90, hjust = 0.25, vjust = 1, gp = gpar(cex = 1.0)))
pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_compare_all.png",sep="")
ggsave(pth,
       width = 8.5,
       height = 10)

### combine A-B

sample_name1 <- "_A"
sample_name2 <- "_B"
p <- ggarrange(plots_no_onepeptide[[sample_name1]]+rremove("ylab"),plots_no_twopeptides[[sample_name1]]+rremove("ylab"), 
          plots_no_onepeptide[[sample_name2]]+rremove("ylab"),plots_no_twopeptides[[sample_name2]]+rremove("ylab"), 
          labels = c("- Extract A", "- Extract A", "- Extract B", "- Extract B"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom",
          hjust = c(-1.6,-1.7,-1.6,-1.7),
          align = "hv",
          font.label = list(size = 14, color = "black", face = "plain"))
annotate_figure(p, left = textGrob("No of protein identifications", rot = 90, hjust = 0.25, vjust = 1, gp = gpar(cex = 1.0)))
pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_compare",sample_name1,sample_name2,".png",sep="")
ggsave(pth,
       width = 7,
       height = 7)  

### combine C-D

sample_name1 <- "_C"
sample_name2 <- "_D"
p <- ggarrange(plots_no_onepeptide[[sample_name1]]+rremove("ylab"),plots_no_twopeptides[[sample_name1]]+rremove("ylab"), 
               plots_no_onepeptide[[sample_name2]]+rremove("ylab"),plots_no_twopeptides[[sample_name2]]+rremove("ylab"), 
               labels = c("- Extract C", "- Extract C", "- Extract D", "- Extract D"),
               ncol = 2, nrow = 2,
               common.legend = TRUE, legend="bottom",
               hjust = c(-1.6,-1.7,-1.6,-1.7),
               align = "hv",
               font.label = list(size = 14, color = "black", face = "plain"))
annotate_figure(p, left = textGrob("No of protein identifications", rot = 90, hjust = 0.25, vjust = 1, gp = gpar(cex = 1.0)))
pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_compare",sample_name1,sample_name2,".png",sep="")
ggsave(pth,
       width = 7,
       height = 7) 

### combine E-F

sample_name1 <- "_E"
sample_name2 <- "_F"
p <- ggarrange(plots_no_onepeptide[[sample_name1]]+rremove("ylab"),plots_no_twopeptides[[sample_name1]]+rremove("ylab"), 
               plots_no_onepeptide[[sample_name2]]+rremove("ylab"),plots_no_twopeptides[[sample_name2]]+rremove("ylab"), 
               labels = c("- Extract E", "- Extract E", "- Extract F", "- Extract F"),
               ncol = 2, nrow = 2,
               common.legend = TRUE, legend="bottom",
               hjust = c(-1.6,-1.7,-1.6,-1.7),
               align = "hv",
               font.label = list(size = 14, color = "black", face = "plain"))
annotate_figure(p, left = textGrob("No of protein identifications", rot = 90, hjust = 0.25, vjust = 1, gp = gpar(cex = 1.0)))
pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_compare",sample_name1,sample_name2,".png",sep="")
ggsave(pth,
       width = 7,
       height = 7) 

############################
#
# compare MBR no MBR
#
############################

# load and arrange data for one peptide
sample_names <- paste("Extract",LETTERS[1:6],sep=" ")
res_2gether_onepeptide <- readRDS("results_onepeptide.rds")
res_2gether_onepeptide <- as.data.frame(do.call(rbind,res_2gether_onepeptide))
res_2gether_onepeptide$sample_name <- rep(sample_names,each=18)

MBR_comparison <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^noMBR"))
summary(MBR_comparison$rsd)
summary(noMBR$rsd)
MBR_comparison$comparison <- ((MBR_comparison$no/noMBR$no)-1)*100
MBR_comparison_onepeptide <- MBR_comparison

# load and arrange data for two peptides
res_2gether_twopeptides <- readRDS("results_twopeptides.rds")
res_2gether_twopeptides <- as.data.frame(do.call(rbind,res_2gether_twopeptides))
res_2gether_twopeptides$sample_name <- rep(sample_names,each=18)

MBR_comparison <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^noMBR"))
summary(MBR_comparison$rsd)
summary(noMBR$rsd)
MBR_comparison$comparison <- ((MBR_comparison$no/noMBR$no)-1)*100
MBR_comparison_twopeptides <- MBR_comparison


# bring them 2gether
MBR_comparison <- rbind(MBR_comparison_onepeptide,MBR_comparison_twopeptides)
MBR_comparison$peptides <- rep(c("one","two"),each=54)

# correlation plot
pairs(MBR_comparison[,c(1,2,6)])
corrplot(cor(MBR_comparison[,c(1,2,6)]),method = 'color', order = 'alphabet')
# how much jitter on the x-axis
dodge <- position_dodge(.3)

# plot results
ggplot(MBR_comparison, aes(x=Type, y=comparison, color=Evaluation, shape=peptides)) +
  geom_point(position = dodge, size = 2.5) +
  scale_shape_manual(values=c(1, 16)) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Proportion of more proteins\nfound using MBR [%]") +
  labs(fill = "Evaluation \n method") + 
  facet_wrap(~sample_name, ncol=2) +
  scale_fill_brewer(palette="Dark2", labels = c("pooled\nMaxQuant","pooled\nMaxQuant oldversion","manually\npooled","manually\npooled")) +
  scale_color_manual(values = c(mypalette_red)) 

pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_comparisonmbr.png")
ggsave(pth,
       width = 7,
       height = 7) 

#################
#
# compare one and two peptides
#
################

# number of proteins

one <- readRDS("results_onepeptide.rds")
two <- readRDS("results_twopeptides.RDS")

one <- as.data.frame(one)
two <- as.data.frame(two)

diff_one_two <- ((one/two)-1)*100
diff_one_two <- diff_one_two[,grepl("no",colnames(diff_one_two))]
diff_one_two$Type <- one$X_A.Type
diff_one_two$Evaluation <- one$X_A.Evaluation

diff_4plot <- melt(diff_one_two[-c(13:18),])
diff_4plot$extract <-paste(rep("Extract",each=12),rep(LETTERS[1:6],each=12),sep=" ")

diff_4plot$Evaluation <- factor(diff_4plot$Evaluation,
                                                 levels = c("MBR_pooled","MBR_pooled\noldversion","MBR","noMBR_pooled","noMBR_pooled\noldversion","noMBR"))

plot_peptides <- ggplot(diff_4plot, aes(x=Type, y=value, color=Evaluation)) +
  geom_point(position = dodge, size = 2.5, alpha = 0.8) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Deviation one vs. two peptides [%]") +
  labs(fill = "Evaluation \n method") +
  facet_wrap(~extract) + 
  scale_color_manual(values = color_manual) + 
  guides(colour = guide_legend(ncol=6))
  
plot_peptides

pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_comparison_onetwopeptides_noproteins.png")

ggsave(pth,
       width = 7,
       height = 7) 



one <- readRDS("results_onepeptide.rds")
two <- readRDS("results_twopeptides.RDS")

one <- as.data.frame(one)
two <- as.data.frame(two)

diff_one_two <- ((one/two)-1)*100
diff_one_two <- diff_one_two[,grepl("rsd",colnames(diff_one_two))]
diff_one_two$Type <- one$X_A.Type
diff_one_two$Evaluation <- one$X_A.Evaluation

diff_4plot <- melt(diff_one_two[-c(13:18),])
diff_4plot$extract <-paste(rep("Extract",each=12),rep(LETTERS[1:6],each=12),sep=" ")

diff_4plot$Evaluation <- factor(diff_4plot$Evaluation,
                                levels = c("MBR_pooled","MBR_pooled\noldversion","MBR","noMBR_pooled","noMBR_pooled\noldversion","noMBR"))

plot_peptides <- ggplot(diff_4plot, aes(x=Type, y=value, color=Evaluation)) +
  geom_point(position = dodge, size = 2.5, alpha = 0.8) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Deviation one vs. two peptides [%]") +
  labs(fill = "Evaluation \n method") +
  facet_wrap(~extract) + 
  scale_color_manual(values = color_manual) +
  guides(colour = guide_legend(ncol=6))

plot_peptides

pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_comparison_onetwopeptides_rsd.png")
ggsave(pth,
       width = 7,
       height = 7) 

###########################################
#
# visualize difference secondpeptides
#
###########################################

# number of proteins

one_raw <- as.data.frame(readRDS("results_secondpeptides_onepeptide.rds"))
two_raw <- as.data.frame(readRDS("results_seconpeptides_twopeptides.rds"))

one <- one_raw
two <- two_raw

one <- one[,grepl("no",colnames(one))]
two <- two[,grepl("no",colnames(two))]

one_diff <- t(data.frame(((one[2,]/one[1,])-1)*100, # 6x mbr 6x razor 6xone
                         ((one[4,]/one[3,])-1)*100,            # 6x no MBR 6x razor 5x one
                         ((one[6,]/one[5,])-1)*100,            # 6x  MBR 6x unique 5x one
                         ((one[8,]/one[7,])-1)*100,            # 5x  no MBR 6x unique 5x one
                         ((one[10,]/one[9,])-1)*100,            # 5x MBR 6x lfq 5x one
                         ((one[12,]/one[11,])-1)*100))            # 5x no MBR 6x lfq 5x one

two_diff <-  t(data.frame(
  (((two[2,]/two[1,])-1)*100),
  (((two[4,]/two[3,])-1)*100),
  (((two[6,]/two[5,])-1)*100),
  (((two[8,]/two[7,])-1)*100),
  (((two[10,]/two[9,])-1)*100),
  (((two[12,]/two[11,])-1)*100)))

compare_secondpeptides <- as.data.frame(rbind(one_diff,two_diff))
colnames(compare_secondpeptides) <-c("results")
compare_secondpeptides$number <- rep(c("one","two"),each = 36)
compare_secondpeptides$MBR <- rep(rep(c("MBR","no_MBR"),each = 6),6)
compare_secondpeptides$Type <- rep(rep(c("Razor","Unique","LFQ"), each = 12),2)
compare_secondpeptides$Concentration <- as.factor(rep(c("A","B","C","D","E","F"), 12))
compare_secondpeptides$shapes <- paste(compare_secondpeptides$MBR,compare_secondpeptides$number
                                       ,sep="_")

compare_secondpeptides$Type <- factor(compare_secondpeptides$Type,
                                      levels = c("Razor", "Unique", "LFQ"))

dodge <- position_dodge(.3)

ggplot(compare_secondpeptides, aes(x=Concentration, y=results, color=shapes)) +
  geom_point(position = dodge, size = 3, aes(shape = shapes), alpha = 0.8) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(shape=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Deviation using second pepitdes [%]") +
  xlab("Extract") +
  labs(fill = "Evaluation \n method") +
  facet_wrap(~Type, ncol = 1) +
  scale_color_manual(name = "Evaluation method", values = c(mypalette_red[1],mypalette_red[1],mypalette_blue[1],mypalette_blue[1]), labels = c("MBR\none peptide ", "MBR\ntwo peptides", "no MBR\none peptide ", "no MBR\ntwo peptides")) + 
  scale_shape_manual(name = "Evaluation method", values=c(1, 16, 2, 17), labels = c("MBR\none peptide ", "MBR\ntwo peptides", "no MBR\none peptide ", "no MBR\ntwo peptides") )

pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_comparison_secondpeptides.png"
ggsave(pth,
       width = 7,
       height = 7)   

# rsd

one_raw <- as.data.frame(readRDS("results_secondpeptides_onepeptide.rds"))
two_raw <- as.data.frame(readRDS("results_seconpeptides_twopeptides.rds"))

one <- one_raw
two <- two_raw

one <- one[,grepl("rsd",colnames(one))]
two <- two[,grepl("rsd",colnames(two))]

one_diff <- t(data.frame(((one[2,]/one[1,])-1)*100, # 6x mbr 6x razor 6xone
                         ((one[4,]/one[3,])-1)*100,            # 6x no MBR 6x razor 5x one
                         ((one[6,]/one[5,])-1)*100,            # 6x  MBR 6x unique 5x one
                         ((one[8,]/one[7,])-1)*100,            # 5x  no MBR 6x unique 5x one
                         ((one[10,]/one[9,])-1)*100,            # 5x MBR 6x lfq 5x one
                         ((one[12,]/one[11,])-1)*100))            # 5x no MBR 6x lfq 5x one

two_diff <-  t(data.frame(
  (((two[2,]/two[1,])-1)*100),
  (((two[4,]/two[3,])-1)*100),
  (((two[6,]/two[5,])-1)*100),
  (((two[8,]/two[7,])-1)*100),
  (((two[10,]/two[9,])-1)*100),
  (((two[12,]/two[11,])-1)*100)))

compare_secondpeptides <- as.data.frame(rbind(one_diff,two_diff))
colnames(compare_secondpeptides) <-c("results")
compare_secondpeptides$number <- rep(c("one","two"),each = 36)
compare_secondpeptides$MBR <- rep(rep(c("MBR","no_MBR"),each = 6),6)
compare_secondpeptides$Type <- rep(rep(c("Razor","Unique","LFQ"), each = 12),2)
compare_secondpeptides$Concentration <- as.factor(rep(c("A","B","C","D","E","F"), 12))
compare_secondpeptides$shapes <- paste(compare_secondpeptides$MBR,compare_secondpeptides$number
                                       ,sep="_")

compare_secondpeptides$Type <- factor(compare_secondpeptides$Type,
                                      levels = c("Razor", "Unique", "LFQ"))

dodge <- position_dodge(.3)

ggplot(compare_secondpeptides, aes(x=Concentration, y=results, color=shapes)) +
  geom_point(position = dodge, size = 3, aes(shape = shapes), alpha = 0.8) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Deviation using second pepitdes [%]") +
  xlab("Extract") +
  labs(fill = "Evaluation \n method") +
  facet_wrap(~Type, ncol = 1) +
  scale_color_manual(name = "Evaluation method", values = c(mypalette_red[1],mypalette_red[1],mypalette_blue[1],mypalette_blue[1]), labels = c("MBR\none peptide ", "MBR\ntwo peptides", "no MBR\none peptide ", "no MBR\ntwo peptides")) + 
  scale_shape_manual(name = "Evaluation method", values=c(1, 16, 2, 17), labels = c("MBR\none peptide ", "MBR\ntwo peptides", "no MBR\none peptide ", "no MBR\ntwo peptides") )


pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_comparison_secondpeptides_rsd.png"
ggsave(pth,
       width = 7,
       height = 7)   

##########################################
#
# ups1
#
##########################################

library(ggpubr)
library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape)

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
# load plots

plots_no_twopeptides <- readRDS("plots_no_twopeptides.rds")
plots_no_onepeptide <- readRDS("plots_no_onepeptide.rds")
sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")


# create single plots for each concentration


for (sample_name in sample_names){
  ggarrange(plots_no_onepeptide[[sample_name]],plots_no_twopeptides[[sample_name]], 
            ncol = 1, nrow = 2,
            common.legend = TRUE, legend="bottom")
  pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/UPS1_rerun_",sample_name,".png",sep="")
  ggsave(pth,
         width = 7,
         height = 7)  
}



####### comparing via concentration two peptides

res_twopeptides <- readRDS("results_twopeptides.rds")
sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")
number_proteins <- c()
name <- c()
conc <- c()
rsd <- c()
type <- c()
evaluation <- c()

for (sample_name in sample_names){
  tmp <- res_twopeptides[[sample_name]]
  number_proteins <- append(number_proteins,tmp[,1])
  name <- append(name,
                 paste(tmp[,3],"_",tmp[,4], sep = ""))
  conc <- c(conc, rep(sample_name,6))
  type <- c(type, paste(tmp[,3]))
  evaluation <- c(evaluation, paste(tmp[,4]))
  rsd <- c(rsd, paste(tmp[,2]))
}


res_2gether <- as.data.frame(cbind(number_proteins,conc,name,evaluation, type, rsd))
colnames(res_2gether) <- c("NumberProteins","Concentration","Name","Evaluation","Type","rsd")
res_2gether$NumberProteins <- as.numeric(res_2gether$NumberProteins)
res_2gether$rsd <- as.numeric(res_2gether$rsd)
res_2gether$Concentration <- rep(c(10, 25, 2, 4, 50), each = 18)

plot_no_two <- ggplot(data=res_2gether, aes(x=Concentration, y=NumberProteins, group=Name, color = Evaluation, linetype = Type )) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("Number of proteins") + 
  xlab("Concentration [fmol]") +
  scale_color_manual(values = c(color_manual)) +
  labs(color = "Evaluation \n method", title = "two peptides - number proteins") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  guides(colour = guide_legend(ncol=6))

plot_rsd_two <- ggplot(data=res_2gether, aes(x=Concentration, y=rsd, group=Name, color = Evaluation, linetype = Type )) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("rsd [%]") + 
  xlab("Concentration [fmol]") +
  labs(color = "Evaluation \n method", title = "two peptides - rsd [%]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  scale_color_manual(values = c(color_manual)) +
  guides(colour = guide_legend(ncol=6))

ggarrange(plot_no_two, plot_rsd_two, 
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_noproteins_rsd_twopeptides.png"
ggsave(pth,
      width = 7,
      height = 7)  

res_2gether_twopeptides <- res_2gether

############################
####### comparing via concentration one peptide
#############################
res_peptides <- readRDS("results_onepeptide.rds")
sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")
number_proteins <- c()
name <- c()
conc <- c()
rsd <- c()
type <- c()
evaluation <- c()

for (sample_name in sample_names){
  tmp <- res_peptides[[sample_name]]
  number_proteins <- append(number_proteins,tmp[,1])
  name <- append(name,
                 paste(tmp[,3],"_",tmp[,4], sep = ""))
  conc <- c(conc, rep(sample_name,6))
  type <- c(type, paste(tmp[,3]))
  evaluation <- c(evaluation, paste(tmp[,4]))
  rsd <- c(rsd, paste(tmp[,2]))
}


res_2gether <- as.data.frame(cbind(number_proteins,conc,name,evaluation, type, rsd))
colnames(res_2gether) <- c("NumberProteins","Concentration","Name","Evaluation","Type","rsd")
res_2gether$NumberProteins <- as.numeric(res_2gether$NumberProteins)
res_2gether$rsd <- as.numeric(res_2gether$rsd)
res_2gether$Concentration <- rep(c(10, 25, 2, 4, 50), each = 18)
res_2gether$Evaluation <- as.character(res_2gether$Evaluation)

plot_no_one <- ggplot(data=res_2gether, aes(x=Concentration, y=NumberProteins, group=Name, color = Evaluation, linetype = Type )) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("Number of proteins") + 
  xlab("Concentration [fmol]") +
  scale_color_manual(values = c(color_manual)) +
  labs(color = "Evaluation \n method", title = "one peptide - number proteins") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) + 
  scale_linetype_manual(values=c("solid","dotted","longdash")) + 
  guides(colour = guide_legend(ncol=6))


plot_rsd_one <- ggplot(data=res_2gether, aes(x=Concentration, y=rsd, group=Name, color = Evaluation, linetype = Type )) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("rsd [%]") + 
  xlab("Concentration [fmol]") +
  labs(color = "Evaluation \n method", title = "one peptide - rsd [%]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  scale_color_manual(values = c(color_manual)) +
  guides(colour = guide_legend(ncol=6))


ggarrange(plot_no_one, plot_rsd_one, 
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_noproteins_rsd_onepeptide.png"
ggsave(pth,
       width = 7,
       height = 7)  


res_2gether_onepeptide <- res_2gether

# no and rsd one and two peptides plot

ggarrange(plot_no_one, plot_rsd_one, plot_no_two, plot_rsd_two,
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_noproteins_rsd.png"
ggsave(pth,
       width = 7,
       height = 7)  

#############################################
#
# compare one and two peptides
#
#############################################

# comparison one and two peptides
comp_one_two <- as.data.frame(((res_2gether_onepeptide$NumberProteins/res_2gether_twopeptides$NumberProteins)-1)*100)
comp_one_two$rsd <- ((res_2gether_onepeptide$rsd/res_2gether_twopeptides$rsd)-1)*100
colnames(comp_one_two) <- c("No","rsd [%]")
comp_one_two$Evaluation <- res_2gether_onepeptide$Evaluation
comp_one_two$Concentration <- res_2gether_onepeptide$Concentration
comp_one_two$conc_4plot <- paste(comp_one_two$Concentration,"fmol",sep="")
comp_one_two$conc_4plot <- factor(comp_one_two$conc_4plot,
                                  levels=c("2fmol","4fmol","10fmol","25fmol","50fmol"))
# melt results for ggplot
comp_one_two_4plot <- melt(comp_one_two, id = c("Evaluation", "Concentration","conc_4plot"))
# remove all 0 values
comp_one_two_4plot[comp_one_two_4plot == 0] <- NA


ggplot(comp_one_two_4plot, aes(x=variable, y=value)) + 
  geom_boxplot() +
  geom_point(aes(color = Evaluation)) +
  scale_color_manual(values = rep(c(mypalette_red,mypalette_blue),4)) +
  scale_alpha(0.8) +
  xlab("") +
  ylab("Deviation between one and two peptides [%]") +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(colour = guide_legend(ncol=2)) +
  facet_wrap(~conc_4plot, ncol = 5)

pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_compare_onetwopeptides.png"
ggsave(pth,
       width = 7,
       height = 7)  

###############################################
#
# compare no human proteins 
#
###############################################
# should be 47

path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/UPS1/tmpfiles/UPS1_oneuniquepeptide_human/"

##########read in all data with pattern .txt
# get file list
file.list <- list.files(path = path, pattern='*.csv')
# create file path
file.path <- paste(path,file.list,sep = "")
# read in files - for each file one variable is created - results are stored in a list
res_raw <- lapply(file.path, function(i){read.csv2(i, row.names = 1)})
names(res_raw) <- file.list

i <- 0
res_number_human_proteins_onepeptide <- c()
for (i in 1:4){
  res_number_human_proteins_onepeptide <- c(res_number_human_proteins_onepeptide,colMeans(res_raw[[i]]))
}


##### two peptides

path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/UPS1/tmpfiles/UPS1_twouniquepeptides_human/"

##########read in all data with pattern .txt
# get file list
file.list <- list.files(path = path, pattern='*.csv')
# create file path
file.path <- paste(path,file.list,sep = "")
# read in files - for each file one variable is created - results are stored in a list
res_raw <- lapply(file.path, function(i){read.csv2(i, row.names = 1)})
names(res_raw) <- file.list

i <- 0
res_number_human_proteins_twopeptides <- c()
for (i in 1:4){
  res_number_human_proteins_twopeptides <- c(res_number_human_proteins_twopeptides,colMeans(res_raw[[i]]))
}

type <- c(rep("Razor", 5),
          rep("Unique", 5),
          rep("LFQ", 5))
type <- rep(type,  4)
type <- rep(type, 2)
conc <- rep(c(10,25,2,4,50),12)
conc <- rep(conc, 2)
peptides <- c(rep("one", 60),
              rep("two",60))
evaluationmethod <- c(rep("MBR_pooled",15),
                      rep("MBR_pooled_oldversion",15),
                      rep("no_MBR_pooled",15),
                      rep("no_MBR_pooled_oldversion",15))
evaluationmethod <- rep(evaluationmethod, 2)

mean_human_proteins <- c(res_number_human_proteins_onepeptide,res_number_human_proteins_twopeptides)


results_numberhumanproteins <- data.frame(mean_human_proteins, type, conc, peptides, evaluationmethod)

evaluationmethod_forloop <- c("MBR_pooled", "MBR_pooled_oldversion", "no_MBR_pooled", "no_MBR_pooled_oldversion")
plots_numberhumanproteins <- list()

colors_human <- c(mypalette_red[c(1:2)],mypalette_blue[c(1:2)])
i = 0

for (i in 1:4){
  plot <- results_numberhumanproteins %>% 
    dplyr::filter(evaluationmethod == evaluationmethod_forloop[i]) %>%
    ggplot(., aes(x=conc, y=mean_human_proteins, color = type, linetype = peptides)) +
    geom_line() +
    geom_point() + 
    ylab("Number of proteins") + 
    xlab("Concentration [fmol]") +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
    labs(color = "Type", linetype = "Peptides") +
    scale_color_manual(values = rep(colors_human[i],3) ) +
    scale_alpha(0.8) +
    facet_wrap(~type)+
    ylim(c(30,47)) +
    ggtitle(gsub("_"," ",evaluationmethod_forloop[i]))
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  plots_numberhumanproteins[[i]] <- plot  
}

ggarrange(plotlist = plots_numberhumanproteins,
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_nohumanproteins.png"
ggsave(pth,
       width = 7,
       height = 7) 

#####################################################
#
# comparing MBR no MBR
#
#####################################################

# one peptide

MBR_comparison <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^noMBR"))
summary(MBR_comparison$rsd)
summary(noMBR$rsd)
MBR_comparison$comparison <- ((MBR_comparison$NumberProteins/noMBR$NumberProteins)-1)*100

plot_compMBR_one <- ggplot(data=MBR_comparison, aes(x=Concentration, y=comparison, group=Name, color = Evaluation, linetype = Type )) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("Proportion of more proteins\nfound using MBR [%]") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_discrete(labels = c("manually\npooled","pooled\nMaxQuant","pooled\nMaxQuant oldversion")) + 
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  scale_color_manual(values = mypalette_red) + 
  ggtitle("one peptide")
  
  

# two peptides

MBR_comparison <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^noMBR"))
summary(MBR_comparison$rsd)
summary(noMBR$rsd)
MBR_comparison$comparison <- ((MBR_comparison$NumberProteins/noMBR$NumberProteins)-1)*100

plot_compMBR_two <- ggplot(data=MBR_comparison, aes(x=Concentration, y=comparison, group=Name, color = Evaluation, linetype = Type)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.8) + 
  ylab("Proportion of more proteins\nfound using MBR [%]") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_discrete(labels = c("manually\npooled","pooled\nMaxQuant","pooled\nMaxQuant oldversion")) + 
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  scale_color_manual(values = mypalette_red) +
  ggtitle("two peptides")

ggarrange(plot_compMBR_one, plot_compMBR_two, 
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_comparisonMBRnoMBR.png"
ggsave(pth,
       width = 7,
       height = 7)  



#####################################################
#
# second peptides
#
#####################################################

# viszulaize difference number of proteins

one_raw <- as.data.frame(readRDS("results_secondpeptides_onepeptide.rds"))
two_raw <- as.data.frame(readRDS("results_secondpeptides_twopeptides.rds"))

one <- one_raw
two <- two_raw

one <- one[,grepl("no",colnames(one))]
two <- two[,grepl("no",colnames(two))]

one_diff <- t(data.frame(((one[2,]/one[1,])-1)*100, # 5x mbr 5x razor 5xone
              ((one[4,]/one[3,])-1)*100,            # 5x no MBR 5x razor 5x one
              ((one[6,]/one[5,])-1)*100,            # 5x  MBR 5x unique 5x one
              ((one[8,]/one[7,])-1)*100,            # 5x  no MBR 5x unique 5x one
              ((one[10,]/one[9,])-1)*100,            # 5x MBR 5x lfq 5x one
              ((one[12,]/one[11,])-1)*100))            # 5x no MBR 5x lfq 5x one

two_diff <-  t(data.frame(
               (((two[2,]/two[1,])-1)*100),
               (((two[4,]/two[3,])-1)*100),
               (((two[6,]/two[5,])-1)*100),
               (((two[8,]/two[7,])-1)*100),
               (((two[10,]/two[9,])-1)*100),
               (((two[12,]/two[11,])-1)*100)))

compare_secondpeptides <- as.data.frame(rbind(one_diff,two_diff))
colnames(compare_secondpeptides) <-c("results")
compare_secondpeptides$number <- rep(c("one","two"),each = 30)
compare_secondpeptides$MBR <- rep(rep(c("MBR","no_MBR"),each = 5),6)
compare_secondpeptides$Type <- rep(rep(c("Razor","Unique","LFQ"), each = 10),2)
compare_secondpeptides$Concentration <- as.factor(rep(c(10,25,2,4,50), 12))
compare_secondpeptides$shapes <- paste(compare_secondpeptides$MBR,compare_secondpeptides$number
                                        ,sep="_")

compare_secondpeptides$Type <- factor(compare_secondpeptides$Type,
                                      levels = c("Razor", "Unique", "LFQ"))

dodge <- position_dodge(.3)

ggplot(compare_secondpeptides, aes(x=Concentration, y=results, color=shapes)) +
  geom_point(position = dodge, size = 3, aes(shape = shapes), alpha = 0.8) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  ylab("Deviation using option second pepitdes [%]") +
  xlab("Concentration[fmol]") +
  labs(fill = "Evaluation \n method") +
  facet_wrap(~Type, ncol = 1) +
  scale_color_manual(name = "Evaluation method", values = c(mypalette_red[1],mypalette_red[1],mypalette_blue[1],mypalette_blue[1]), labels = c("MBR - one peptide ", "MBR - two peptides", "no MBR - one peptide ", "no MBR - two peptides")) + 
  scale_shape_manual(name = "Evaluation method", values=c(1, 16, 2, 17), labels = c("MBR - one peptide ", "MBR - two peptides", "no MBR - one peptide ", "no MBR - two peptides") )

 
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_comparison_secondpeptides.png"
ggsave(pth,
       width = 7,
       height = 7)   

# visualize differences using rsd

# viszulaize difference number of proteins

one_raw <- as.data.frame(readRDS("results_secondpeptides_onepeptide.rds"))
two_raw <- as.data.frame(readRDS("results_secondpeptides_twopeptides.rds"))

one <- one_raw
two <- two_raw

one <- one[,grepl("rsd",colnames(one))]
two <- two[,grepl("rsd",colnames(two))]

one_diff <- t(data.frame(((one[2,]/one[1,])-1)*100, # 5x mbr 5x razor 5xone
                         ((one[4,]/one[3,])-1)*100,            # 5x no MBR 5x razor 5x one
                         ((one[6,]/one[5,])-1)*100,            # 5x  MBR 5x unique 5x one
                         ((one[8,]/one[7,])-1)*100,            # 5x  no MBR 5x unique 5x one
                         ((one[10,]/one[9,])-1)*100,            # 5x MBR 5x lfq 5x one
                         ((one[12,]/one[11,])-1)*100))            # 5x no MBR 5x lfq 5x one

two_diff <-  t(data.frame(
  (((two[2,]/two[1,])-1)*100),
  (((two[4,]/two[3,])-1)*100),
  (((two[6,]/two[5,])-1)*100),
  (((two[8,]/two[7,])-1)*100),
  (((two[10,]/two[9,])-1)*100),
  (((two[12,]/two[11,])-1)*100)))

compare_secondpeptides <- as.data.frame(rbind(one_diff,two_diff))
colnames(compare_secondpeptides) <-c("results")
compare_secondpeptides$number <- rep(c("one","two"),each = 30)
compare_secondpeptides$MBR <- rep(rep(c("MBR","no_MBR"),each = 5),6)
compare_secondpeptides$Type <- rep(rep(c("Razor","Unique","LFQ"), each = 10),2)
compare_secondpeptides$Concentration <- as.factor(rep(c(10,25,2,4,50), 12))
compare_secondpeptides$shapes <- paste(compare_secondpeptides$MBR,compare_secondpeptides$number
                                       ,sep="_")

compare_secondpeptides$Type <- factor(compare_secondpeptides$Type,
                                      levels = c("Razor", "Unique", "LFQ"))

dodge <- position_dodge(.3)

ggplot(compare_secondpeptides, aes(x=Concentration, y=results, color=shapes)) +
  geom_point(position = dodge, size = 3, aes(shape = shapes), alpha = 0.8) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  ylab("Deviation using option second pepitdes [%]") +
  xlab("Concentration[fmol]") +
  labs(fill = "Evaluation \n method") +
  facet_wrap(~Type, ncol = 1) +
  scale_color_manual(name = "Evaluation method", values = c(mypalette_red[1],mypalette_red[1],mypalette_blue[1],mypalette_blue[1]), labels = c("MBR - one peptide ", "MBR - two peptides", "no MBR - one peptide ", "no MBR - two peptides")) + 
  scale_shape_manual(name = "Evaluation method", values=c(1, 16, 2, 17), labels = c("MBR - one peptide ", "MBR - two peptides", "no MBR - one peptide ", "no MBR - two peptides") )


pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_comparison_secondpeptides_rsd.png"
ggsave(pth,
       width = 7,
       height = 7)   


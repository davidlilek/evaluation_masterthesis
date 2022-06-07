##########################################
#
# ups1
#
##########################################

library(ggpubr)

plots_no_twopeptides <- readRDS("plots_no_twopeptides.rds")
plots_rsd_twopeptides <- readRDS("plots_rsd_twopeptides.rds")
plots_no_onepeptide <- readRDS("plots_no_onepeptide.rds")
plots_rsd_onepeptide <- readRDS("plots_rsd_onepeptide.rds")
sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")




for (sample_name in sample_names){
  ggarrange(plots_no_onepeptide[[sample_name]],plots_no_twopeptides[[sample_name]], plots_rsd_onepeptide[[sample_name]] ,plots_rsd_twopeptides[[sample_name]], 
            labels = c("A", "B", "C", "D"),
            ncol = 2, nrow = 2,
            common.legend = TRUE, legend="bottom")
  pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/UPS1_rerun_",sample_name,".png",sep="")
  ggsave(pth,
         width = 7,
         height = 7)  
}


####### comparing via concentration
res_twopeptides <- readRDS("results_twopeptides.rds")
#change also file name for png file if you run it with one peptide
#res_twopeptides <- readRDS("results_onepeptide.rds")
sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")
number_proteins <- c()
name <- c()
conc <- c()
rsd <- c()
type <- c()
evaluation <- c()
for (sample_name in sample_names){
  tmp <- res_twopeptides[[sample_name]]
  number_proteins <- append(number_proteins,tmp[c(1,4,7,10,13,16),1])
  name <- append(name,
                 paste(tmp[c(1,4,7,10,13,16),3],"_",tmp[c(1,4,7,10,13,16),4], sep = ""))
  conc <- c(conc, rep(sample_name,6))
  type <- c(type, paste(tmp[c(1,4,7,10,13,16),3]))
  evaluation <- c(evaluation, paste(tmp[c(1,4,7,10,13,16),4]))
  rsd <- c(rsd, paste(tmp[c(1,4,7,10,13,16),2]))
}


res_2gether <- as.data.frame(cbind(number_proteins,conc,name,evaluation, type, rsd))
colnames(res_2gether) <- c("NumberProteins","Concentration","Name","Evaluation","Type","rsd")
res_2gether$NumberProteins <- as.numeric(res_2gether$NumberProteins)
res_2gether$rsd <- as.numeric(res_2gether$rsd)
res_2gether$Concentration <- rep(c(10, 25, 2, 4, 50), each = 6)

ggplot(data=res_2gether, aes(x=Concentration, y=NumberProteins, group=Name, color = Type, linetype = Evaluation )) +
  geom_line() +
  geom_point() + 
  ylab("Number of Proteins") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())


ggplot(data=res_2gether, aes(x=Concentration, y=rsd, group=Name, color = Type, linetype = Evaluation )) +
  geom_line() +
  geom_point() + 
  ylab("Number of Proteins") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

########with all values
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

plot_no <- ggplot(data=res_2gether, aes(x=Concentration, y=NumberProteins, group=Name, color = Evaluation, linetype = Type )) +
  geom_line() +
  geom_point() + 
  ylab("Number of Proteins") + 
  xlab("Concentration [fmol]") +
  scale_fill_brewer(palette="Dark2") +
  labs(color = "Evaluation \n method") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())


plot_rsd <- ggplot(data=res_2gether, aes(x=Concentration, y=rsd, group=Name, color = Evaluation, linetype = Type )) +
  geom_line() +
  geom_point() + 
  ylab("rsd [%]") + 
  xlab("Concentration [fmol]") +
  labs(color = "Evaluation \n method") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())



ggarrange(plot_no, plot_rsd, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_noproteins_rsd.png"
ggsave(pth,
      width = 7,
      height = 7)  


######################## compare no human proteins 
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
                      rep("noMBR_pooled",15),
                      rep("noMBR_pooled_oldversion",15))
evaluationmethod <- rep(evaluationmethod, 2)

mean_human_proteins <- c(res_number_human_proteins_onepeptide,res_number_human_proteins_twopeptides)


results_numberhumanproteins <- data.frame(mean_human_proteins, type, conc, peptides, evaluationmethod)

evaluationmethod_forloop <- c("MBR_pooled", "MBR_pooled_oldversion", "noMBR_pooled", "noMBR_pooled_oldversion")
plots_numberhumanproteins <- list()
for (evaluation in evaluationmethod_forloop){
  plot <- results_numberhumanproteins %>% 
    dplyr::filter(evaluationmethod == evaluation) %>%
    ggplot(., aes(x=conc, y=mean_human_proteins, color = type, linetype = peptides)) +
    geom_line() +
    geom_point() + 
    ylab("Number of Proteins") + 
    xlab("Concentration [fmol]") +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())+
    labs(color = "Type", linetype = "Peptides") +
    facet_wrap(~type)+
    ylim(c(30,47)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  plots_numberhumanproteins[[evaluation]] <- plot  
}

ggarrange(plotlist = plots_numberhumanproteins,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_nohumanproteins.png"
ggsave(pth,
       width = 7,
       height = 7) 


  
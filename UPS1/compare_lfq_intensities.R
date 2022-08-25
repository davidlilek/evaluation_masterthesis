######################################
# compare lfq intensities
#
# reads in log2 lfq intensities from rerun experiment @ different concentration
# evaluation was performed with two unique peptides 
# for each quantified protein a linear model was fitted (conc was log2 transformed before)
# r2, intercept and slope were saved an visualized
# res_raw[[3]][[8]] was removed because only 46 human proteins occured
#
# A MBR_pooled
# B MBR_pooled_oldversion
# C noMBR_pooled
# D noMBR_pooled_oldversion
#
#######################################

# load libraries
library(dplyr)
library(ggplot2)
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


sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")
onepeptide <- list()
twopeptides <- list()

########################### one peptide

path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/UPS1/tmpfiles/UPS1_oneuniquepeptide_human_lfq_comparison/"

##########read in all data with pattern .txt
# get file list
file.list <- list.files(path = path, pattern='*.rds')
# remove second peptides
file.list <- file.list[-c(3,6)]
# create file path
file.path <- paste(path,file.list,sep = "")
# read in files - for each file one variable is created - results are stored in a list
res_raw <- lapply(file.path, function(i){readRDS(i)})
names(res_raw) <- file.list

#remove run 1 because it has only 46 human proteins
res_raw[[1]][[1]] <- NA
res_raw[[3]][[8]] <- NA

i <- 0 
for (i in 1:4){
  results_lfq_human <- do.call(cbind, res_raw[[i]])
  results_lfq_human_mean <- c()
  for (sample_name in sample_names){
    tmp <- results_lfq_human %>%
      dplyr::select(contains(sample_name))
    results_lfq_human_mean <- c(results_lfq_human_mean, rowMeans(tmp))
    
  }
  conc <- rep(c(10,25,2,4,50),each=47)
  protein <- rep(1:47,5)
  results_lfq_human_mean <- data.frame(results_lfq_human_mean, conc, protein)
  colnames(results_lfq_human_mean) <- c("mean", "conc", "protein")
  
  res <- results_lfq_human_mean 
  res$protein <- as.factor(res$protein)
  fitted_models <- results_lfq_human_mean %>% 
    group_by(protein) %>% 
    filter(!any(is.na(mean))) %>%
    do(model = lm(mean ~ log2(conc), data = .))
  
  
  onepeptide[[i]] <- fitted_models %>% do(data.frame(coef_intercept = coef(.$mod)[[1]],
                                                     coef_slope = coef(.$mod)[[2]],
                                                      r2 = summary(.$mod)$r.squared,
                                                      group = .[[1]]))
  
  
}


# merge results of the linear models for each evaluation method
lfq_2gether <- merge(onepeptide[[1]], onepeptide[[2]], 
      all = TRUE,
      by = "group")

lfq_2gether <- merge(lfq_2gether, onepeptide[[3]], all = TRUE, by = "group")
lfq_2gether <- merge(lfq_2gether, onepeptide[[4]], all = TRUE, by = "group")
# get all lfq_2gether in a vector
lfq_2gether_vector <- unlist(lfq_2gether[,-1], use.names = FALSE)


# create variable names and value factors

name <- c(
  paste(rep("MBR\n",39*3),"pooled",sep=""),
  paste(rep("MBR\n",39*3),"pooled\noldversion",sep=""),
  paste(rep("no MBR\n",39*3),"pooled",sep=""),
  paste(rep("no MBR\n",39*3),"pooled\noldversion",sep="")
)

value <- c(rep("intercept",39),
  rep("slope", 39),
           rep("r2",39))
value <- rep(value,4)

res <- data.frame(lfq_2gether_vector, name, value)

ggplot(res, aes(x=name, y=lfq_2gether_vector, fill=name)) +
  geom_boxplot(alpha=0.8) +
  facet_wrap(.~value, scales = "free") + 
  labs(y="", x="Evaluation method") +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank()) +
  scale_fill_manual(values = c(mypalette_red[c(1:2)],mypalette_blue[c(1:2)])) + 
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())



pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_linearmodels_lfqintensity.png"
ggsave(pth,
        width = 10,
        height = 7)  

# additional code
MBR <- subset(res, name == "onepetide1" | name == "onepetide2" )
noMBR <- subset(res, name == "onepetide3" | name == "onepetide4" )

boxplot(MBR$lfq_2gether_vector,noMBR$lfq_2gether_vector)

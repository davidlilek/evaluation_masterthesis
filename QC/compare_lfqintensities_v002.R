############################################################
#
# script for comparing lfq values of different experiments
#
#
# file path are defined and lfq values of all experiments are extracted
# then each experiment is compared to each other
# threshold for NA values can be defined
# for all comparisons mean, sd, difference mean and a t.test is performed
#   t.test
#     if sd(A) or sd(B) = 0 -> check whether means are equal (p=0) or not (p=1)
# for each comparison volcano plot is created 
# all volcano plots are visualized together and saved as png file
# one stacked barplot for all comparisons is drawn
#
#
############################################################

#load libraries
library(tidyverse)
library(dplyr)
library(gplots)
library(ggfortify)
library(gtools)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)

# load functions
filtering <- function(data_raw){
  data_raw %>% filter( 
    Potential.contaminant == "",
    Reverse == "",
    Only.identified.by.site == "")
}

# get colors
mypalette_grey <-brewer.pal(9,"Greys")


####################################
#
# define file path and create grid
#
####################################

# define project path
# for comparison project path and folder name of result file will be combined 
# e.g. N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/3_20220406_TR/20220406_TR_QC_MBR_pooled_rerun_4masterthesis
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/3_20220406_TR/"

# define folder name of result file
files <- c("20220406_TR_QC_MBR_pooled_rerun_4masterthesis",
           "20220406_TR_QC_MBR_pooled_rerun_oldversion_4masterthesis",
           "20220406_TR_QC_noMBR_pooled_rerun_4masterthesis",
           "20220406_TC_QC_noMBR_pooled_rerun_oldversion_4masterthesis")

#create grid
grid <- expand.grid(files,files)
grid <- as.data.frame(grid)
grid$Var1 <- as.character(grid$Var1)
grid$Var2 <- as.character(grid$Var2)
# remove duplicates e.g. A-A or B-B
keep <- apply(grid, 1, function(x) length(unique(x[!is.na(x)])) != 1)
grid <- grid[keep, ]
# get only unique combinations
grid <- unique(t(apply(grid, 1, sort)))


#######################################
#
# make comparisons
#
#######################################

# create needed variables
res_comparison <- list()
plots <- list()

# set values to 0
grd = 0
i = 0

# loop over of all possible combinations
for (grd in 1:nrow(grid)){
  ######################
  # files A comparison
  ######################
  path_A = paste(pth, grid[grd,1],"/", sep="")
  # get all txt files in folder
  file.list <- list.files(path = path_A, pattern='*.txt')
  # create file path
  file.path <- paste(path_A,file.list,sep = "")
  # read in files - for each file one variable is created - results are stored in a list
  res_raw <- lapply(file.path, read.delim)
  names(res_raw) <- file.list
  
  d.lfq.2gether <- list()
  
  #loop over all files in folder and store in in a list of list d.lfq.2gether
  for (x in file.list){
    data_raw <- as.data.frame(res_raw[[x]])
    d.filtered <- filtering(data_raw)
    d.lfq <- d.filtered %>% dplyr::select(starts_with("LFQ"),starts_with("FASTA"))
    d.lfq$LFQ.intensity.QC[d.lfq$LFQ.intensity.QC == 0] <- NA
    d.lfq.2gether[[x]] <- d.lfq
  }
  
  # merge lfq values for run 1 to 10
  # therefore first run 1 and 2 have to be merged 
  # for run 3-10 a for loop can be used
  m_A <- merge(d.lfq.2gether[[1]],d.lfq.2gether[[2]], 
               all = TRUE,
               by = "Fasta.headers")
  
  for (i in 3:10){
    m_A <- merge(m_A, d.lfq.2gether[[i]],
                 all = TRUE,
                 by = "Fasta.headers")
  }
  
  # define colnames of merged A (m_A)
  colnames(m_A) <- c("FASTA",paste(rep("A_",10),1:10,sep=""))
  
  #######################
  # files b comparison
  #######################
  
  # same code as for files A comparison -> no comments in this section
  path_B = paste(pth, grid[grd,2] , "/", sep = "")
  file.list <- list.files(path = path_B, pattern='*.txt')
  file.path <- paste(path_B,file.list,sep = "")
  res_raw <- lapply(file.path, read.delim)
  names(res_raw) <- file.list
  d.lfq.2gether <- list()
  for (x in file.list){
    data_raw <- as.data.frame(res_raw[[x]])
    d.filtered <- filtering(data_raw)
    d.lfq <- d.filtered %>% dplyr::select(starts_with("LFQ"),starts_with("FASTA"))
    d.lfq$LFQ.intensity.QC[d.lfq$LFQ.intensity.QC == 0] <- NA
    d.lfq.2gether[[x]] <- d.lfq
  }
  
  m_B <- merge(d.lfq.2gether[[1]],d.lfq.2gether[[2]], 
               all = TRUE,
               by = "Fasta.headers")
  
  for (i in 3:10){
    m_B <- merge(m_B, d.lfq.2gether[[i]],
                 all = TRUE,
                 by = "Fasta.headers")
  }
  colnames(m_B) <- c("FASTA",paste(rep("B_",10),1:10,sep=""))
  
  m_2gether <- merge(m_A,m_B,
                     all = TRUE,
                     by = "FASTA")
  
  ###################################
  # prepare data for comparison A-B
  ###################################
  
  # define column index for A & B
  A <- 2:11
  B <- 12:21
  
  # create data.frame for mean/sd/p value
  result <- data.frame(matrix(NA,   
                              nrow = nrow(m_2gether),
                              ncol = 0))
  rownames(result) <- m_2gether$FASTA
  result[c("pvalue","mean_A","sd_A","mean_B","sd_B")] <- NA
  
  # define threshold for NA values of A & B
  b.A <- rowSums(is.na(m_2gether[,A])) <= 5
  b.B <- rowSums(is.na(m_2gether[,B])) <= 5
  # calculate row wise mean for A and B
  mean.A <- rowMeans(m_2gether[rowSums(is.na(m_2gether))<= 2,A],
                 na.rm = TRUE)
  mean.B <- rowMeans(m_2gether[rowSums(is.na(m_2gether))<= 2,B],
                     na.rm = TRUE)
  
  i=0
  
  # perform calculation
  # first check if A and B fulfill threshold criteria -> if OK perform mean and sd calculation
  # if not sd(A) or sd(B) = 0 perform t test else set p value to 0 (means different) or 1 (means equal)
  # if too many NA values for A or B set p value to 0
  for (i in seq(1:nrow(m_2gether))){
    if (b.A[i] && b.B[i]){
      result[i,"mean_A"] <- mean(as.numeric(m_2gether[i,A]), na.rm = TRUE)
      result[i,"sd_A"] <- sd(as.numeric(m_2gether[i,A]), na.rm = TRUE)
      result[i,"mean_B"] <- mean(as.numeric(m_2gether[i,B]), na.rm = TRUE)
      result[i,"sd_B"] <- sd(as.numeric(m_2gether[i,B]), na.rm = TRUE)
      
      if (!(result[i,"sd_A"] == 0 || result[i,"sd_B"] == 0)){
        ttest <- t.test(m_2gether[i,A],m_2gether[i,B])
        result[i,"pvalue"] <- ttest$p.value
      } else {
        if (result[i, "mean_A"] == result[i, "mean_B"]){
          result[i,"pvalue"] <- 1
        } else {result[i,"pvalue"] <- 0}
      }
      
    } else if (!(b.A[i] && b.B[i])){
      result[i,"pvalue"] <- NA
      result[i,"mean_A"] <- mean(as.numeric(m_2gether[i,A]), na.rm = TRUE)
      result[i,"sd_A"] <- sd(as.numeric(m_2gether[i,A]), na.rm = TRUE)
      result[i,"mean_B"] <- mean(as.numeric(m_2gether[i,B]), na.rm = TRUE)
      result[i,"sd_B"] <- sd(as.numeric(m_2gether[i,B]), na.rm = TRUE)
    }  
  }
  
  # remove all rows which only contain NA values
  result.clear <- result[rowSums(is.na(result[,-1])) != ncol(result[,-1]), ]
  # remove all results with p value NA
  result.clear <- result.clear[!is.na(result.clear$pvalue),]
  # p value adjustment
  p <- result.clear$pvalue
  result.clear$pvalue <- p.adjust(p, method = "BH", n = length(p))
  # set all NaN values to 0
  result.clear$mean_A[is.nan(result.clear$mean_A)] <- 0
  result.clear$mean_B[is.nan(result.clear$mean_B)] <- 0
  # calculate mean difference
  # therefore log2 of LFQ values is calculated
  tmp <- log2(result.clear$mean_A)
  tmp[tmp == -Inf] <- NA
  mean_A <- tmp
  tmp <- log2(result.clear$mean_B)
  tmp[tmp == -Inf] <- NA
  mean_B <- tmp
  
  diff_mean <- mean_A - mean_B
  diff_mean <- as.data.frame(diff_mean)
  diff_mean$treatment <- rep("A",nrow(diff_mean))
  
  # violin plot for mean difference
  p <- ggplot(diff_mean, aes(x=treatment, y=diff_mean)) + 
    geom_violin()
  p
  # get number of p values  below 0.05, equal to 1 and all
  length(result.clear$pvalue[result.clear$pvalue<0.05])
  length(result.clear$pvalue[result.clear$pvalue==1])
  length(result.clear$pvalue)
  
  
  # volcano plots 
  tmp <- foldchange(result.clear[,2],result.clear[,4])
  result.clear$log2FC <- foldchange2logratio(tmp)
  result.clear$col <- cut(result.clear$pvalue,
                          breaks = c(-Inf, 0.05, Inf),
                          labels = c("significant", "not significant"))
  
  # set column significant
  result.clear$significant <- cut(result.clear$pvalue,
                                  breaks = c(-Inf, 0.05, Inf),
                                  labels = c("significant", "not significant"))
  # calculate difference
  result.clear$diff <- result.clear$mean_A-result.clear$mean_B
  # set column for difference
  result.clear$difference <- cut(result.clear$diff,
                                 breaks = c(-Inf,0,Inf),
                                 labels = c("lower","higher"))
  levels(result.clear$difference) <- c(levels(result.clear$difference),"equal")
  result.clear$difference[which(result.clear$diff == 0)] <- "equal"
  
  result.clear$group <- paste(result.clear$significant,result.clear$difference)
  result.clear$group[result.clear$group == "not significant equal"] <- "equal"
  result.clear$group[result.clear$group == "significant equal"] <- "equal"
  
  result.clear$col <- as.character(result.clear$col)
  p <- ggplot(data = result.clear, aes(x = log2FC, y = -1*log10(pvalue), colour = col))+
    geom_point()
  
  # store volcano plots and results of comparison
  plots[[grd]] <- p
  res_comparison[[grd]] <- result.clear
  
}

# volcanoplot of all comparison
ggarrange(plotlist = plots, nrow = 3, ncol = 2,
                   labels = LETTERS[1:6],
          common.legend = TRUE, legend = "bottom")
ggsave("../pics/volcanoplots.png", width = 7, height = 7)

#########################################
# create stacked barplot for comparisons
#########################################

# result vector
res_plot <- c()

i = 0

# loop over all comparisons and 
for (i in 1:6){
  tmp <- as.data.frame(table(res_comparison[[i]]$group))
  # remove number 6 -> NA
  res_plot <- c(res_plot, tmp[-6,2])
}

# get names for data.frame
names <- tmp[,1]
# create group names
group <- rep(as.character(names[-6]), 6)
# create comparison letters
comparison <- rep(LETTERS[1:6], each = 5)
# create data.frame
res_plot_group <- data.frame(cbind(group, comparison, res_plot))
colnames(res_plot_group) <- c("Group", "Comparison", "counts")
res_plot_group$counts <- as.double(res_plot_group$counts)
# create stacked barplot

#reverse color to get the more intense for MBR
mypalette_red <- mypalette_grey
mypalette_red <- mypalette_red[c(9,7,5,3,2)]

res_plot_group$Group <- factor(res_plot_group$Group,
                  levels = c("equal",
                            "significant higher",
                            "significant lower",
                            "not significant higher",
                            "not significant lower"))
b_plot <- ggplot(res_plot_group, aes(fill = Group, y = counts, x=Comparison, label = counts)) +
  geom_bar( position = "stack",stat="identity", alpha = 0.8, color = "black") + 
  scale_fill_manual(values = mypalette_red) +
  scale_linetype_manual(values = c("dashed")) +
  theme_bw() + 
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.40, 'cm'))
  
b_plot
ggsave("../pics/barplot_lfqintensities.png",
       width = 7,
       height = 3.5)


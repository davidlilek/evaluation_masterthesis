library(dplyr)
library(ggplot2)

sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")
onepeptide <- list()
twopeptides <- list()

########################### one peptide

path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/UPS1/tmpfiles/UPS1_oneuniquepeptide_human_lfq_comparison/"

##########read in all data with pattern .txt
# get file list
file.list <- list.files(path = path, pattern='*.rds')
# create file path
file.path <- paste(path,file.list,sep = "")
# read in files - for each file one variable is created - results are stored in a list
res_raw <- lapply(file.path, function(i){readRDS(i)})
names(res_raw) <- file.list

#remove run 7 because it has only 46 human proteins
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
    do(model = lm(mean ~ conc, data = .))
  
  
  onepeptide[[i]] <- fitted_models %>% do(data.frame(coef = coef(.$mod)[[2]],
                                                      r2 = summary(.$mod)$r.squared,
                                                      group = .[[1]]))
  
  
}


########################### two peptides

path = "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/UPS1/tmpfiles/UPS1_twouniquepeptides_human_lfq_comparison/"

##########read in all data with pattern .txt
# get file list
file.list <- list.files(path = path, pattern='*.rds')
# create file path
file.path <- paste(path,file.list,sep = "")
# read in files - for each file one variable is created - results are stored in a list
res_raw <- lapply(file.path, function(i){readRDS(i)})
names(res_raw) <- file.list

#remove run 7 because it has only 46 human proteins
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
    do(model = lm(mean ~ conc, data = .))
  
  
  twopeptides[[i]] <- fitted_models %>% do(data.frame(coef = coef(.$mod)[[2]],
                                  r2 = summary(.$mod)$r.squared,
                                  group = .[[1]]))
  
  
}

######### create plots
lfq_2gether <- merge(onepeptide[[1]], onepeptide[[2]], 
      all = TRUE,
      by = "group")

lfq_2gether <- merge(lfq_2gether, onepeptide[[3]], all = TRUE, by = "group")
lfq_2gether <- merge(lfq_2gether, onepeptide[[4]], all = TRUE, by = "group")
lfq_2gether <- merge(lfq_2gether, twopeptides[[1]], all = TRUE, by = "group")
lfq_2gether <- merge(lfq_2gether, twopeptides[[2]], all = TRUE, by = "group")
lfq_2gether <- merge(lfq_2gether, twopeptides[[3]], all = TRUE, by = "group")
lfq_2gether <- merge(lfq_2gether, twopeptides[[4]], all = TRUE, by = "group")

lfq_2gether_vector <- unlist(lfq_2gether[,-1], use.names = FALSE)
lfq_2gether_vector

name <- c(
  paste(rep("onepetide",39*2),"1",sep=""),
  paste(rep("onepetide",39*2),"2",sep=""),
  paste(rep("onepetide",39*2),"3",sep=""),
  paste(rep("onepetide",39*2),"4",sep=""),
  paste(rep("twopetides",39*2),"1",sep=""),
  paste(rep("twopetides",39*2),"2",sep=""),
  paste(rep("twopetides",39*2),"3",sep=""),
  paste(rep("twopetides",39*2),"4",sep="")
)

value <- c(rep("slope", 39),
           rep("r2",39))
value <- rep(value,8)

res <- data.frame(lfq_2gether_vector, name, value)

ggplot(res, aes(x=name, y=lfq_2gether_vector, fill=value)) +
  geom_boxplot() +
  scale_y_continuous(
    "Slope", 
    sec.axis = sec_axis(~ . /20, name = "r2")
  )

ggplot() +
  geom_boxplot(aes(x=name, y=lfq_2gether_vector, fill=value), data = subset(res, value == "slope")) +
  geom_boxplot(aes(x=name, y=lfq_2gether_vector, fill=value), data = subset(res, value == "r2")) +
  scale_y_continuous("slope",
                     sec.axis = sec_axis(~ . /20, name = "r2")
  )


data_slope <- subset(res, value == "slope")
data_r2 <- subset(res, value == "r2")

data_combined <- cbind(data_slope,data_r2)
colnames(data_combined) <- c("slope","name","value","r2","name","value")
data_combined <- data_combined[,-c(2,3)]
ggplot() +
  geom_boxplot(aes(x=name, y=slope), data = data_combined) +
  geom_boxplot(aes(x=name, y=r2), data = data_combined) +
  scale_y_continuous("slope",
                     sec.axis = sec_axis(~ . /20, name = "r2")
  )

ggplot(res, aes(x=name, y=lfq_2gether_vector, fill=value)) +
  geom_boxplot() +
  facet_wrap(.~value, scales = "free") + 
  labs(y="", x="Evaluation method") +
  scale_x_discrete(labels= c("one A",
                               "one B",
                               "one C",
                               "one D",
                               "two A",
                               "two B",
                               "two C",
                               "two D")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


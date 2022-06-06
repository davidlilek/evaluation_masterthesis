test <- results_razor[[1]]
view(test)
test_filtered <- dplyr::filter(test, grepl("Human",test$FASTA))
sapply(dplyr::filter(test, grepl("Human",test$FASTA))[-6],function(x)sum(x, na.rm = TRUE))



lfq_filtered<-dplyr::filter(results_lfq[[1]],grepl("Human",results_lfq[[1]]$FASTA))
lfq_filtered__7<-dplyr::filter(results_lfq[[7]],grepl("Human",results_lfq[[7]]$FASTA))

lfq_filtered$FASTA == lfq_filtered__2$FASTA



results_lfq_human
colMeans(do.call(rbind, results_lfq_human),na.rm = TRUE)

lapply(results_lfq_human,mean)

apply(simplify2array(results_lfq_human), 1:2, mean)

sample_names <- c("10fmol","25fmol","2fmol","4fmol","50fmol")
results_lfq_human <- do.call(cbind, results_lfq_human)

##### here it starts evluation mean lfq intensities

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

ggplot(results_lfq_human_mean, aes(x=as.factor(protein), y=mean, fill=as.factor(conc))) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal() +
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("No. of protein identifications") +
  labs(fill = "Evaluation \n method") 

ggplot(data=results_lfq_human_mean, aes(x=conc, y=mean, color = protein)) +
  geom_line() +
  geom_point() + 
  ylab("Number of Proteins") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
 
res <- results_lfq_human_mean 
res$protein <- as.factor(res$protein)
fitted_models <- results_lfq_human_mean %>% 
  group_by(protein) %>% 
  filter(!any(is.na(mean))) %>%
  do(model = lm(mean ~ conc, data = .))

  
fitted_models %>% do(data.frame(coef = coef(.$mod)[[1]],
                                r2 = summary(.$mod)$r.squared,
                                group = .[[1]]))

library(dplyr)
library(tidyr)
library(broom)
fitted_models_coeff <- glance(fitted_models, model)

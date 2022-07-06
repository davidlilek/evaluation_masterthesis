##########################################
#
# extracts
#
##########################################

#(1) run 20220406_TR_extract_summary_rerun_two peptides
#(2) run 20220406_TR_extract_summary_rerun_two peptides

library(ggpubr)
library(stringr)
library(dplyr)
library(ggplot2)

plots_no_twopeptides <- readRDS("plots_no_twopeptides.rds")
plots_rsd_twopeptides <- readRDS("plots_rsd_twopeptides.rds")
plots_no_onepeptide <- readRDS("plots_no_onepeptide.rds")
plots_rsd_onepeptide <- readRDS("plots_rsd_onepeptide.rds")
sample_names <- c("_A","_B","_C","_D","_E","_F")




for (sample_name in sample_names){
  ggarrange(plots_no_onepeptide[[sample_name]],plots_no_twopeptides[[sample_name]], plots_rsd_onepeptide[[sample_name]] ,plots_rsd_twopeptides[[sample_name]], 
            labels = c("A", "B", "C", "D"),
            ncol = 2, nrow = 2,
            common.legend = TRUE, legend="bottom")
  pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun",sample_name,".png",sep="")
  ggsave(pth,
         width = 7,
         height = 7)  
}

######## for second peptide analysis

plots_no_twopeptides <- readRDS("plots_no_secondpeptides_twopeptides.rds")
plots_rsd_twopeptides <- readRDS("plots_rsd_secondpeptides_twopeptides.rds")
plots_no_onepeptide <- readRDS("plots_no_secondpeptides_onepeptide.rds")
plots_rsd_onepeptide <- readRDS("plots_rsd_secondpeptides_onepeptide.rds")
sample_names <- c("_A","_B","_C","_D","_E","_F")


for (sample_name in sample_names){
  ggarrange(plots_no_onepeptide[[sample_name]],plots_no_twopeptides[[sample_name]], plots_rsd_onepeptide[[sample_name]] ,plots_rsd_twopeptides[[sample_name]], 
            labels = c("A", "B", "C", "D"),
            ncol = 2, nrow = 2,
            common.legend = TRUE, legend="bottom")
  pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_secondpeptides",sample_name,".png",sep="")
  ggsave(pth,
         width = 7,
         height = 7)  
}




### combine A-B

sample_name1 <- "_C"
sample_name2 <- "_D"
ggarrange(plots_no_onepeptide[[sample_name1]],plots_no_twopeptides[[sample_name1]], plots_rsd_onepeptide[[sample_name1]] ,plots_rsd_twopeptides[[sample_name1]], 
         plots_no_onepeptide[[sample_name2]],plots_no_twopeptides[[sample_name2]], plots_rsd_onepeptide[[sample_name2]] ,plots_rsd_twopeptides[[sample_name2]], 
         labels = c("A", "B", "C", "D","a","b","c","d"),
         ncol = 2, nrow = 4,
         common.legend = TRUE, legend="top",
         hjust = c(rep(-4.5,4),rep(-5.25,4)),
         align = "hv")
pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_secondpeptides_compare",sample_name1,sample_name2,".png",sep="")
ggsave(pth,
       width = 7,
       height = 7)  



################# compare MBR no MBR

# one peptide
sample_names <- LETTERS[1:6]
res_2gether_onepeptide <- readRDS("results_onepeptide.rds")
res_2gether_onepeptide <- as.data.frame(do.call(rbind,res_2gether_onepeptide))
res_2gether_onepeptide$sample_name <- rep(sample_names,each=18)

MBR_comparison <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^noMBR"))
summary(MBR_comparison$rsd)
summary(noMBR$rsd)
MBR_comparison$comparison <- ((MBR_comparison$no/noMBR$no)-1)*100
MBR_comparison_onepeptide <- MBR_comparison

# two peptides
res_2gether_twopeptides <- readRDS("results_twopeptides.rds")
res_2gether_twopeptides <- as.data.frame(do.call(rbind,res_2gether_twopeptides))
res_2gether_twopeptides$sample_name <- rep(sample_names,each=18)

MBR_comparison <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^noMBR"))
summary(MBR_comparison$rsd)
summary(noMBR$rsd)
MBR_comparison$comparison <- ((MBR_comparison$no/noMBR$no)-1)*100
MBR_comparison_twopeptides <- MBR_comparison


# bring the 2gether
MBR_comparison <- rbind(MBR_comparison_onepeptide,MBR_comparison_twopeptides)
MBR_comparison$peptides <- rep(c("one","two"),each=54)


 ggplot(MBR_comparison, aes(x=Type, y=comparison, color=Evaluation)) +
  geom_boxplot() +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Proportion of more proteins\nfound using MBR [%]") +
  labs(fill = "Evaluation \n method") + 
  facet_wrap(~sample_name, ncol=2) +
  scale_fill_brewer(palette="Dark2", labels = c("pooled\nMaxQuant","pooled\nMaxQuant oldversion","manually\npooled","manually\npooled")) 

 
 p_plot_one
p_plot_one <- ggplot(MBR_comparison, aes(x=Type, y=comparison, fill=Evaluation)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  scale_fill_brewer(palette="Dark2", labels = c("pooled\nMaxQuant","pooled\nMaxQuant oldversion","manually\npooled","manually\npooled")) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Proportion of more proteins\nfound using MBR [%]") +
  labs(fill = "Evaluation \n method") 

p_plot_one





p_plot_two <- ggplot(MBR_comparison, aes(x=Type, y=comparison, fill=Evaluation)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  scale_fill_brewer(palette="Dark2", labels = c("pooled\nMaxQuant","pooled\nMaxQuant oldversion","manually\npooled","manually\npooled")) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.15, 'cm')) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  ylab("Proportion of more proteins\nfound using MBR [%]") +
  labs(fill = "Evaluation \n method")

p_plot_two


# mbr no mbr
ggarrange(
          labels = c("A", "B"),
          ncol = 2, nrow = 4,
          common.legend = TRUE, legend="top",
          hjust = c(rep(-4.5,4),rep(-5.25,4)),
          align = "hv")
pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun_secondpeptides_compare",sample_name1,sample_name2,".png",sep="")
ggsave(pth,
       width = 7,
       height = 7) 
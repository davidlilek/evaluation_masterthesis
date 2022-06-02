##########################################
#
# extracts
#
##########################################

#(1) run 20220406_TR_extract_summary_rerun_two peptides
#(2) run 20220406_TR_extract_summary_rerun_two peptides

library(ggpubr)

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





###############################################################
#
#  OLD STUFF
#
###################################

# for (sample_name in sample_names){
#   ggarrange(plots_no_twopeptides[[sample_name]], plots_rsd_twopeptides[[sample_name]], 
#             labels = c("A", "B"),
#             ncol = 2, nrow = 1,
#             common.legend = TRUE, legend="bottom")
#   pth <- paste("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_rerun",sample_name,".png",sep="")
#   ggsave(pth,
#          width = 7,
#          height = 7)  
# }
# 
# figure <- ggarrange(plotlist=c(plots_no_twopeptides, plots_rsd_twopeptides), 
#           labels = c(1:12),
#           ncol = 2, nrow = 6,
#           common.legend = TRUE, legend="bottom")
# 
# remove_y <- theme(
#   axis.text.y = element_blank(),
#   axis.ticks.y = element_blank(),
#   axis.title.y = element_blank()
# )
# 
# annotate_figure(figure,
#                 top = text_grob("Visualizing Tooth Growth", color = "red", face = "bold", size = 14),
#                 bottom = text_grob("Data source: \n ToothGrowth data set", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10),
#                 left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
#                 right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90),
#                 fig.lab = "Figure 1", fig.lab.face = "bold"
# )
# 
# 
# pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/extracts_try.png"
# ggsave(pth,
#        width = 14,
#        height = 14)  
# 
# 
# 
# 
# 
# # run first qc summary rerun one & two peptides
# p <- ggarrange(p_plot_one, p_plot_two, rsd_plot_one, rsd_plot_two, 
#                labels = c("A", "B", "C", "D"),
#                ncol = 2, nrow = 2,
#                common.legend = TRUE, legend="bottom")
# ggsave("../pics/extracts_rerun.png",
#        width = 7,
#        height = 7)
# 
# #run first qc summary rerun one peptide
# df$rsd <- round(df$rsd,3)
# knitr::kable(df, booktabs = T, "latex")
# #run first qc summary rerun two peptides
# df$rsd <- round(df$rsd,3)
# knitr::kable(df, booktabs = T, "latex")
# 
# 

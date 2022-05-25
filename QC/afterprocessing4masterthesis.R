# run first qc summary rerun one & two peptides
p <- ggarrange(p_plot_one, p_plot_two, rsd_plot_one, rsd_plot_two, 
               labels = c("A", "B", "C", "D"),
               ncol = 2, nrow = 2,
               common.legend = TRUE, legend="bottom")
ggsave("./pics/QC_rerun.png",
       width = 7,
       height = 7)

#run first qc summary rerun one peptide
df$rsd <- round(df$rsd,3)
knitr::kable(df, booktabs = T, "latex")
#run first qc summary rerun two peptides
df$rsd <- round(df$rsd,3)
knitr::kable(df, booktabs = T, "latex")



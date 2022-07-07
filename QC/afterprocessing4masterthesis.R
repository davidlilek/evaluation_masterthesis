# protein identfications and rsd
plot_one <- readRDS("plot_qc_onepeptide.rds")
plot_two <- readRDS("plot_qc_twopeptides.rds")
p <- ggarrange(plot_one, plot_two, 
               ncol = 1, nrow = 2,
               common.legend = TRUE, legend="bottom")
ggsave("../pics/QC_rerun.png",
       width = 7,
       height = 7)

# for second peptides
plot_one <- readRDS("plot_qc_secondpeptides_onepeptide.rds")
plot_two <- readRDS("plot_qc_secondpeptides_twopeptides.rds")
p <- ggarrange(plot_one, plot_two, 
               ncol = 1, nrow = 2,
               common.legend = TRUE, legend="bottom")
ggsave("../pics/QC_rerun_secondpeptides.png",
       width = 7,
       height = 7)


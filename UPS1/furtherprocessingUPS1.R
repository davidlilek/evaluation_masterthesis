library(stringr)
library(dplyr)


########## comparing one and two peptides 
res_2gether_onepeptide$diffonetwo <- res_2gether_onepeptide$NumberProteins/res_2gether_twopeptides$NumberProteins

razor <- subset(res_2gether_onepeptide, diffonetwo > 1 & Type == "Razor")
unique <-subset(res_2gether_onepeptide, diffonetwo > 1 & Type == "Unique")

mean(razor$diffonetwo)
mean(unique$diffonetwo)


plot(razor$Concentration,razor$diffonetwo)
plot(unique$Concentration,unique$diffonetwo)


razor <- subset(res_2gether_onepeptide, diffonetwo > 1 & Type == "Razor")
unique <-subset(res_2gether_onepeptide, diffonetwo > 1 & Type == "Unique")


razor_MBR <- razor %>% filter(str_detect(Evaluation, "^MBR"))
razor_noMBR <- razor %>% filter(grepl("noMBR", Evaluation))

unique_MBR <- unique %>% filter(str_detect(Evaluation, "^MBR"))
unique_noMBR <- unique %>% filter(grepl("noMBR", Evaluation))

plot(razor_MBR$Concentration, razor_MBR$diffonetwo)
plot(unique_MBR$Concentration, unique_MBR$diffonetwo)

plot(razor_noMBR$Concentration, razor_noMBR$diffonetwo)
plot(unique_noMBR$Concentration, unique_noMBR$diffonetwo)

mean(razor_MBR$diffonetwo)
mean(razor_noMBR$diffonetwo)
mean(unique_MBR$diffonetwo)
mean(unique_noMBR$diffonetwo)


############# compare MBR and no MBR


# one peptide

MBR_comparison <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_onepeptide %>% filter(str_detect(Evaluation, "^noMBR"))

MBR_comparison$comparison <- ((MBR_comparison$NumberProteins/noMBR$NumberProteins)-1)*100

plot_compMBR_one <- ggplot(data=MBR_comparison, aes(x=Concentration, y=comparison, group=Name, color = Type, linetype = Evaluation )) +
  geom_line() +
  geom_point() + 
  ylab("Proportion of more proteins\nfound using MBR [%]") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_discrete(labels = c("manually\npooled","pooled\nMaxQuant","pooled\nMaxQuant oldversion"))

# two peptides

MBR_comparison <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^MBR"))
noMBR <- res_2gether_twopeptides %>% filter(str_detect(Evaluation, "^noMBR"))

MBR_comparison$comparison <- ((MBR_comparison$NumberProteins/noMBR$NumberProteins)-1)*100

plot_compMBR_two <- ggplot(data=MBR_comparison, aes(x=Concentration, y=comparison, group=Name, color = Type, linetype = Evaluation )) +
  geom_line() +
  geom_point() + 
  ylab("Proportion of more proteins\nfound using MBR [%]") + 
  xlab("Concentration [fmol]") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  scale_linetype_discrete(labels = c("manually\npooled","pooled\nMaxQuant","pooled\nMaxQuant oldversion"))

ggarrange(plot_compMBR_one, plot_compMBR_two, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend="bottom")
pth <- "N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/ups1_comparisonMBRnoMBR.png"
ggsave(pth,
       width = 7,
       height = 7)  


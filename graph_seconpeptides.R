library(ggplot2)

identified <- c(9946,10525,10916,9804,7250,3698)
second_peptides <- c(304,597,1167,1419,1425,1048)
isolation_window <- c(1,2,4,8,16,32)

df <- as.data.frame(c(identified, second_peptides))
df$Peptides <- c(rep("conventional",6),rep("second",6))
df$window <- as.factor(rep(isolation_window, 2))

colnames(df) <- c("counts", "Peptides", "Window" )


# create stacked barplot
b_plot <- ggplot(df, aes(fill = Peptides, y = counts , x= Window, label = counts)) +
  geom_bar( position = "stack",stat="identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.6)) +
  xlab("Isolation window") +
  labs(fill = "Peptides")
b_plot
ggsave("N:/1_A_Bachelor_Master_Intern/00_M_2022/David/Data/evaluation_masterthesis/pics/graph_secondpeptides.png",
       width = 7,
       height = 3.5)

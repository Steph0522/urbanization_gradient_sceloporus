# Alpha diversity

source("Code/0_loading_data.R")

# effect of depth

alpha_hill_corrplot(table_taxa2)

ggsave("Plots/correlation-depth.png", width = 12, height = 6)

# alpha diversity 
alphahill <- alpha_hill_plot(
  table    = table_taxa2r,
  metadata = metas2r%>% filter(!MUESTRA %in% muestras),
  type     = "boxplot",
  x_col    = "estado2",
  fill_col = "estado2",
  free_y   = TRUE, stat = "wilcox.test",
  save_table = FALSE,
) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.6, size = 2, inherit.aes = TRUE
  )+ theme(legend.position = "none")

alphahill


ggsave(alphahill, file = "Plots/alpha_hill.png", width = 10, height = 6)


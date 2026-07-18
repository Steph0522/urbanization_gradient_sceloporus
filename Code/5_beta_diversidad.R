#beta diversity

source("Code/0_loading_data.R")
library(cowplot)

# nmds with geographical area
beta_estado <- beta_div_plot(
  table      = table_taxa2,
  metadata   = metas2, 
  distance   = "aitchison",
  ordination = "NMDS",
  group_col  = "estado2", legend_title = "Geographical area",
  title = "NMDS - all samples"
)


# continuos nmds
metas3 <- metas2 %>%
  mutate(
    dist_km2 = log(dist_km),
    across(LHC:Elevación, as.numeric)
  )

colors <- viridis::viridis(n = 8)

beta_km <- beta_div_plot(
  table      = table_taxa2,
  metadata   = metas3,
  distance   = "aitchison",
  ordination = "NMDS",
  palette = colors,
  group_col  = "dist_km")+
  ggtitle("NMDS - all samples")



# mds por estado
estados <- unique(metas3$estado2) %>% na.omit()

beta_plots <- lapply(estados, function(est) {
  meta_est <- metas3 %>% filter(estado2 == est)
  tab_est  <- table_taxa2[, c(intersect(meta_est$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

  beta_div_plot(
    table      = tab_est,
    metadata   = meta_est,
    distance   = "aitchison",
    ordination = "NMDS",
    group_col  = "dist_km",
    palette = colors,
  ) + ggtitle(paste("NMDS - ", est))
})

beta_plots[[1]]
beta_plots[[2]]



plot_grid(beta_estado, beta_km, beta_plots[[1]], beta_plots[[2]], align = "hv", labels = c("A", "B", "C", "D"))
ggsave(file = "Plots/beta_nmds.png", width = 16, height = 8, dpi = 300)



# permanova all
MicroBioMeta::beta_test_table(
  table        = table_taxa2,
  metadata     = metas2,
  formula_str  = "dist_km",
  method       = "aitchison",
  test         = "permanova",
  permutations = 999
)

MicroBioMeta::beta_test_table(
  table        = table_taxa2,
  metadata     = metas2,
  formula_str  = "estado2",
  method       = "aitchison",
  test         = "permanova",
  permutations = 999
)

metas3_df <- as.data.frame(metas3)

MicroBioMeta::beta_test_table(
  table        = table_taxa2,
  metadata     = metas3_df,
  formula_str  = "dist_km * estado2",
  method       = "aitchison",
  test         = "permanova",
  permutations = 999
)


# by categorical area

permanova_results <- lapply(estados, function(est) {
  meta_est <- metas3 %>% filter(estado2 == est)
  tab_est  <- table_taxa2[, c(intersect(meta_est$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

  cat("\n====", est, "====\n")
  MicroBioMeta::beta_test_table(
    table        = tab_est,
    metadata     = meta_est,
    formula_str  = "dist_km",
    method       = "aitchison",
    test         = "permanova",
    permutations = 999
  )
})
permanova_results


#differential abundance analysis

source("Code/0_loading_data.R")
library(cowplot)

# ---- Todas las muestras: estado2 (categorico) ----
da_estado2 <- ancombc_plot(
  table        = table_taxa2,
  metadata     = metas2,
  col_cond     = "estado2",
  formula      = "estado2",
  p_adj_method = "BH"
) + ggtitle("All samples - Geographical area")

# ---- Todas las muestras: dist_km (continuo) ----
da_dist_km <- ancombc_plot(
  table        = table_taxa2,
  metadata     = metas2,
  col_cond     = "dist_km",
  formula      = "dist_km",
  p_adj_method = "BH",
  prv_cut      = 0.1
) + ggtitle("All samples - Distance to urban center")

# ---- Puebla and Mexico City: dist_km ----
meta_pue <- metas2 %>% filter(estado2 == "Puebla and Mexico City")
tab_pue  <- table_taxa2[, c(intersect(meta_pue$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

da_pue <- ancombc_plot(
  table        = tab_pue,
  metadata     = meta_pue,
  col_cond     = "dist_km",
  formula      = "dist_km",
  p_adj_method = "BH",
  prv_cut      = 0.1
) + ggtitle("Puebla and Mexico City - Distance to urban center")

# ---- Tlaxcala: dist_km ----
meta_tlax <- metas2 %>% filter(estado2 == "Tlaxcala")
tab_tlax  <- table_taxa2[, c(intersect(meta_tlax$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

da_tlax <- ancombc_plot(
  table        = tab_tlax,
  metadata     = meta_tlax,
  col_cond     = "dist_km",
  formula      = "dist_km",
  p_adj_method = "BH",
  prv_cut      = 0.1
) + ggtitle("Tlaxcala - Distance to urban center")

da_estado2
da_dist_km
da_pue
da_tlax

plot_grid(
  da_estado2, da_dist_km, da_pue, da_tlax,
  ncol   = 2,
  labels = c("A", "B", "C", "D")
)
ggsave(file = "Plots/differential_abundance.png", width = 16, height = 14, dpi = 300)

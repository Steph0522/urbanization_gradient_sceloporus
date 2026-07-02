# ABUNDANCIAS - SCEL_MIGUEL
# Barplot de abundancia relativa (género) y heatmap
# Ejecutar desde el directorio raíz del proyecto (SCEL_MIGUEL.Rproj)

source("Code/00_data_loading.R")

# ---- Barplot de abundancia por categoría de distancia (Puebla y CDMX) ----
barplot_genus <- abundance_barplot(
  table        = table_taxa2,
  metadata     = metas2 %>%
    filter(SAMPLEID != "23") %>%
    filter(estado2 == "Puebla y CDMX"),
  x_axis_title = "Categoría de distancia",
  level        = "genus",
  x_col        = "dist_cat2",
  label        = "Genus",
  width_equal  = FALSE,
  top_n        = 30,
  add_remained = TRUE
)

barplot_genus$data$dist_cat2 <- factor(
  barplot_genus$data$dist_cat2,
  levels = c("Short", "Medium", "Long")
)

barplot_genus
ggsave(barplot_genus, file = "barplot_cat_phylum_2states.png", width = 14, height = 6)

# ---- Heatmap por categoría de distancia ----
metas22      <- metas2 %>% arrange(dist_cat2)
taxonomy     <- table_taxa2$taxonomy
table_taxa22 <- table_taxa2[, match(metas22$SAMPLEID, colnames(table_taxa2)), drop = FALSE]
table_taxa22$taxonomy <- taxonomy

heat <- MicroBioMeta::abundance_heatmap_plot(
  table      = table_taxa22,
  metadata   = metas22,
  condition1 = "dist_cat2",
  top_n      = 30,
  cluster    = TRUE
)
heat
ggsave(heat, file = "heat.png", width = 10, height = 6)

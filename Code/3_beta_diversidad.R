# DIVERSIDAD BETA - SCEL_MIGUEL
# NMDS (Aitchison), PERMANOVA, beta-partición, Venn, dbRDA
# Ejecutar desde el directorio raíz del proyecto (SCEL_MIGUEL.Rproj)

source("Code/00_data_loading.R")
source("beta_div_plot2.R")
library(cowplot)

# ---- NMDS Aitchison por categoría de distancia (Puebla y CDMX) ----
beta_estado <- beta_div_plot(
  table      = table_taxa2,
  metadata   = metas2 %>% filter(estado2 == "Puebla y CDMX"),
  distance   = "aitchison",
  ordination = "NMDS",
  group_col  = "dist_cat2"
)
beta_estado
ggsave(beta_estado, file = "beta_stat2.png", width = 8, height = 4)

# ---- NMDS continuo (dist_km en log) ----
metas3 <- metas2 %>%
  mutate(
    dist_km2 = log(dist_km),
    across(Masa:Elevación, as.numeric)
  )

beta_km <- beta_div_plot2(
  table      = table_taxa2,
  metadata   = metas3,
  distance   = "aitchison",
  ordination = "NMDS",
  group_col  = "dist_km"
)
beta_km
ggsave(beta_km, file = "beta_cat2_km.png", width = 8, height = 4)

# ---- NMDS por estado (log dist_km) ----
estados <- unique(metas3$estado2) %>% na.omit()

beta_plots <- lapply(estados, function(est) {
  meta_est <- metas3 %>% filter(estado2 == est)
  tab_est  <- table_taxa2[, c(intersect(meta_est$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

  beta_div_plot2(
    table      = tab_est,
    metadata   = meta_est,
    distance   = "aitchison",
    ordination = "NMDS",
    group_col  = "dist_km"
  ) + ggtitle(est)
})

beta_plots[[1]]
beta_plots[[2]]

beta_estados <- plot_grid(plotlist = beta_plots, ncol = 2)
ggsave(beta_estados, file = "beta_elevation_estados.png", width = 12, height = 4, dpi = 300)

# ---- PERMANOVA global (dist_km) ----
MicroBioMeta::beta_test_table(
  table        = table_taxa2,
  metadata     = metas2,
  formula_str  = "dist_km",
  method       = "euclidean",
  test         = "permanova",
  permutations = 999
)

# ---- PERMANOVA solo Tlaxcala ----
samples_tlx     <- metas2 %>% filter(estado2 == "Tlaxcala") %>% pull(SAMPLEID)
table_taxa2_tlx <- table_taxa2 %>% dplyr::select(taxonomy, dplyr::all_of(samples_tlx))
metas2_tlx      <- metas2 %>% filter(estado2 == "Tlaxcala")

MicroBioMeta::beta_test_table(
  table        = table_taxa2_tlx,
  metadata     = metas2_tlx,
  formula_str  = "dist_km",
  method       = "euclidean",
  test         = "permanova",
  permutations = 999
)

# ---- PERMANOVA interacción dist_km * estado2 ----
metas3_df <- as.data.frame(metas3)

MicroBioMeta::beta_test_table(
  table        = table_taxa2,
  metadata     = metas3_df,
  formula_str  = "dist_km * estado2",
  method       = "euclidean",
  test         = "permanova",
  permutations = 999
)

# ---- PERMANOVA por estado (loop) ----
permanova_results <- lapply(estados, function(est) {
  meta_est <- metas3 %>% filter(estado2 == est)
  tab_est  <- table_taxa2[, c(intersect(meta_est$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

  cat("\n====", est, "====\n")
  MicroBioMeta::beta_test_table(
    table        = tab_est,
    metadata     = meta_est,
    formula_str  = "dist_km",
    method       = "euclidean",
    test         = "permanova",
    permutations = 999
  )
})

# ---- Beta-partición ----
betapart <- MicroBioMeta::beta_partition_plot(
  table     = table_taxa2,
  metadata  = metas2,
  group_col = "dist_cat2"
)
betapart
ggsave(betapart, file = "betapart.png", width = 10, height = 4)

# ---- Diagrama de Venn ----
MicroBioMeta::venn_diagram_plot(
  table          = table_taxa2,
  metadata       = metas2,
  merge_by       = "dist_cat2",
  min_prevalence = 0,
  method         = "ggvenndiagram"
)

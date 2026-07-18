# DISTANCIA-DECAIMIENTO (Distance Decay)
# Mantel test + decay de similitud (Jaccard, Bray-Curtis, Horn) vs distancia
# geografica, global y por area geografica (estado2).
#beta decay

source("Code/0_loading_data.R")
library(cowplot)


decay_jaccard <- beta_decay_plot(
  table    = table_taxa2r,
  metadata = metas2r,
  lat_col  = "lat",
  lon_col  = "lon",
  distance = "jaccard"
) + ggtitle("All samples")

decay_jaccard

# ---- Decay (Jaccard) por estado: un subset y un plot independiente por estado ----
estados <- unique(metas2r$estado) %>% na.omit()

decay_estado_plots <- lapply(estados, function(est) {
  meta_est <- metas2r %>% filter(estado == est)
  tab_est  <- table_taxa2r[, c(intersect(meta_est$SAMPLEID, colnames(table_taxa2r)), "taxonomy")]

  beta_decay_plot(
    table    = tab_est,
    metadata = meta_est,
    lat_col  = "lat",
    lon_col  = "lon",
    distance = "jaccard"
  ) + ggtitle(est)
})

decay_estado_plots[[1]]
decay_estado_plots[[2]]
decay_estado_plots[[3]]

plot_grid(
  plotlist = c(list(decay_jaccard), decay_estado_plots),
  ncol   = 2,
  labels = c("A", "B", "C", "D")
)
ggsave("Plots/decay_distance_jaccard.png", width = 12, height = 8, dpi = 300)

# ===========================================================================
# PARTE 2: Decay (Bray-Curtis) - todas las muestras
# ===========================================================================

decay_bray <- beta_decay_plot(
  table    = table_taxa2,
  metadata = metas2,
  lat_col  = "lat",
  lon_col  = "lon",
  distance = "bray"
) + ggtitle("Distance decay - Bray-Curtis")

decay_bray

ggsave("Plots/bray_decay.png", width = 7, height = 4, dpi = 300, plot = decay_bray)

# ===========================================================================
# PARTE 3: Decay (Horn / Hill q=1) - global + por estado
# ===========================================================================

decay_horn <- beta_decay_plot(
  table    = table_taxa2,
  metadata = metas2,
  lat_col  = "lat",
  lon_col  = "lon",
  distance = "horn"
) + ggtitle("All samples")

decay_horn

estados_horn <- unique(metas2$estado) %>% na.omit()

decay_horn_estado_plots <- lapply(estados_horn, function(est) {
  meta_est <- metas2 %>% filter(estado == est)
  tab_est  <- table_taxa2[, c(intersect(meta_est$SAMPLEID, colnames(table_taxa2)), "taxonomy")]

  beta_decay_plot(
    table    = tab_est,
    metadata = meta_est,
    lat_col  = "lat",
    lon_col  = "lon",
    distance = "horn"
  ) + ggtitle(est)
})

decay_horn_estado_plots[[1]]
decay_horn_estado_plots[[2]]
decay_horn_estado_plots[[3]]

plot_grid(
  plotlist = c(list(decay_horn), decay_horn_estado_plots),
  ncol   = 2,
  labels = c("A", "B", "C", "D")
)
ggsave("Plots/decay_distance_horn.png", width = 12, height = 8, dpi = 300)

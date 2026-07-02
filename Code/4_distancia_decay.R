# DISTANCIA-DECAIMIENTO (Distance Decay) - SCEL_MIGUEL
# Mantel test, Bray-Curtis/Horn decay global y por estado
# Ejecutar desde el directorio raíz del proyecto (SCEL_MIGUEL.Rproj)

source("Code/00_data_loading.R")
source("functions_betadiv.R")
library(vegan)
library(geosphere)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(broom)

# ===========================================================================
# PARTE 1: Mantel con vegan (Bray-Curtis rarefied vs distancia geográfica)
# ===========================================================================

# ---- Preparar tabla transpuesta ----
tab <- table_taxa2r[, colnames(table_taxa2r) != "taxonomy"]
tab <- tab[, colnames(tab) %in% metas2r$SAMPLEID]
tab <- t(tab)

# ---- Distancia Bray-Curtis ----
bc_dist <- vegdist(tab, method = "bray")

# ---- Matriz de distancias geográficas ----
coords_mat <- metas2r %>%
  filter(SAMPLEID %in% rownames(tab)) %>%
  arrange(match(SAMPLEID, rownames(tab))) %>%
  dplyr::select(lon, lat) %>%
  as.matrix()

geo_dist <- distm(coords_mat) / 1000

sample_names <- metas2r %>%
  filter(SAMPLEID %in% rownames(tab)) %>%
  arrange(match(SAMPLEID, rownames(tab))) %>%
  pull(SAMPLEID)

rownames(geo_dist) <- sample_names
colnames(geo_dist) <- sample_names
geo_dist_vegan     <- as.dist(geo_dist)

# ---- Mantel test global ----
mantel_result <- mantel(bc_dist, geo_dist_vegan, method = "pearson", permutations = 999)
mantel_result

# ---- Plot decay global ----
bc_long <- melt(as.matrix(bc_dist), varnames = c("s1", "s2")) %>%
  mutate(s1 = as.character(s1), s2 = as.character(s2)) %>%
  filter(s1 < s2) %>%
  mutate(Similarity = 1 - value) %>%
  dplyr::select(s1, s2, Similarity)

geo_long <- melt(geo_dist, varnames = c("s1", "s2")) %>%
  mutate(s1 = as.character(s1), s2 = as.character(s2)) %>%
  filter(s1 < s2) %>%
  rename(SpatialDistance = value)

plot_data <- left_join(bc_long, geo_long, by = c("s1", "s2"))

ann_text <- data.frame(
  SpatialDistance = max(plot_data$SpatialDistance) * 0.3,
  Similarity      = 0.9,
  label           = paste0(
    "Mantel: r = ",           signif(mantel_result$statistic, 3),
    ", p-value = ",           signif(mantel_result$signif, 3),
    "\nRegression: slope = ", signif(coef(lm(Similarity ~ SpatialDistance, plot_data))[2], 3)
  )
)

dist_plot <- ggplot(plot_data, aes(x = SpatialDistance, y = Similarity)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  ylab("Bray-Curtis similarity") +
  xlab("Spatial Distance (km)") +
  theme_linedraw() +
  theme(
    legend.position  = "none",
    axis.text        = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
    panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
  ) +
  geom_text(
    data        = ann_text,
    aes(x = SpatialDistance, y = Similarity, label = label),
    size        = 3,
    inherit.aes = FALSE
  )

dist_plot
ggsave(dist_plot, file = "decay_distance.png", width = 4, height = 4, dpi = 300)

# ---- Decay por estado ----
estados <- unique(metas2r$estado) %>% na.omit()

plots_estado <- lapply(estados, function(est) {
  meta_est    <- metas2r %>% filter(estado == est)
  tab_est     <- tab[rownames(tab) %in% meta_est$SAMPLEID, ]

  bc_dist_est <- vegdist(tab_est, method = "bray")

  coords_est <- meta_est %>%
    arrange(match(SAMPLEID, rownames(tab_est))) %>%
    dplyr::select(lon, lat) %>%
    as.matrix()

  geo_dist_est <- distm(coords_est) / 1000
  rownames(geo_dist_est) <- meta_est %>%
    arrange(match(SAMPLEID, rownames(tab_est))) %>% pull(SAMPLEID)
  colnames(geo_dist_est) <- rownames(geo_dist_est)

  mantel_est <- mantel(bc_dist_est, as.dist(geo_dist_est), method = "pearson", permutations = 999)

  bc_long_est <- melt(as.matrix(bc_dist_est), varnames = c("s1", "s2")) %>%
    mutate(s1 = as.character(s1), s2 = as.character(s2)) %>%
    filter(s1 < s2) %>%
    mutate(Similarity = 1 - value) %>%
    dplyr::select(s1, s2, Similarity)

  geo_long_est <- melt(geo_dist_est, varnames = c("s1", "s2")) %>%
    mutate(s1 = as.character(s1), s2 = as.character(s2)) %>%
    filter(s1 < s2) %>%
    rename(SpatialDistance = value)

  plot_data_est <- left_join(bc_long_est, geo_long_est, by = c("s1", "s2"))

  ann_est <- data.frame(
    SpatialDistance = max(plot_data_est$SpatialDistance) * 0.3,
    Similarity      = 0.9,
    label           = paste0(
      "Mantel: r = ",           signif(mantel_est$statistic, 3),
      ", p-value = ",           signif(mantel_est$signif, 3),
      "\nRegression: slope = ", signif(coef(lm(Similarity ~ SpatialDistance, plot_data_est))[2], 3)
    )
  )

  ggplot(plot_data_est, aes(x = SpatialDistance, y = Similarity)) +
    geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    ylab("Bray-Curtis similarity") +
    xlab("Spatial Distance (km)") +
    ggtitle(est) +
    theme_linedraw() +
    theme(
      legend.position  = "none",
      axis.text        = element_text(size = 12),
      panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
      panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
    ) +
    geom_text(
      data        = ann_est,
      aes(x = SpatialDistance, y = Similarity, label = label),
      size        = 3,
      inherit.aes = FALSE
    )
})

plots_estado[[1]]
plots_estado[[2]]
plots_estado[[3]]

ggsave("decay_distance_by_estado.png",
       plot_grid(plotlist = plots_estado, ncol = 2),
       width = 8, height = 4, dpi = 300)

# ===========================================================================
# PARTE 2: Bray-Curtis y Horn con functions_betadiv.R
# ===========================================================================

set.seed(12343)
otu <- table_taxa2 %>% dplyr::select(-taxonomy)
map <- metas2 %>% rename(SampleID = SAMPLEID) %>% mutate_at(c("SampleID"), as.integer)

coords_mat_full <- metas2 %>%
  dplyr::select(SAMPLEID, lon, lat) %>%
  column_to_rownames(var = "SAMPLEID") %>%
  as.matrix()

distance_full <- distm(coords_mat_full) / 1000
colnames(distance_full) <- rownames(coords_mat_full)
rownames(distance_full) <- rownames(coords_mat_full)
distance_complete <- distance_full
distance_full[upper.tri(distance_full)] <- NA

otu_match  <- otu.match(otu)
otu_single <- otu.single(otu_match)
otu_norm   <- otu.norm(otu_single)

bc.dist  <- beta_div_dist_bray(otu_norm)
bc.dist2 <- beta_div_dist_hill(otu_norm, q = 1)

bc.dist.tidy.filt  <- bc.dist.tidy.filter(bc.dist)
bc.dist.tidy.filt2 <- bc.dist.tidy.filter.hill(bc.dist2)

cor_test    <- cor.b(bc.dist.tidy.filt)
lm_test     <- lm.b(bc.dist.tidy.filt)
mantel_test <- mantel.b(bc.dist, distance_complete)

cor_test2    <- cor.b(bc.dist.tidy.filt2)
lm_test2     <- lm.b(bc.dist.tidy.filt2)
mantel_test2 <- mantel.b(bc.dist2, distance_complete)

stats <- data.frame(
  label = paste0(
    "Mantel: r = ",          signif(mantel_test$statistic, 3),
    ", p-value = ",          signif(mantel_test$signif, 3),
    "\nRegression: slope = ", signif(lm_test$estimate, 3)
  )
)

stats2 <- data.frame(
  label = paste0(
    "Mantel: r = ",          signif(mantel_test2$statistic, 3),
    ", p-value = ",          signif(mantel_test2$signif, 3),
    "\nRegression: slope = ", signif(lm_test2$estimate, 3)
  )
)

d1 <- bc.dist.tidy.filt %>%
  ggplot(aes(x = SpatialDistance, y = Similarity)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  ylab("Bray-curtis similarity") +
  xlab("Spatial Distance") +
  theme_linedraw() +
  theme(
    legend.position  = "none",
    axis.text        = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
    panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
  ) +
  annotate("text", label = stats, x = 20, y = 0.7, size = 3.5, colour = "black")

d2 <- bc.dist.tidy.filt2 %>%
  ggplot(aes(x = SpatialDistance, y = Similarity)) +
  geom_point(shape = 16, size = 1, alpha = 0.5, color = "#566573") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  ylab("Horn similarity") +
  xlab("Spatial Distance") +
  theme_linedraw() +
  theme(
    legend.position  = "none",
    axis.text        = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
    panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
  ) +
  annotate("text", label = stats2, x = 20, y = 0.8, size = 3.5, colour = "black")

d1
d2

ggsave("bray.png", width = 7, height = 4, dpi = 300, plot = d1, device = "png")
ggsave("horn.png", width = 7, height = 4, dpi = 300, plot = d2, device = "png")

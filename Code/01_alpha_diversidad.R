# DIVERSIDAD ALFA - SCEL_MIGUEL
# Boxplots de números de Hill por categoría de distancia y por estado
# Ejecutar desde el directorio raíz del proyecto (SCEL_MIGUEL.Rproj)

source("Code/00_data_loading.R")

# ---- Alpha Hill por distancia (rarefied, coloreado por estado) ----
alphahill <- alpha_hill_plot(
  table    = table_taxa2r,
  metadata = metas2r, #%>% filter(!MUESTRA %in% muestras),
  type     = "boxplot",
  x_col    = "estado2",
  fill_col = "dist_cat2",
  free_y   = TRUE
) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.6, size = 2, inherit.aes = TRUE
  )

alphahill$data$dist_cat2 <- factor(
  alphahill$data$dist_cat2,
  levels = c("Short", "Medium", "Long")
)
alphahill + stat_compare_means(label.x = 1.5)

# ---- Alpha Hill por estado (rarefied) ----
alphahill_estado <- alpha_hill_plot(
  table    = table_taxa2r,
  metadata = metas2r %>% filter(!MUESTRA %in% muestras),
  type     = "boxplot",
  x_col    = "estado2",
  fill_col = "estado2",
  free_y   = TRUE
) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.6, size = 2, inherit.aes = TRUE
  ) +
  stat_compare_means(label.x = 1.5) +
  geom_text_repel(
    aes(label = MUESTRA),
    position = position_jitter(width = 0.15, height = 0),
    size = 3,
    show.legend = FALSE
  )

alphahill_estado
ggsave(alphahill_estado, file = "alpha_cat2_rar.png", width = 14, height = 6)

# ---- Alpha sin rarificar - por estado ----
alpha_estado <- alpha_hill_plot(
  table               = table_taxa2,
  metadata            = metas2 %>% filter(!MUESTRA %in% muestras),
  type                = "boxplot",
  fill_col            = "estado2",
  x_col               = "estado2",
  facet_orientation   = "horizontal",
  free_y              = TRUE
) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.6, size = 2, inherit.aes = TRUE
  ) +
  stat_compare_means()

alpha_estado

# ---- Alpha sin rarificar - por categoría de distancia ----
alpha_dist_cat <- alpha_hill_plot(
  table               = table_taxa2,
  metadata            = metas2 %>% filter(!MUESTRA %in% muestras),
  type                = "boxplot",
  fill_col            = "dist_cat2",
  x_col               = "dist_cat2",
  facet_orientation   = "horizontal",
  free_y              = TRUE
) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.6, size = 2, inherit.aes = TRUE
  ) +
  stat_compare_means(label.x = 1.3)

alpha_dist_cat
ggsave(alpha_dist_cat, file = "alpha_cat2_2state.png", width = 14, height = 6)

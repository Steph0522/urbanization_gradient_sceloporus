# correlations 
source("Code/0_loading_data.R")

# load hill data
hills_wide <- read.delim("Data/hill.txt")

hills_long <- hills_wide %>%
  pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "hill")

hills_long3 <- hills_long %>%
  group_by(estado2) %>%
  mutate(dist_km_scaled = scale(dist_km)) %>%
  ungroup()
p_tlax <- alpha_hill_gradient_plot(
  table        = table_taxa2r,
  metadata     = metas2r %>% filter(estado2 == "Tlaxcala"),
  cont_var     = "dist_km",
  method       = "spearman",
  figure_title = "Tlaxcala"
)+
  ggh4x::facetted_pos_scales(
    y = list(
      q == "q0" ~ scale_y_continuous(limits = c(0, 550)),
      q == "q1" ~ scale_y_continuous(limits = c(0, 270)),
      q == "q2" ~ scale_y_continuous(limits = c(0, 160))
      
    )
  )

p_pue <- alpha_hill_gradient_plot(
  table        = table_taxa2r,
  metadata     = metas2r %>% filter(estado2 == "Puebla and Mexico City"),
  cont_var     = "dist_km",
  method       = "spearman",
  figure_title = "Puebla and Mexico City"
)+
  ggh4x::facetted_pos_scales(
    y = list(
      q == "q0" ~ scale_y_continuous(limits = c(0, 550)),
      q == "q1" ~ scale_y_continuous(limits = c(0, 270)),
      q == "q2" ~ scale_y_continuous(limits = c(0, 160))
      
    )
  )

cowplot::plot_grid(p_tlax+theme(axis.title.x = element_blank()), 
                   p_pue+ xlab("Distance to urban center (Km)"), ncol = 1)


ggsave("Plots/alpha_cor.png", width = 14, height = 10)

alpha_hill_gradient_plot(
  table     = table_taxa2r,
  metadata  = metas2r,
  cont_var  = "dist_km",
  group_col = "estado2",
  method    = "spearman"

  ) + theme(legend.position = "none")+ xlab("Distance to urban center (Km)")
ggsave("Plots/alpha_cor2.png", width = 14, height = 6)

alpha_hill_gradient_plot(
  table    = table_taxa2r,
  metadata = metas2r,
  cont_var = "dist_km",
  method   = "spearman"
)+
  ggh4x::facetted_pos_scales(
    y = list(
      q == "q0" ~ scale_y_continuous(limits = c(0, 550)),
      q == "q1" ~ scale_y_continuous(limits = c(0, 270)),
      q == "q2" ~ scale_y_continuous(limits = c(0, 160))
      
    )
  )

# Función auxiliar: calcula stats de regresión y arma el plot
plot_hill_dist <- function(data, titulo = NULL) {
  stats <- data %>%
    group_by(q) %>%
    do({
      m <- lm(hill ~ dist_km, data = .)
      data.frame(
        intercept = coef(m)[1],
        slope     = coef(m)[2],
        r2        = summary(m)$r.squared,
        pval      = summary(m)$coefficients[2, 4],
        corr      = cor(.$hill, .$dist_km, use = "complete.obs")
      )
    }) %>%
    mutate(
      label = sprintf("y = %.3f + %.3f·x\nR² = %.3f\np = %.4f\nr = %.3f",
                      intercept, slope, r2, pval, corr),
      x_pos = -Inf,
      y_pos =  Inf
    )

  p <- ggplot(data, aes(x = dist_km, y = hill)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, color = "#2166AC") +
    facet_wrap(~q, scales = "free_y") +
    labs(
      x     = "Distancia (km)",
      y     = "Número efectivo de ASVs",
      title = if (is.null(titulo)) "Relación entre diversidad alfa (Hill) y distancia" else titulo
    ) +
    geom_text(
      data        = stats,
      aes(x = x_pos, y = y_pos, label = label),
      hjust       = -0.05, vjust = 1.1,
      size        = 3.5,
      inherit.aes = FALSE
    ) +
    theme_bw(base_size = 14)
  p
}

# ---- Plot 1: Todos los datos ----
p_todos <- plot_hill_dist(hills_long, titulo = "Todos los datos")
p_todos
ggsave(p_todos, file = "alpha_dist_todos.png", width = 14, height = 6)

# ---- Plot 2: Solo Tlaxcala ----
p_tlax <- plot_hill_dist(
  hills_long %>% filter(estado2 == "Tlaxcala"),
  titulo = "Tlaxcala"
)
p_tlax
ggsave(p_tlax, file = "alpha_dist_tlaxcala.png", width = 14, height = 6)

# ---- Plot 3: Solo Puebla y CDMX ----
p_pue <- plot_hill_dist(
  hills_long %>% filter(estado2 == "Puebla y CDMX"),
  titulo = "Puebla y CDMX"
)
p_pue
ggsave(p_pue, file = "alpha_dist_puebla.png", width = 14, height = 6)

# ===========================================================================
# REGRESIÓN GLOBAL CON COLOR POR ESTADO (versión con anotaciones por grupo)
# ===========================================================================

# ---- Regresión Hill vs distancia (global, por métrica) ----
stats_df <- hills_long %>%
  group_by(q) %>%
  do({
    m <- lm(hill ~ dist_km, data = .)
    data.frame(
      intercept = coef(m)[1],
      slope     = coef(m)[2],
      r2        = summary(m)$r.squared,
      pval      = summary(m)$coefficients[2, 4],
      corr      = cor(.$hill, .$dist_km, use = "complete.obs")
    )
  }) %>%
  mutate(
    label = sprintf("y = %.3f + %.3f·x\nR² = %.3f\np = %.4f\nr = %.3f",
                    intercept, slope, r2, pval, corr),
    x_pos = -Inf,
    y_pos =  Inf
  )

alpha_dist <- ggplot(hills_long, aes(x = dist_km, y = hill, color = estado2)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  facet_wrap(~q, scales = "free_y") +
  labs(
    x     = "Distancia (km)",
    y     = "Número efectivo de ASVs",
    title = "Relación entre diversidad alfa (Hill) y distancia"
  ) +
  geom_text(
    data        = stats_df,
    aes(x = x_pos, y = y_pos, label = label),
    hjust       = -0.05, vjust = 1.1,
    size        = 4,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 14)

alpha_dist
ggsave(alpha_dist, file = "alpha_dist_rar_cats2.png", width = 14, height = 6)

# ---- Tabla de correlaciones por estado ----
cors_estado <- hills_wide %>%
  group_by(estado2) %>%
  summarise(
    r_q0    = round(cor(dist_km, q0, use = "complete.obs"), 3),
    pval_q0 = round(cor.test(dist_km, q0)$p.value, 4),
    r_q1    = round(cor(dist_km, q1, use = "complete.obs"), 3),
    pval_q1 = round(cor.test(dist_km, q1)$p.value, 4),
    r_q2    = round(cor(dist_km, q2, use = "complete.obs"), 3),
    pval_q2 = round(cor.test(dist_km, q2)$p.value, 4),
    n       = n(),
    .groups = "drop"
  ) %>%
  rename(
    Estado = estado2,
    "r q0" = r_q0, "p q0" = pval_q0,
    "r q1" = r_q1, "p q1" = pval_q1,
    "r q2" = r_q2, "p q2" = pval_q2,
    "n"    = n
  )

tabla <- ggtexttable(cors_estado, rows = NULL, theme = ttheme("mBlue"))
tabla

# ---- Plot correlación por estado y métrica (distancia estandarizada) ----
stats_estado <- hills_long3 %>%
  group_by(estado2, q) %>%
  summarise(
    r    = round(cor(dist_km, hill, use = "complete.obs"), 3),
    pval = round(cor.test(dist_km, hill)$p.value, 4),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, "\np = ", pval))

ggplot(hills_long3, aes(x = dist_km_scaled, y = hill)) +
  geom_point(aes(color = estado2, shape = estado2), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, aes(color = estado2)) +
  geom_text(
    data = stats_estado %>%
      group_by(estado2) %>%
      mutate(dist_km_scaled = 1.5),
    aes(x = dist_km_scaled, y = Inf, label = label),
    vjust = 1.3, size = 3, inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("Puebla y CDMX" = "#E07B54", "Tlaxcala" = "#4A90A4")) +
  scale_shape_manual(values = c("Puebla y CDMX" = 16, "Tlaxcala" = 17)) +
  labs(
    title    = "Alpha diversity and relative urban distance",
    subtitle = "Distance standardized within each state",
    x        = "Distance to urban area (scaled within state)",
    y        = "Hill number",
    color    = "Estado",
    shape    = "Estado"
  ) +
  theme_linedraw() +
  theme(
    legend.position      = "bottom",
    plot.title.position  = "plot",
    axis.text            = element_text(size = 11),
    strip.text           = element_text(size = 11, face = "bold.italic"),
    panel.grid.major     = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
    panel.grid.minor     = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
  ) +
  facet_grid(q ~ estado2, scales = "free_y")

# ---- Corrplot: Hill q0 vs variables morfométricas ----
cors_q0 <- hills_long %>%
  filter(q == "q0") %>%
  mutate(across(TC:Elevación, as.numeric)) %>%
  dplyr::select(Masa:Elevación, dist_km, hill)

testRes <- cor.mtest(cors_q0, conf.level = 0.95)
M       <- cor(cors_q0, use = "pairwise.complete.obs")

png(filename = "cors1.png", height = 450)
corrplot(M,
         method      = "circle",
         type        = "lower",
         insig       = "blank",
         addCoef.col = "black",
         number.cex  = 0.8,
         order       = "AOE",
         diag        = FALSE,
         p.mat       = testRes$p,
         sig.level   = 0.05)
dev.off()

# ---- Alpha Hill corrplot (MicroBioMeta) ----


# ---- Correlación con variables ambientales (phylum level) ----
env <- metas2 %>%
  dplyr::select(SAMPLEID, dist_km, Alto_cabeza, Regeneración, LHC, Cola, Elevación, TC, TA, TS) %>%
  tibble::column_to_rownames("SAMPLEID") %>%
  as.data.frame()

table_taxa2 <- as.data.frame(table_taxa2)
metas2      <- as.data.frame(metas2)
env         <- as.data.frame(env)

MicroBioMeta::corr_env_abund_plot(
  table     = table_taxa2,
  env_table = env,
  metadata  = metas2,
  level     = "phylum"
)

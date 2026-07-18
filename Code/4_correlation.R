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

#plots separated

p_tlax <- alpha_decay_plot(
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

p_pue <- alpha_decay_plot(
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


# plots joined
alpha_decay_plot(
  table     = table_taxa2r,
  metadata  = metas2r,
  cont_var  = "dist_km",
  group_col = "estado2",
  method    = "spearman"

  ) + theme(legend.position = "none")+ xlab("Distance to urban center (Km)")
ggsave("Plots/alpha_cor2.png", width = 14, height = 6)

alpha_decay_plot(
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



# correlation with other
env <- meta_morfo2 %>%
  dplyr::select(MUESTRA,Alto_cabeza, Regeneración, LHC, Cola, Elevación, TC, TA, TS) %>%
  tibble::column_to_rownames("MUESTRA") %>%
  as.data.frame()

table_taxa2 <- as.data.frame(table_taxa2)
metas2      <- as.data.frame(metas2)
env         <- as.data.frame(env)

MicroBioMeta::corr_env_abund_plot(
  table     = table_taxa2,
  env_table = env,
  metadata  = meta_morfo2,
  level     = "phylum"
)

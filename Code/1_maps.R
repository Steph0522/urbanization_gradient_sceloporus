library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)
library(ggspatial)


# sampling and sites map
bbox <- c(left = -99.5, bottom = 18.8, right = -97.8, top = 19.7)

mapa_base <- get_stadiamap(
  bbox  = bbox,
  zoom  = 10,
  maptype = "stamen_terrain" 
)

mapa <- ggmap(mapa_base,  darken = c(0.3, "white")) +
  geom_point(
    data = meta %>% filter(!is.na(lon)),
    aes(x = lon, y = lat, fill = estado2),
    shape  = 21,     
    size   = 3,
    alpha  = 0.9,
    color  = "black", 
    stroke = 0.8      
  ) +
  scale_fill_manual(
    values = c("Tlaxcala" = "#E69F00", "Puebla and Mexico City" = "#56B4E9"),
    name   = "Geographical area"
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", axis.text = element_text(colour = "black"))

mapa
ggsave(mapa, file = "Plots/map_samples.png", width = 10, height = 6, dpi = 300)



# distribution and urban points map
metas2r_filt <- metas2r %>%
  group_by(estado2) %>%
  mutate(
    lon_centro = mean(lon, na.rm = TRUE),
    lat_centro = mean(lat, na.rm = TRUE),
    x_km  = distHaversine(cbind(lon_centro, lat),  cbind(lon,  lat))  * ifelse(lon  > lon_centro, 1, -1),
    y_km  = distHaversine(cbind(lon,  lat_centro), cbind(lon,  lat))  * ifelse(lat  > lat_centro, 1, -1),
    x_km2 = distHaversine(cbind(lon_centro, lat2), cbind(lon2, lat2)) * ifelse(lon2 > lon_centro, 1, -1),
    y_km2 = distHaversine(cbind(lon2, lat_centro), cbind(lon2, lat2)) * ifelse(lat2 > lat_centro, 1, -1)
  ) %>%
  ungroup()

ggplot() +
  geom_segment(
    data = metas2r_filt,
    aes(x = lon, y = lat, xend = lon2, yend = lat2, color = dist_km),
    alpha = 0.5, linewidth = 0.5
  ) +
  geom_point(
    data  = metas2r_filt,
    aes(x = lon, y = lat, color = dist_km, size = dist_km),
    alpha = 0.85
  ) +
  geom_point(
    data  = metas2r_filt,
    aes(x = lon2, y = lat2),
    shape = 4, size = 3, color = "black"
  ) +
  scale_color_gradientn(
    colors = c("#2166AC", "#74ADD1", "#FEE090", "#F46D43", "#A50026"),
    name   = "Distance to\nurban center (km)"
  ) +
  facet_wrap(~estado2, scales = "free") +
  labs(x = "lon", y = "lat") +
  theme_linedraw() +
  theme(
    legend.position  = "right",
    strip.text       = element_text(size = 12, face = "bold.italic"),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
    panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
  ) +
  geom_text_repel(
    data  = metas2r_filt,
    aes(x = lon, y = lat, label = MUESTRA),
    size  = 3, color = "black"
  )

ggsave("Plots/map_km_urban.png", width = 10, height = 5, dpi = 300)

# boxplot states
metas2r %>%
  ggplot(aes(x = estado2, y = dist_km, fill = estado2)) +
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Distrito Federal" = "#E07B54", "Tlaxcala" = "#4A90A4")) +
  labs(
    x    = NULL,
    y    = "Distance to nearest urban center (km)",
    fill = NULL
  ) +
  theme_linedraw() +
  theme(
    legend.position  = "none",
    axis.text        = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "#E5E8E8"),
    panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "#E5E8E8")
  ) +
  stat_compare_means(method = "wilcox.test") +
  ggrepel::geom_text_repel(aes(label = MUESTRA), size = 3)

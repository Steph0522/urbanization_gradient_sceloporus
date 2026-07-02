# MAPAS - SCEL_MIGUEL
# Mapa interactivo leaflet + mapa de segmentos + distribución de distancias por estado
# Ejecutar desde el directorio raíz del proyecto (SCEL_MIGUEL.Rproj)

source("Code/00_data_loading.R")
library(leaflet)
library(geosphere)

# ---- Resumen de muestras por grupo ----
metas  %>% drop_na() %>% group_by(dist_cat)  %>% count_()
metas  %>% drop_na() %>% group_by(dist_cat2) %>% count_()
meta   %>% group_by(estado2) %>% count_()
meta   %>% group_by(estado2, dist_cat2) %>% count_()
meta   %>% group_by(estado) %>% count_()

# ---- Mapa interactivo (leaflet) - animal vs punto urbano ----
leaflet(meta) %>%
  addTiles() %>%
  addCircleMarkers(
    lng  = ~lon,
    lat  = ~lat,
    color = ~case_when(
      dist_cat2 == "Short"  ~ "blue",
      dist_cat2 == "Medium" ~ "orange",
      dist_cat2 == "Long"   ~ "red",
      TRUE                  ~ "gray"
    ),
    popup = ~paste(
      "ID:", ID, "<br>",
      "Estado:", estado, "<br>",
      "Distancia:", round(dist_km, 2), "km"
    ),
    radius      = 6,
    fillOpacity = 0.8,
    group       = "Animales"
  ) %>%
  addCircleMarkers(
    lng  = ~lon2,
    lat  = ~lat2,
    color       = "green",
    popup       = ~paste("Punto urbano de:", ID, "<br>", "Estado:", estado),
    radius      = 6,
    fillOpacity = 0.8,
    group       = "Comunidades urbanas"
  ) %>%
  addLegend(
    position = "bottomright",
    colors   = c("blue", "orange", "red", "green"),
    labels   = c("Short", "Medium", "Long", "Comunidad urbana"),
    title    = "Distancia / Tipo de punto"
  ) %>%
  addLayersControl(
    overlayGroups = c("Animales", "Comunidades urbanas"),
    options       = layersControlOptions(collapsed = FALSE)
  )

# ---- Mapa con segmentos animal → punto urbano (por estado) ----
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

ggsave("map_km_estados.png", width = 10, height = 5, dpi = 300)

# ---- Distribución de distancias por estado (boxplot + puntos) ----
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

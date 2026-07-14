# Load data and pre-processing

# Load libraries
#library(devtools)
#install_github("Steph0522/MicroBioMeta")
library(qiime2R)
library(tidyverse)
library(readxl)
library(dplyr)
library(sp)
library(sf)
library(geosphere)
library(ggpubr)
library(ggrepel)
library(corrplot)
library(leaflet)
library(rnaturalearth)
library(devtools)
library(MicroBioMeta)

#load data
table  <- read_qza("Data/run_f250_r230_feature-table.qza")$data %>% as.data.frame()
tabler <- read_qza("Data/rarefied_table.qza")$data %>% as.data.frame()


# Pre-processing

colnames(table)  <- sub(".*\\.", "", colnames(table))
colnames(tabler) <- sub(".*\\.", "", colnames(tabler))

depth <- read.delim("Data/depth.csv", sep = ";")
depth$Sample.ID <- colnames(table)

metadata <- read.csv("Data/metadata_ids.csv", sep = ";") %>%
  mutate_all(as.character) %>%
  dplyr::select(ID, everything(), -ID.CAM)

meta_morfo <- read_excel("Data/metadata_final.xlsx") %>%
  mutate(sex = case_when(Sex == "HG" ~ "H", TRUE ~ as.character(Sex))) %>%
  mutate_all(as.character)

metadatas <- meta_morfo %>%
  filter(MUESTRA != "23") %>%
  right_join(metadata, by = c("MUESTRA", "ID")) %>%
  mutate(
    across(
      c(ID, lat, lon, lat2, lon2, TC:Elevación),
      ~as.numeric(na_if(trimws(.x), "NA"))
    )
  ) %>%
  arrange(ID) %>%
  rowwise() %>%
  mutate(dist_km = distHaversine(c(lon, lat), c(lon2, lat2)) / 1000) %>%
  ungroup() %>%
  mutate(
    dist_category = case_when(
      dist_km < 1                   ~ "High",
      dist_km >= 1 & dist_km <= 100 ~ "Medium",
      dist_km > 100                 ~ "Low",
      TRUE                          ~ NA_character_
    )
  ) %>%
  mutate_at(c("ID"), as.character) %>%
  rename(SAMPLEID = INDEX) %>%
  dplyr::select(SAMPLEID, everything()) %>%
  filter(!MUESTRA == "23")

metadatas$dist_cat <- cut(
  metadatas$dist_km,
  breaks = seq(min(metadatas$dist_km, na.rm = TRUE),
               max(metadatas$dist_km, na.rm = TRUE),
               length.out = 4),
  labels = c("Short", "Medium", "Long"),
  include.lowest = TRUE
)

metadatas$dist_cat2 <- cut(
  metadatas$dist_km,
  breaks = quantile(metadatas$dist_km, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  labels = c("Short", "Medium", "Long"),
  include.lowest = TRUE
)

metadatas$dist_cat  <- factor(metadatas$dist_cat,  levels = c("Short", "Medium", "Long"))
metadatas$dist_cat2 <- factor(metadatas$dist_cat2, levels = c("Short", "Medium", "Long"))

metadata_depth <- metadata %>%
  full_join(depth, by = c("ID" = "Sample.ID")) %>%
  filter(!is.na(Frequency) & Frequency > 5000)

metas  <- metadatas[match(sub(".*\\.", "", colnames(read_qza("Data/run_f250_r230_feature-table.qza")$data)),
                          metadatas$SAMPLEID), ]
metasr <- metadatas[match(colnames(tabler), metadatas$SAMPLEID), ]

meta_sf    <- st_as_sf(metas %>% drop_na(), coords = c("lon", "lat"), crs = 4326)
mex_states <- ne_states(country = "Mexico", returnclass = "sf")
meta_sf    <- st_join(meta_sf, mex_states["name"])

meta <- metas %>%
  drop_na() %>%
  mutate(
    estado  = meta_sf$name,
    estado2 = factor(case_when(
      estado == "Puebla"           ~ "Puebla and Mexico City",
      estado == "Distrito Federal" ~ "Puebla and Mexico City",
      estado == "Tlaxcala"         ~ "Tlaxcala"
    ))
  )

taxonomy_gg2 <- read_qza("Data/taxonomy_gg2_weighted.qza")$data %>%
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)

table_taxa  <- merge_feature_taxonomy(table,  taxonomy_gg2) %>%
  dplyr::select(-any_of(c("32", "38", "117")))

table_taxar <- merge_feature_taxonomy(tabler, taxonomy_gg2) %>%
  dplyr::select(-any_of(c("32", "38", "117")))

EXCLUDE <- c("32", "38", "117")

metas2 <- meta %>%
  filter(!is.na(SAMPLEID), !SAMPLEID %in% EXCLUDE) %>%
  mutate(
    dist_cat  = factor(dist_cat,  levels = c("Short", "Medium", "Long")),
    dist_cat2 = factor(dist_cat2, levels = c("Short", "Medium", "Long"))
  )

align_table <- function(tbl, ids) {
  cols_ok <- c(intersect(ids, colnames(tbl)), "taxonomy")
  tbl[, cols_ok]
}

table_taxa2  <- align_table(table_taxa,  metas2$SAMPLEID)
table_taxa2r <- align_table(table_taxar, metas2$SAMPLEID)

samples_ok <- intersect(metas2$SAMPLEID, colnames(table_taxa2r))
metas2r    <- metas2 %>% filter(SAMPLEID %in% samples_ok)

muestras <- c("54", "53")

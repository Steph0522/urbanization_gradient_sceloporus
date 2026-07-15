# relative abundance

source("Code/0_loading_data.R")


#barplot phylum and genus

barplot_phylum <- abundance_barplot(
  table        = table_taxa2,
  metadata     = metas2 %>%
   filter(SAMPLEID != "23"), #%>%
  #    filter(estado2 == "Puebla y CDMX"),
  x_axis_title = "Geographical area",
  level        = "phylum",
  x_col        = "estado2",
  label        = "Phylum",
  width_equal  = FALSE,
#  top_n        = 30,
  add_remained = TRUE, 
save_table = FALSE,
)


barplot_phylum



barplot_genus <- abundance_barplot(
  table        = table_taxa2,
  metadata     = metas2, 
  x_axis_title = "Geographical area",
  level        = "genus",
  x_col        = "estado2",
  label        = "Genus",
  width_equal  = FALSE,
  top_n        = 20,
  add_remained = TRUE
)


barplot_genus


both <- cowplot::plot_grid(barplot_phylum, barplot_genus, labels = c("A", "B"))

both
ggsave(both, file = "Plots/barplot.png", width = 14, height = 6)


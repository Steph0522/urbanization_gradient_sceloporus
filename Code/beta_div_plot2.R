#' Beta diversity plot with multiple distance and ordination methods
#'
#' This function computes beta diversity using several distance metrics and ordination methods
#' (PCA, PCoA, NMDS). It requires an abundance table with taxonomy, metadata, and allows customization
#' of color and shape aesthetics. It also supports compositional transformation via ALDEx2.
#'
#' @param table A data frame with abundances. The last column must contain taxonomy information.
#' @param metadata A data frame with sample metadata. The first column must contain the sample IDs.
#' @param distance Distance method: one of "euclidean", "bray", "jaccard", "sorensen",
#'  "compositional" (default), "aitchison", or "robust.aitchison".
#' @param ordination Ordination method: one of "PCA" (default), "PCoA", or "NMDS".
#' @param group_col Column in `metadata` to fill color points.
#' @param group_colors Optional named vector of colors.
#' @param shape_col Optional column in `metadata` to shape points.
#' @param legend_title Optional legend title.
#' @param n_taxa Number of top contributing taxa to display as arrows in PCA.
#'
#' @return A `ggplot2` object.
#' @export

beta_div_plot2 <- function(table, metadata, 
                           distance = "compositional",
                           ordination = "PCA", 
                           group_col = NULL,
                           group_colors = NULL,
                           shape_col = NULL,
                           legend_title = NULL,
                           arrows_size = 10,
                           n_taxa = 5) {
  
  requireNamespace("vegan")
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
  requireNamespace("ALDEx2")
  requireNamespace("dplyr")
  requireNamespace("stringr")
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  if (ordination == "PCA" && distance != "compositional") {
    stop("PCA is only available with 'compositional' (Aitchison). Use PCoA or NMDS instead.")
  }
  
  tax_col <- names(table)[ncol(table)]
  taxonomy <- table[[tax_col]]
  abund_table <- table[, -ncol(table)]
  
  if ("Feature.ID" %in% names(abund_table)) {
    feature_ids <- abund_table$Feature.ID
    abund_table <- abund_table[, !(names(abund_table) == "Feature.ID")]
  } else {
    feature_ids <- rownames(abund_table)
  }
  rownames(abund_table) <- feature_ids
  otu_table <- as.data.frame(lapply(abund_table, as.numeric), check.names = FALSE)
  rownames(otu_table) <- feature_ids
  
  
  metadata_ids <- trimws(as.character(metadata[[1]]))
  sample_ids <- trimws(colnames(otu_table))
  colnames(otu_table) <- sample_ids
  metadata[[1]] <- metadata_ids
  
  common_samples <- intersect(sample_ids, metadata_ids)
  if (length(common_samples) == 0) stop("No matching samples between table and metadata.")
  otu_table <- otu_table[, common_samples]
  metadata <- metadata[metadata_ids %in% common_samples, ]
  
  if (distance == "compositional") {
    set.seed(123)
    aldex_obj <- ALDEx2::aldex.clr(otu_table, mc.samples = 128,
                                   denom = "all", verbose = FALSE, useMC = FALSE)
    otu_trans <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))
    dist_matrix <- dist(otu_trans, method = "euclidean")
  } else if (distance %in% c("aitchison", "robust.aitchison")) {
    dist_matrix <- vegan::vegdist(t(otu_table), method = distance, pseudocount = 0.5)
    otu_trans <- NULL
  } else {
    dist_matrix <- vegan::vegdist(t(otu_table), method = distance)
    otu_trans <- NULL
  }
  
  expl_var <- NULL
  res <- switch(ordination,
                "PCA" = {
                  pca_input <- if (!is.null(otu_trans)) otu_trans else t(otu_table)
                  pca <- prcomp(pca_input)
                  list(ord = pca, expl_var = round(100 * summary(pca)$importance[2, 1:2], 1))
                },
                "PCoA" = {
                  pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
                  eigs <- pcoa$eig
                  list(ord = pcoa, expl_var = round(100 * eigs[1:2] / sum(eigs[eigs > 0]), 1))
                },
                "NMDS" = list(ord = vegan::metaMDS(dist_matrix, k = 2, trymax = 100), expl_var = NULL),
                stop("Invalid ordination method.")
  )
  
  ord_res <- res$ord
  expl_var <- res$expl_var
  
  
  ord_df <- switch(ordination,
                   "PCA" = as.data.frame(ord_res$x),
                   "PCoA" = {
                     df <- as.data.frame(ord_res$points)
                     colnames(df) <- c("PCoA1", "PCoA2")  
                     df
                   },
                   "NMDS" = as.data.frame(ord_res$points))
  
  ord_df$SampleID <- rownames(ord_df)
  colnames(metadata)[1] <- "SampleID"
  merged <- dplyr::inner_join(ord_df, metadata, by = "SampleID")
  if (nrow(merged) == 0) stop("Ninguna muestra en común entre tabla y metadata.")
  
  x_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[1], " (", expl_var[1], "%)") else names(ord_df)[1]
  y_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[2], " (", expl_var[2], "%)") else names(ord_df)[2]
  
  if (is.null(group_colors)) {
    group_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  }
  legend_name <- ifelse(is.null(legend_title), group_col, legend_title)
  fill_scale <- ggplot2::scale_fill_manual(name = legend_name, values = group_colors)
  
  if (is.null(shape_col)) {
    p <- ggplot2::ggplot(merged, ggplot2::aes_string(
      x = names(ord_df)[1],
      y = names(ord_df)[2],
      fill= group_col
    )) +
      ggplot2::geom_point(size = 4, shape = 21)
  } else {
    p <- ggplot2::ggplot(merged, ggplot2::aes_string(
      x = names(ord_df)[1],
      y = names(ord_df)[2],
      color = group_col,
      shape = shape_col
    )) +
      ggplot2::geom_point(size = 4)
  }
  
  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black", size = 12, family = "serif"),
      axis.text.y = ggplot2::element_text(color = "black", size = 12,  family = "serif"),
      axis.title.x = ggplot2::element_text(color = "black", size = 14, family = "serif"),
      axis.title.y = ggplot2::element_text(color = "black", size = 14, family = "serif"),
      legend.text = ggplot2::element_text(color = "black", size = 12, family = "serif"),
      legend.title = ggplot2::element_text(color = "black", size = 14, family = "serif", face = "bold"),
      plot.title = ggplot2::element_text(color = "black", size = 16, family = "serif", face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(paste(ordination, "-", distance)) +
    
    # 👇👇 NUEVA CAPA: muestra el nombre del SampleID junto a cada punto
   # ggplot2::geom_text(ggplot2::aes(label = SampleID), vjust = -1, size = 3, family = "serif")
  
  if (ordination == "PCA") {
    rot_df <- as.data.frame(ord_res$rotation)
    rot_df$Feature.ID <- rownames(rot_df)
    rot_df$mag <- sqrt(rot_df$PC1^2 + rot_df$PC2^2)
    rot_df <- rot_df[order(rot_df$mag, decreasing = TRUE), ][1:n_taxa, ]
    rot_df$PC1 <- rot_df$PC1 * arrows_size
    rot_df$PC2 <- rot_df$PC2 * arrows_size
    rot_df$Taxon <- taxonomy[match(rot_df$Feature.ID, feature_ids)]
    rot_df$label <- stringr::str_extract(rot_df$Taxon, "(?<=__)[^;]*$") 
    rot_df$label <- gsub(" ", "\n", rot_df$label)
    
    p <- p +
      ggplot2::geom_segment(data = rot_df,
                            aes(x = 0, y = 0, xend = PC1, yend = PC2),
                            arrow = ggplot2::arrow(length = unit(0.7, "cm")),
                            color = "gray30",
                            inherit.aes = FALSE) #+
      #ggrepel::geom_label_repel(data = rot_df,
       #                         aes(x = PC1, y = PC2, label = label),
        #                        fill = "white", color = "black",
         #                       fontface = "italic", size = 4, family= "serif",
          #                      inherit.aes = FALSE)
  }
  
  # --- Centrar ejes simétricamente ---
  x_limits <- range(merged[[names(ord_df)[1]]], na.rm = TRUE)
  y_limits <- range(merged[[names(ord_df)[2]]], na.rm = TRUE)
  max_range <- max(abs(x_limits), abs(y_limits))
  p <- p + ggplot2::coord_cartesian(xlim = c(-max_range, max_range),
                                    ylim = c(-max_range, max_range)) +
    ggplot2::coord_fixed()
  
  return(p)
}

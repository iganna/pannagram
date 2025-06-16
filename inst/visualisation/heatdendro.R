heatdendro <- function(mx,
                       show_xticks = FALSE,
                       n_clust = 0,
                       col_clust = c("#B56E6B", "grey70")) {
  # Check rownames and colnames, assign default names if NULL
  if (is.null(rownames(mx))) {
    rownames(mx) <- paste0("Row", seq_len(nrow(mx)))
  }
  if (is.null(colnames(mx))) {
    colnames(mx) <- paste0("Col", seq_len(ncol(mx)))
  }
  
  nx <- ncol(mx) * 2
  ny <- nrow(mx) * 2
  
  # Perform clustering
  row_dend <- as.dendrogram(hclust(dist(mx)))
  col_dend <- as.dendrogram(hclust(dist(t(mx))))
  
  if(n_clust > 1){
    # Color branches of the row dendrogram and customize appearance
    row_dend <- color_branches(row_dend, k = 2, col = col_clust)
    row_dend <- set(row_dend, "labels", rep("", length(labels(row_dend))))
    row_dend <- set(row_dend, "branches_lwd", 0.5)  
  }
  
  # Get order of rows and columns from dendrograms
  row_order <- order.dendrogram(row_dend)
  col_order <- order.dendrogram(col_dend)
  
  mx.ordered <- mx[row_order, col_order]
  
  # Prepare data for ggplot
  df <- as.data.frame(mx.ordered)
  df$row <- factor(rownames(df), levels = rownames(mx.ordered))
  df_long <- pivot_longer(df, -row, names_to = "col", values_to = "value")
  df_long$col <- factor(df_long$col, levels = colnames(mx.ordered))
  
  # Heatmap with y-axis labels on the right side
  heatmap_gg <- ggplot(df_long, aes(x = col, y = row, fill = as.factor(value))) +
    geom_tile() +
    scale_fill_manual(
      values = c("0" = "#F8F2DE", "1" = "#A31D1D"),
      name = "ORF",
      labels = c("0" = "absent", "1" = "present")
    ) +
    theme_void() +  
    theme(
      axis.text.x = if (show_xticks) element_text(angle = 90, hjust = 1, vjust = 0.5) else element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    scale_y_discrete(expand = c(1/ny, 1/ny)) +
    scale_x_discrete(expand = c(1/nx, 1/nx))
  
  # Column dendrogram (on top)
  dend_col_gg <- ggdendrogram(col_dend, rotate = FALSE, theme_dendro = FALSE) +
    theme_void() +
    scale_x_continuous(expand = c(1/nx, 1/nx))  # Remove padding
  
  # Row dendrogram (on left)
  
  if(n_clust > 1){
    ggdend <- as.ggdend(row_dend)
    dend_row_gg <- ggplot(ggdend, theme = theme_void()) +
      coord_flip() +
      scale_x_continuous(expand = c(1/ny, 1/ny)) +
      scale_y_reverse() +
      theme_void()
  } else {
    dend_row_gg <- ggdendrogram(row_dend, rotate = TRUE, theme_dendro = FALSE) +
      theme_void() +
      scale_x_continuous(expand = c(1/ny, 1/ny)) + 
      scale_y_reverse()
  }
  
  
  
  empty_bottom <- plot_spacer()
  
  # Compose the final layout
  final_plot <- empty_bottom + plot_spacer() + dend_col_gg + empty_bottom +
    empty_bottom + dend_row_gg + heatmap_gg + empty_bottom +
    empty_bottom + plot_spacer() + empty_bottom + empty_bottom +
    plot_layout(
      ncol = 4, nrow = 3,
      widths = c(0.1, 1, 4, 0.1),
      heights = c(0.6, 4, 0.1)
    ) &
    theme(plot.margin = unit(rep(0, 4), "pt"))
  
  return(final_plot)
}

library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

countplot <- function(data, sv.class, show.legend = F, colormap = NULL) {

  # Collapse each row into a unique key and count occurrences
  keys <- apply(data, 1, paste, collapse = "_")
  counts <- table(keys)
  
  # Reconstruct the unique rows from the keys
  data.uni <- do.call(rbind, strsplit(names(counts), "_"))
  data.uni <- apply(data.uni, 2, as.numeric)
  colnames(data.uni) <- colnames(data)
  names(counts) <- NULL
  
  # Order columns by sv.class
  sv.class <- sv.class[colnames(data.uni)]
  species.order <- names(sv.class)[order(sv.class)]
  data.uni <- data.uni[, species.order]
  sv.class <- sv.class[species.order]
  
  # Cluster and order rows by similarity
  dist_rows <- dist(data.uni)
  hc <- hclust(dist_rows)
  row.order <- hc$order
  data.uni <- data.uni[row.order, ]
  counts <- counts[row.order]
  
  # Convert to long format for ggplot
  df <- as.data.frame(t(data.uni))
  df$row <- seq_len(nrow(df))
  df$class <- sv.class  
  df.long <- reshape2::melt(df, id.vars = c("row", "class"), 
                            variable.name = "column", value.name = "value")
  df.dot <- df.long[df.long$value == 1, ]
  df.dot$column <- as.numeric(gsub("[^0-9]", "", df.dot$column))
  df.dot$row.label <- colnames(data.uni)[df.dot$row]
  
  # Ensure row.label is a factor with correct ordering
  df.dot$row.label <- factor(df.dot$row.label, levels = rev(species.order))
  
  # Calculate y-ranges for vertical segments
  segments_df <- df.dot %>%
    group_by(column) %>%
    summarise(
      y_min = min(as.numeric(row.label)),
      y_max = max(as.numeric(row.label))
    )
  
  # Create alternating background stripes
  y_levels <- levels(df.dot$row.label)
  y_numeric <- seq_along(y_levels)
  strip_rows <- y_numeric[y_numeric %% 2 == 1]
  
  background_df <- data.frame(
    ymin = strip_rows - 0.5,
    ymax = strip_rows + 0.5,
    fill = rep(#c("gray95", "#f0fcfd"), 
               c("gray95", "gray95"),
               length.out = length(strip_rows))  # чередующиеся цвета
  )
  
  
  # All possible positions for light gray background dots
  all_positions <- expand.grid(
    column = unique(df.dot$column),
    row.label = levels(df.dot$row.label)
  )
  all_positions$y_numeric <- as.numeric(factor(all_positions$row.label, levels = levels(df.dot$row.label)))
  
  x_vals <- sort(unique(df.dot$column))
  x_min <- min(x_vals) - 0.5
  x_max <- max(x_vals) + 0.5
  
  # Main dot plot
  # legend_required <- !all(names(sv.class) == sv.class) & show.legend
  legend_required = show.legend
  main_plot <- ggplot(df.dot, aes(x = column, y = as.numeric(row.label), color = class)) +
    geom_rect(data = background_df,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill),
              color = NA, inherit.aes = FALSE) +
    scale_fill_identity() +
    geom_point(data = all_positions,
               aes(x = column, y = y_numeric),
               color = "gray85", size = 1, inherit.aes = FALSE) +
    geom_segment(data = segments_df, aes(x = column, xend = column, y = y_min, yend = y_max),
                 color = "gray40", inherit.aes = FALSE, size = 1) +
    geom_point(size = 3) +
    scale_y_continuous(breaks = y_numeric, labels = y_levels) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    theme_minimal() +
    labs(x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      legend.position = if (legend_required) "bottom" else "none",
      legend.title = element_blank(),
      legend.box = "horizontal"
    )
  
  if(!is.null(colormap)){
    main_plot = main_plot + scale_color_manual(values = colormap)
  }
  
  # Barplot for counts
  df.counts <- data.frame(
    column = 1:length(counts), 
    count = as.numeric(unname(counts))
  )
  
  bar_plot <- ggplot(df.counts, aes(x = column, y = count)) +
    geom_col(fill = "grey50", width = 0.6) +
    geom_text(aes(label = count), vjust = -0.5, size = 3) +
    theme_minimal() +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max(df.counts$count) * 1.4), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  # Combine both plots vertically
  combined_plot <- plot_grid(
    bar_plot + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, -5, 0)
    ),
    main_plot + theme(
      plot.margin = margin(-5, 0, 0, 0)
    ) + ggtitle(NULL),
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(1, 4)
  )
  
  return(combined_plot)
}

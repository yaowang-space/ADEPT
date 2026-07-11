# ADEPT: Automated Depth Profiling Technique
# plotting.R — ggplot2 visualization
#
# Publication-ready depth profile plots with plateau overlays.

#' Create a depth profile plot
#'
#' Generates a ggplot showing raw ages, LOESS-smoothed trend, and
#' final plateau segments with uncertainty bands.
#'
#' @param df data.frame with Time, Raw_Age, loess_Age columns
#' @param segments Segment data.frame with Start, End, Segment_Mean,
#'   Total_uncertainty, and Filter_4 columns
#' @param title Plot title (defaults to Analysis name)
#' @return A ggplot object
#' @keywords internal
plot_depth_profile <- function(df, segments, title = NULL) {
  if (is.null(title)) {
    title <- paste("Analysis:", df[1, "Analysis"])
  }

  final_color <- ifelse(!is.na(segments$Filter_4), "red", NA)

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df,
                        ggplot2::aes(x = .data$Time, y = .data$Raw_Age),
                        color = "#f46f20", shape = 16, size = 1.1, alpha = 0.3) +
    ggplot2::geom_line(data = df,
                       ggplot2::aes(x = .data$Time, y = .data$loess_Age),
                       color = "black", linewidth = 1, alpha = 1.5) +
    ggplot2::geom_rect(
      data = segments,
      ggplot2::aes(xmin = .data$Start, xmax = .data$End,
                   ymin = .data$Segment_Mean - .data$Total_uncertainty,
                   ymax = .data$Segment_Mean + .data$Total_uncertainty),
      fill = final_color, alpha = 0.15) +
    ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = .data$Start, xend = .data$End,
                   y = .data$Segment_Mean, yend = .data$Segment_Mean),
      linewidth = 1.1, color = final_color) +
    ggplot2::ggtitle(title) +
    ggplot2::labs(x = "Time (s)", y = "Age (Ma)") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line        = ggplot2::element_line(color = "black"),
      axis.ticks       = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = ggplot2::unit(-0.3, "cm"),
      plot.title       = ggplot2::element_text(size = 10),
      axis.text        = ggplot2::element_text(size = 8),
      axis.title       = ggplot2::element_text(size = 10)
    )

  return(p)
}

#' Save depth profile plots to PDF
#'
#' @param plots List of ggplot objects
#' @param output_dir Directory to save PDFs
#' @param prefix Filename prefix (e.g., "Plot")
#' @keywords internal
save_plots <- function(plots, output_dir = ".", prefix = "Plot") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  for (i in seq_along(plots)) {
    label <- if (!is.null(plots[[i]]$label)) plots[[i]]$label else i
    fname <- file.path(output_dir, paste0(prefix, "_", label, ".pdf"))
    ggplot2::ggsave(
      filename = fname, plot = plots[[i]],
      device = "pdf", width = 2.6 * 25.4, height = (36 / 21) * 25.4,
      units = "mm"
    )
  }
}

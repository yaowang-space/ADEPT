# ADEPT: Automated Depth Profiling Technique
# adept.R — Main exported function
#
# Wang et al. (2026) method for automated identification of age plateaus
# from LA-ICP-MS depth profiling data.

#' Automated Depth Profiling Technique
#'
#' Quantitatively identifies and extracts age plateaus from LA-ICP-MS
#' U-Pb depth profiling data. Supports three input formats: direct ages
#' (Age68/Age75/Age76), raw isotope counts (Pb206/Pb207/U238), or
#' isotopic ratios (Pb206_U238/Pb207_U235/Pb207_Pb206).
#'
#' @param file_path Path to the input Excel file (.xlsx)
#' @param chunk_size Number of rows per processing chunk (default 411)
#' @param lower_ablation_time Minimum effective ablation time in seconds (default 29)
#' @param upper_ablation_time Maximum effective ablation time in seconds (default 58)
#' @param max_age_limit Maximum valid age in Ma (default 4540)
#' @param min_age_limit Minimum valid age in Ma (default 0)
#' @param min_plateau_resolution Minimum plateau duration in seconds.
#'   \code{NULL} or values < 5 default to 5 seconds.
#' @param variance_threshold Maximum allowed intra-plateau variance (default 0.1192)
#' @param filter_direction \code{"Forward"} keeps ascending age sequences;
#'   \code{"Reverse"} keeps descending age sequences.
#' @param mcmc Logical. Run Bayesian MCMC posterior analysis? (default \code{FALSE})
#' @param plot Logical. Generate and save depth profile plots? (default \code{TRUE})
#' @param output_path Output Excel file path. \code{NULL} auto-generates from input name.
#' @param plot_dir Directory for plot PDFs. \code{NULL} uses current directory.
#'
#' @return Invisibly returns a list with \code{summary}, \code{full}, and \code{plots}.
#'   Also writes a two-sheet Excel file and (optionally) PDF plots to disk.
#'
#' @export
#' @importFrom stats loess predict var sd residuals lm coef complete.cases na.omit
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' result <- adept("Input.xlsx")
#' result <- adept("Input.xlsx", filter_direction = "Reverse")
#' result <- adept("Input.xlsx", mcmc = TRUE)  # enable MCMC
#' }
adept <- function(
    file_path,
    chunk_size                = 411,
    lower_ablation_time       = 29,
    upper_ablation_time       = 58,
    max_age_limit             = 4540,
    min_age_limit             = 0,
    min_plateau_resolution    = NULL,
    variance_threshold        = 0.1192,
    filter_direction          = c("Forward", "Reverse"),
    mcmc                      = FALSE,
    plot                      = TRUE,
    output_path               = NULL,
    plot_dir                  = NULL
) {
  filter_direction <- match.arg(filter_direction)

  # ---- Read input ----
  data_list <- read_input(file_path)
  sheet_names <- names(data_list)

  # ---- Initialize ----
  output <- list()
  extra_mean_cols <- character(0)
  plots <- list()

  # ---- Auto output path ----
  if (is.null(output_path)) {
    base <- tools::file_path_sans_ext(basename(file_path))
    output_path <- file.path(dirname(file_path), paste0(base, "_Output.xlsx"))
  }
  if (is.null(plot_dir)) {
    plot_dir <- dirname(file_path)
  }

  # ---- Main processing loop ----
  sheet_idx <- 0
  for (sheet_name in sheet_names) {
    sheet_idx <- sheet_idx + 1
    segment_data_raw <- data_list[[sheet_name]]
    n_groups <- ceiling(nrow(segment_data_raw) / chunk_size)
    group_counter <- 1

    message(sprintf("Processing: '%s' (%d zircon(s))",
                    sheet_name, n_groups))

    for (i in seq_len(n_groups)) {
      start_row <- (i - 1) * chunk_size + 1
      end_row   <- min(i * chunk_size, nrow(segment_data_raw))

      # Parse data format
      parsed <- parse_segment_data(segment_data_raw, start_row, end_row)
      segment.data <- parsed$data
      extra_names  <- parsed$extra_names
      extra_raw    <- parsed$extra_data

      # Track rows for extra column alignment
      segment.data$.ROWID. <- seq_len(nrow(segment.data))
      if (length(extra_names) > 0 && !is.null(extra_raw)) {
        extra_raw$.ROWID. <- seq_len(nrow(extra_raw))
        extra_mean_cols <- paste0(extra_names, "_Mean")
      }

      # Check for missing ages
      if (all(is.na(segment.data$Age68)) ||
          all(is.na(segment.data$Age75)) ||
          all(is.na(segment.data$Age76))) {
        output[[length(output) + 1]] <- c(
          segment.data[1, 1], group_counter, 0,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, 0,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA
        )
        group_counter <- group_counter + 1
        next
      }

      # Numeric conversion & NA removal
      segment.data[, 2:ncol(segment.data)] <- lapply(
        segment.data[, 2:ncol(segment.data), drop = FALSE], as.numeric)
      segment.data <- na.omit(segment.data)

      # Ablation time filter
      subset_data <- subset(segment.data,
                            Time >= lower_ablation_time &
                            Time <= upper_ablation_time)
      subset_data$Raw_Age <- ifelse(subset_data$Age68 < 1000,
                                    subset_data$Age68,
                                    ifelse(subset_data$Age76 > 1000,
                                           subset_data$Age76, NA))

      if (all(is.na(subset_data$Age68)) ||
          all(is.na(subset_data$Age75)) ||
          all(is.na(subset_data$Age76))) {
        output[[length(output) + 1]] <- c(
          segment.data[1, 1], group_counter, nrow(subset_data),
          NA, NA, NA, NA, NA, NA, NA, NA, NA, 0,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA
        )
        group_counter <- group_counter + 1
        next
      }

      # ---- Preprocessing ----
      subset_data$Age68 <- arima_outlier(subset_data$Age68)
      subset_data$Age75 <- arima_outlier(subset_data$Age75)
      subset_data$Age76 <- arima_outlier(subset_data$Age76)
      subset_data <- discordance_filter(subset_data)
      subset_data$subset_Age68 <- mean_fill(subset_data$subset_Age68, subset_data$Age68)
      subset_data$subset_Age76 <- mean_fill(subset_data$subset_Age76, subset_data$Age76)
      subset_data$subset_Age <- ifelse(subset_data$subset_Age68 < 1000,
                                       subset_data$subset_Age68,
                                       ifelse(subset_data$subset_Age76 > 1000,
                                              subset_data$subset_Age76, NA))
      if (sum(!is.na(subset_data$subset_Age)) < 10) {
        output[[length(output) + 1]] <- c(
          segment.data[1, 1], group_counter, nrow(subset_data),
          NA, NA, NA, NA, NA, NA, NA, NA, NA, 0,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA
        )
        group_counter <- group_counter + 1
        next
      }

      subset_data <- subset_data[complete.cases(subset_data$subset_Age), ]
      subset_data$Row_Number <- seq_len(nrow(subset_data))

      # Merge extra columns (convert to numeric, non-numeric → NA)
      if (length(extra_names) > 0 && !is.null(extra_raw)) {
        matched_rows <- match(subset_data$.ROWID., extra_raw$.ROWID.)
        for (col in extra_names) {
          subset_data[[col]] <- suppressWarnings(as.numeric(as.character(extra_raw[[col]][matched_rows])))
        }
      }

      # ---- LOESS smoothing ----
      subset_data <- loess_segment(subset_data, span = 0.15)

      if (sum(!is.na(subset_data$standardized_loess)) < 10) {
        output[[length(output) + 1]] <- c(
          segment.data[1, 1], group_counter, nrow(subset_data),
          NA, NA, NA, NA, NA, NA, NA, NA, NA, 0,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
          NA, NA, NA
        )
        group_counter <- group_counter + 1
        next
      }

      # ---- PELT segmentation ----
      pelt <- pelt_segmentation(subset_data$standardized_loess, nrow(subset_data))
      changepoints <- pelt$changepoints
      subset_data$Change_loess <- ifelse(subset_data$Row_Number %in% changepoints,
                                         "YES", "NO")

      seg_result <- build_segments(subset_data, changepoints)
      segments  <- seg_result$segments
      subset_data <- seg_result$df
      seg_starts <- seg_result$segment_starts
      seg_ends   <- seg_result$segment_ends

      # ---- Plateau statistics ----
      segments <- calc_slopes(subset_data, segments, seg_starts, seg_ends)
      segments$Variance <- calc_variance(subset_data, seg_starts, seg_ends)
      segments <- calc_uncertainty(segments, subset_data, seg_starts, seg_ends)

      # Extra column means & age-specific means
      segments <- calc_extra_means(segments, subset_data, extra_names)
      segments <- calc_age_means(segments, subset_data)
      segments <- calc_concordance(segments)

      # ---- Filtering ----
      segments <- apply_filters(segments,
                                min_age       = min_age_limit,
                                max_age       = max_age_limit,
                                var_threshold = variance_threshold,
                                min_res       = min_plateau_resolution,
                                direction     = filter_direction)

      # ---- Progress ----
      n_confirmed <- sum(!is.na(segments$Filter_4))
      message(sprintf("[%d/%d] %s — %d plateau(s) confirmed (%s)",
                      i, n_groups,
                      as.character(subset_data[1, "Analysis"]),
                      n_confirmed, filter_direction))

      # ---- MCMC (optional) ----
      if (mcmc) {
        segments <- run_mcmc(subset_data, segments)
      }

      # ---- Plot ----
      if (plot) {
        p <- plot_depth_profile(subset_data, segments)
        p$label <- paste0(group_counter, "_", subset_data[1, "Analysis"])
        plots[[length(plots) + 1]] <- p
      }

      # ---- Build output row ----
      analysis_name <- as.character(subset_data[1, "Analysis"])

      # Extra column values
      extra_vals <- if (length(extra_mean_cols) > 0) {
        unname(sapply(extra_mean_cols, function(nm) {
          if (nm %in% colnames(segments)) segments[[nm]][1] else NA
        }))
      }

      for (k in seq_len(nrow(segments))) {
        # Collect extra values per row
        extra_vals_row <- if (length(extra_mean_cols) > 0) {
          unname(sapply(extra_mean_cols, function(nm) {
            if (nm %in% colnames(segments)) segments[[nm]][k] else NA
          }))
        }

        if (mcmc) {
          output[[length(output) + 1]] <- c(
            analysis_name, group_counter, nrow(subset_data),
            segments$Start[k], segments$End[k],
            segments$Time_step[k], segments$Max_step[k], segments$Min_step[k],
            segments$standardized_Mean[k], segments$Segment_Mean[k],
            segments$standardized_slope[k], segments$loess_slope[k],
            segments$Variance[k],
            segments$Calibration_uncertainty[k],
            segments$Plateau_uncertainty[k],
            segments$Total_uncertainty[k],
            segments$Number[k],
            segments$Filter_1[k], segments$Filter_2[k],
            segments$Filter_3[k], segments$Filter_4[k],
            segments$Final_Serial_Number[k],
            segments$Plateau_Numbers[k],
            segments$mean[k], segments$lower[k], segments$upper[k],
            segments$sigma[k],
            segments$Rhat[k], segments$n.eff[k], segments$mcp_step[k],
            segments$Filter_mcp[k], segments$Final_step_Numbers[k],
            segments$Final_Age[k], segments$Final_total_uncertainty[k],
            segments$Concordance[k],
            segments$Age68_Mean[k], segments$Age68_Total_uncertainty[k],
            segments$Age75_Mean[k], segments$Age75_Total_uncertainty[k],
            segments$Age76_Mean[k], segments$Age76_Total_uncertainty[k],
            extra_vals_row
          )
        } else {
          output[[length(output) + 1]] <- c(
            analysis_name, group_counter, nrow(subset_data),
            segments$Start[k], segments$End[k],
            segments$Time_step[k], segments$Max_step[k], segments$Min_step[k],
            segments$standardized_Mean[k], segments$Segment_Mean[k],
            segments$Variance[k],
            segments$Calibration_uncertainty[k],
            segments$Plateau_uncertainty[k],
            segments$Total_uncertainty[k],
            segments$Filter_1[k], segments$Filter_2[k],
            segments$Filter_3[k], segments$Filter_4[k],
            segments$Plateau_Numbers[k], segments$Number[k],
            segments$Final_step_Numbers[k], segments$Final_Serial_Number[k],
            segments$Filter_4[k], segments$Final_total_uncertainty[k],
            segments$Concordance[k],
            segments$Age68_Mean[k], segments$Age68_Total_uncertainty[k],
            segments$Age75_Mean[k], segments$Age75_Total_uncertainty[k],
            segments$Age76_Mean[k], segments$Age76_Total_uncertainty[k],
            extra_vals_row
          )
        }
      }

      group_counter <- group_counter + 1
    }
  }

  # ---- Pad all rows to uniform width, then rbind ----
  if (length(output) > 0) {
    max_width <- max(sapply(output, length))
    output <- do.call(rbind, lapply(output, function(row) {
      if (length(row) < max_width) c(row, rep(NA, max_width - length(row)))
      else row
    }))
  }

  # ---- Set column names ----
  if (mcmc) {
    colnames(output) <- c(
      "Analysis", "Group", "Points",
      "Start", "End",
      "Integration time", "Max integration time", "Min integration time",
      "Standardized mean", "Segmentation mean",
      "Standardized slope", "Loess slope",
      "Variance",
      "Calibration uncertainty", "Plateau uncertainty", "Total uncertainty",
      "Plateau serial numbers",
      "Filter 1", "Filter 2", "Filter 3", "Filter 4",
      "Final serial number", "Total integration time Numbers",
      "MCMC mean", "MCMC lower", "MCMC upper", "MCMC sigma",
      "Rhat", "MCMC n.eff", "MCMC integration time",
      "Filter MCMC", "Final integration time numbers",
      "Final age (Ma)", "Final total uncertainty (Ma)", "Concordance (%)",
      "Pb206/U238 age mean (Ma)", "Pb206/U238 total uncertainty (Ma)",
      "Pb207/U235 age mean (Ma)", "Pb207/U235 total uncertainty (Ma)",
      "Pb207/Pb206 age mean (Ma)", "Pb207/Pb206 total uncertainty (Ma)",
      extra_mean_cols
    )
  } else {
    colnames(output) <- c(
      "Analysis", "Group", "Points",
      "Start", "End",
      "Integration time", "Max integration time", "Min integration time",
      "Standardized mean", "Segmentation mean",
      "Variance",
      "Calibration uncertainty", "Plateau uncertainty", "Total uncertainty",
      "Filter 1", "Filter 2", "Filter 3", "Filter 4",
      "Total integration time numbers", "Plateau serial numbers",
      "Final integration time numbers", "Final Serial Number",
      "Final age (Ma)", "Final total uncertainty (Ma)", "Concordance (%)",
      "Pb206/U238 age mean (Ma)", "Pb206/U238 total uncertainty (Ma)",
      "Pb207/U235 age mean (Ma)", "Pb207/U235 total uncertainty (Ma)",
      "Pb207/Pb206 age mean (Ma)", "Pb207/Pb206 total uncertainty (Ma)",
      extra_mean_cols
    )
  }

  # ---- Convert to data.frame ----
  output <- as.data.frame(output, stringsAsFactors = FALSE)

  # ---- Build simplified summary ----
  simplified_cols <- c("Analysis", "Group", "Points",
                       "Final serial number", "Final Serial Number",
                       "Integration time",
                       "Final age (Ma)", "Final total uncertainty (Ma)",
                       "Concordance (%)",
                       "Pb206/U238 age mean (Ma)", "Pb206/U238 total uncertainty (Ma)",
                       "Pb207/U235 age mean (Ma)", "Pb207/U235 total uncertainty (Ma)",
                       "Pb207/Pb206 age mean (Ma)", "Pb207/Pb206 total uncertainty (Ma)")
  if (length(extra_mean_cols) > 0) {
    simplified_cols <- c(simplified_cols, extra_mean_cols)
  }
  simplified_cols <- intersect(simplified_cols, colnames(output))
  output_simple <- output[, simplified_cols, drop = FALSE]
  output_simple <- output_simple[!is.na(output_simple[, "Final age (Ma)"]), , drop = FALSE]

  # ---- Export ----
  writexl::write_xlsx(
    list(Summary = output_simple, Full_Results = output),
    output_path
  )

  if (plot && length(plots) > 0) {
    save_plots(plots, plot_dir)
  }

  invisible(list(summary = output_simple, full = output, plots = plots))
}

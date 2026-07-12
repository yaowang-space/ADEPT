# ADEPT: Automated Depth Profiling Technique
# processing.R — Core data processing pipeline
#
# Internal functions for reading, preprocessing, LOESS smoothing,
# PELT segmentation, plateau statistics, and uncertainty calculation.

#' Read all sheets from an Excel input file
#'
#' @param file_path Path to the input Excel file
#' @return A named list of data.frames, one per sheet
#' @keywords internal
read_input <- function(file_path) {
  sheet_names <- readxl::excel_sheets(file_path)
  data_list <- list()
  for (i in seq_along(sheet_names)) {
    data_list[[sheet_names[i]]] <- openxlsx::read.xlsx(file_path, sheet = sheet_names[i])
  }
  return(data_list)
}

#' Parse raw segment data and extract age columns
#'
#' Detects the input format (Age68, Pb206, or Pb206_U238) and creates
#' a standardized data.frame with Analysis, Time, Age68, Age75, Age76.
#'
#' @param raw Raw data.frame from one sheet
#' @param start_row Starting row index for this chunk
#' @param end_row Ending row index for this chunk
#' @return A list with `data` (standardized data.frame) and `extra_names` (character vector of additional numeric column names)
#' @keywords internal
parse_segment_data <- function(raw, start_row, end_row) {
  colnames(raw) <- gsub(" |/", "_", colnames(raw))

  if ("Age68" %in% colnames(raw) &&
      !all(is.na(raw$Age68[start_row:end_row]))) {
    segment.data <- data.frame(
      Analysis = raw$Analysis[start_row:end_row],
      Time     = raw$Time[start_row:end_row],
      Age68    = raw$Age68[start_row:end_row],
      Age75    = raw$Age75[start_row:end_row],
      Age76    = raw$Age76[start_row:end_row],
      stringsAsFactors = FALSE
    )
    segment.data[, 2:5] <- lapply(segment.data[, 2:5], as.numeric)

  } else if ("Pb206" %in% colnames(raw) &&
             !all(is.na(raw$Pb206[start_row:end_row]))) {
    segment.data <- data.frame(
      Analysis = raw$Analysis[start_row:end_row],
      Time     = raw$Time[start_row:end_row],
      Pb206    = raw$Pb206[start_row:end_row],
      Pb207    = raw$Pb207[start_row:end_row],
      U238     = raw$U238[start_row:end_row],
      stringsAsFactors = FALSE
    )
    segment.data[, 2:5] <- lapply(segment.data[, 2:5], as.numeric)
    segment.data$U235        <- segment.data$U238 / 137.88
    segment.data$Pb206U238   <- segment.data$Pb206 / segment.data$U238
    segment.data$Pb207U235   <- segment.data$Pb207 / segment.data$U235
    segment.data$Pb207Pb206  <- segment.data$Pb207 / segment.data$Pb206
    segment.data$Age68 <- sapply(segment.data$Pb206U238,  calc_age_ratio, method = "U238-Pb206")
    segment.data$Age75 <- sapply(segment.data$Pb207U235,  calc_age_ratio, method = "U235-Pb207")
    segment.data$Age76 <- sapply(segment.data$Pb207Pb206, calc_age_ratio, method = "Pb207-Pb206")

  } else if ("Pb206_U238" %in% colnames(raw) &&
             !all(is.na(raw$Pb206_U238[start_row:end_row]))) {
    segment.data <- data.frame(
      Analysis   = raw$Analysis[start_row:end_row],
      Time       = raw$Time[start_row:end_row],
      Pb206U238  = raw$Pb206_U238[start_row:end_row],
      Pb207U235  = raw$Pb207_U235[start_row:end_row],
      Pb207Pb206 = raw$Pb207_Pb206[start_row:end_row],
      stringsAsFactors = FALSE
    )
    segment.data[, 2:5] <- lapply(segment.data[, 2:5], as.numeric)
    segment.data$Age68 <- sapply(segment.data$Pb206U238,  calc_age_ratio, method = "U238-Pb206")
    segment.data$Age75 <- sapply(segment.data$Pb207U235,  calc_age_ratio, method = "U235-Pb207")
    segment.data$Age76 <- sapply(segment.data$Pb207Pb206, calc_age_ratio, method = "Pb207-Pb206")
  } else {
    stop("Unrecognized data format. Expected Age68, Pb206, or Pb206_U238 columns.")
  }

  # Extract additional numeric columns
  extra_names <- setdiff(colnames(raw), colnames(segment.data))
  extra_data <- NULL
  if (length(extra_names) > 0) {
    extra_data <- raw[start_row:end_row, extra_names, drop = FALSE]
    extra_is_num <- sapply(extra_names, function(cn) {
      col_vals <- extra_data[[cn]]
      if (is.numeric(col_vals)) return(TRUE)
      converted <- suppressWarnings(as.numeric(col_vals))
      all(is.na(col_vals) | !is.na(converted))
    })
    extra_names <- extra_names[extra_is_num]
  }

  list(data = segment.data, extra_names = extra_names, extra_data = extra_data)
}

#' Calculate U-Pb age from isotopic ratio
#'
#' @param ratio Isotopic ratio value
#' @param method One of "U238-Pb206", "U235-Pb207", "Pb207-Pb206"
#' @return Age in Ma, or NA
#' @keywords internal
calc_age_ratio <- function(ratio, method) {
  if (is.finite(ratio) && ratio > 0) {
    input <- c(ratio, 0)
    result <- IsoplotR::age(input, method = method, exterr = FALSE)
    return(result[1])
  } else {
    return(NA)
  }
}

#' ARIMA-based outlier detection
#'
#' Fits an automatic ARIMA model, computes residuals, and marks values
#' beyond 2 SD as outliers (NA).
#'
#' @param series Numeric vector
#' @return Numeric vector with outliers replaced by NA
#' @keywords internal
arima_outlier <- function(series) {
  arima_model <- forecast::auto.arima(series)
  residuals   <- stats::residuals(arima_model)
  outlier_idx <- which(abs(residuals) > 2 * sd(residuals))
  series[outlier_idx] <- NA
  return(series)
}

#' Discordance filter for U-Pb ages
#'
#' For ages < 1000 Ma: marks Age68 as NA if |Age68 - Age75| / Age75 > 0.1
#'
#' @param df data.frame with Age68, Age75, Age76 columns
#' @return Modified data.frame with subset_Age68 and subset_Age columns
#' @keywords internal
discordance_filter <- function(df) {
  df$Age <- ifelse(df$Age68 < 1000, df$Age68,
                   ifelse(df$Age76 > 1000, df$Age76, NA))
  df$subset_Age68 <- ifelse(abs(df$Age68 - df$Age75) / df$Age75 <= 0.1,
                            df$Age68, NA)
  df$subset_Age76 <- df$Age76
  return(df)
}

#' Sliding window mean fill for NA values
#'
#' @param na_series Series with NA values to fill
#' @param original_series Original complete series for reference
#' @param window_size Half-window size (default 5)
#' @return Series with NAs filled by local mean
#' @keywords internal
mean_fill <- function(na_series, original_series, window_size = 5) {
  calc_local_mean <- function(series, index) {
    start_idx <- max(1, index - window_size)
    end_idx   <- min(length(series), index + window_size)
    mean(series[start_idx:end_idx], na.rm = TRUE)
  }
  na_idx <- which(is.na(na_series))
  for (idx in na_idx) {
    na_series[idx] <- calc_local_mean(original_series, idx)
  }
  return(na_series)
}

#' LOESS smoothing and standardization
#'
#' @param df data.frame with subset_Age and Time columns
#' @param span LOESS span parameter (default 0.15)
#' @return data.frame with added loess_Age, standardized_loess, standardized_Age columns
#' @keywords internal
loess_segment <- function(df, span = 0.15) {
  loess_model <- stats::loess(subset_Age ~ Time, data = df, span = span)
  df$loess_Age <- stats::predict(loess_model)
  df <- df[complete.cases(df$loess_Age), ]

  df$log_loess <- df$loess_Age
  df$log_Age   <- df$Age

  minloess <- min(df$log_loess)
  maxloess <- max(df$log_loess)

  df$standardized_loess <- (df$log_loess - minloess) / (maxloess - minloess)
  df$standardized_Age   <- (df$log_Age   - minloess) / (maxloess - minloess)
  df <- df[complete.cases(df$standardized_loess), ]

  attr(df, "minloess") <- minloess
  attr(df, "maxloess") <- maxloess
  return(df)
}

#' PELT changepoint detection
#'
#' @param series Numeric vector of standardized LOESS values
#' @param n Number of observations (for SAIC penalty)
#' @return A list with `changepoints` (integer vector) and `SAIC` (penalty value)
#' @keywords internal
pelt_segmentation <- function(series, n) {
  AIC_result    <- changepoint::cpt.mean(series, method = "PELT",
                                         penalty = "AIC", minseglen = 1)
  aic_pen_value <- attr(AIC_result, "pen.value")
  SAIC <- (aic_pen_value / 100) * log(n)

  cpt_result <- changepoint::cpt.mean(series, method = "PELT",
                                      penalty = "Manual", pen.value = SAIC,
                                      minseglen = 1)
  list(changepoints = cpt_result@cpts, SAIC = SAIC)
}

#' Build plateau segment definitions from changepoints
#'
#' @param df data.frame with Time and standardized_loess columns
#' @param changepoints Integer vector of changepoint row indices
#' @return data.frame of segment start/end times and properties
#' @keywords internal
build_segments <- function(df, changepoints) {
  segment_starts <- c(1, utils::head(changepoints, -1) + 1)
  segment_ends   <- changepoints

  segments <- data.frame(
    Start = df$Time[segment_starts],
    End   = df$Time[segment_ends],
    stringsAsFactors = FALSE
  )

  # Standardized mean per segment
  segments$standardized_Mean <- sapply(seq_len(nrow(segments)), function(i) {
    mean(df$standardized_loess[segments$Start[i] <= df$Time &
                               df$Time <= segments$End[i]])
  })

  # Intra-segment normalization for variance
  normalize_seg <- function(series, s, e) {
    seg <- series[s:e]
    (seg - min(seg)) / (max(seg) - min(seg))
  }

  normalized <- lapply(seq_along(segment_starts), function(i) {
    normalize_seg(df$loess_Age, segment_starts[i], segment_ends[i])
  })
  df$IS_loess <- unlist(normalized)

  # Restore from standardization
  minloess <- attr(df, "minloess")
  maxloess <- attr(df, "maxloess")
  df$Restored_loess <- df$standardized_loess * (maxloess - minloess) + minloess
  df$Restored_age   <- df$standardized_Age   * (maxloess - minloess) + minloess

  # Segment mean (original scale)
  segments$Segment_Mean <- sapply(seq_len(nrow(segments)), function(i) {
    mean(df$Restored_loess[segments$Start[i] <= df$Time &
                           df$Time <= segments$End[i]])
  })

  segments$Number <- seq_len(nrow(segments))
  segments$Time_step <- segments$End - segments$Start
  segments$Max_step  <- max(segments$Time_step)
  segments$Min_step  <- min(segments$Time_step)

  list(segments = segments, df = df, segment_starts = segment_starts,
       segment_ends = segment_ends)
}

#' Calculate slope and intercept for each plateau segment
#'
#' @param df data.frame with Time and age columns
#' @param segments Segment definitions data.frame
#' @param seg_starts Integer vector of segment start row indices
#' @param seg_ends Integer vector of segment end row indices
#' @return Updated segments data.frame with slope/intercept columns
#' @keywords internal
calc_slopes <- function(df, segments, seg_starts, seg_ends) {
  calc_slope_int <- function(series, time, start, end) {
    model <- stats::lm(series[start:end] ~ time[start:end])
    stats::coef(model)
  }

  segments$standardized_intercept <- sapply(seq_along(seg_starts), function(i) {
    calc_slope_int(df$standardized_loess, df$Time, seg_starts[i], seg_ends[i])[1]
  })
  segments$standardized_slope <- sapply(seq_along(seg_starts), function(i) {
    calc_slope_int(df$standardized_loess, df$Time, seg_starts[i], seg_ends[i])[2]
  })
  segments$loess_intercept <- sapply(seq_along(seg_starts), function(i) {
    calc_slope_int(df$loess_Age, df$Time, seg_starts[i], seg_ends[i])[1]
  })
  segments$loess_slope <- sapply(seq_along(seg_starts), function(i) {
    calc_slope_int(df$loess_Age, df$Time, seg_starts[i], seg_ends[i])[2]
  })

  return(segments)
}

#' Calculate plateau variance
#'
#' @param df data.frame with IS_loess column
#' @param seg_starts Integer vector of segment start row indices
#' @param seg_ends Integer vector of segment end row indices
#' @return Numeric vector of variance per segment
#' @keywords internal
calc_variance <- function(df, seg_starts, seg_ends) {
  sapply(seq_along(seg_starts), function(i) {
    stats::var(df$IS_loess[seg_starts[i]:seg_ends[i]], na.rm = TRUE)
  })
}

#' Calculate plateau uncertainties
#'
#' Calibration uncertainty = Segment_Mean * 0.03
#' Plateau uncertainty = sd(Raw_Age) / sqrt(n) within segment
#' Total = sqrt(cal^2 + plateau^2)
#'
#' @param segments Segment definitions data.frame
#' @param df data.frame with Raw_Age column
#' @param seg_starts Segment start row indices
#' @param seg_ends Segment end row indices
#' @return Updated segments data.frame with uncertainty columns
#' @keywords internal
calc_uncertainty <- function(segments, df, seg_starts, seg_ends) {
  segments$Calibration_uncertainty <- segments$Segment_Mean * 0.03

  segments$Plateau_uncertainty <- sapply(seq_along(seg_starts), function(i) {
    seg <- df$Raw_Age[seg_starts[i]:seg_ends[i]]
    seg <- seg[!is.na(seg)]
    n  <- length(seg)
    if (n < 2) return(NA)
    stats::sd(seg) / sqrt(n)
  })

  segments$Total_uncertainty <- sqrt(
    segments$Calibration_uncertainty^2 + segments$Plateau_uncertainty^2
  )

  return(segments)
}

#' Calculate extra (non-age) column means per plateau segment
#'
#' @param segments Segment definitions data.frame
#' @param df data.frame with extra columns
#' @param extra_names Character vector of extra column names
#' @return Updated segments data.frame with *_Mean columns
#' @keywords internal
calc_extra_means <- function(segments, df, extra_names) {
  if (length(extra_names) == 0) return(segments)

  safe_mean <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else mean(x)
  }

  for (col_name in extra_names) {
    mean_name <- paste0(col_name, "_Mean")
    if (col_name %in% colnames(df)) {
      segments[[mean_name]] <- sapply(seq_len(nrow(segments)), function(k) {
        vals <- df[[col_name]][segments$Start[k] <= df$Time &
                               df$Time <= segments$End[k]]
        safe_mean(vals)
      })
    } else {
      segments[[mean_name]] <- NA
    }
  }
  return(segments)
}

#' Calculate Age68, Age75, Age76 means and total uncertainties per plateau
#'
#' @param segments Segment definitions data.frame
#' @param df data.frame with Age68, Age75, Age76 columns
#' @return Updated segments data.frame with age mean/uncertainty columns
#' @keywords internal
calc_age_means <- function(segments, df) {
  safe_mean <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else mean(x)
  }

  for (age_col in c("Age68", "Age75", "Age76")) {
    if (!age_col %in% colnames(df)) next

    mean_name <- paste0(age_col, "_Mean")
    unc_name  <- paste0(age_col, "_Total_uncertainty")

    segments[[mean_name]] <- sapply(seq_len(nrow(segments)), function(k) {
      vals <- df[[age_col]][segments$Start[k] <= df$Time &
                            df$Time <= segments$End[k]]
      safe_mean(vals)
    })

    segments[[unc_name]] <- sapply(seq_len(nrow(segments)), function(k) {
      vals <- df[[age_col]][segments$Start[k] <= df$Time &
                            df$Time <= segments$End[k]]
      vals <- vals[!is.na(vals)]
      n <- length(vals)
      if (n < 2) return(NA)
      plateau_unc <- stats::sd(vals) / sqrt(n)
      cal_unc <- mean(vals) * 0.03
      sqrt(cal_unc^2 + plateau_unc^2)
    })
  }
  return(segments)
}

#' Calculate U-Pb concordance (Age68 / Age75 * 100) per plateau
#'
#' @param segments Segment definitions data.frame
#' @return Updated segments data.frame with Concordance column
#' @keywords internal
calc_concordance <- function(segments) {
  segments$Concordance <- (segments$Age68_Mean / segments$Age75_Mean) * 100
  return(segments)
}

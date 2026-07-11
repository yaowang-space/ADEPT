# ADEPT: Automated Depth Profiling Technique
# filtering.R — Four-step plateau filtering
#
# Internal functions implementing the multi-step filtering cascade.

#' Step 1: Filter by age range
#'
#' Keeps plateaus whose Segment_Mean falls within [min_age, max_age].
#'
#' @param segments Segment data.frame
#' @param min_age Minimum valid age (Ma)
#' @param max_age Maximum valid age (Ma)
#' @return Numeric vector (Filter_1): age value or NA
#' @keywords internal
filter_age_range <- function(segments, min_age = 0, max_age = 4540) {
  ifelse(as.numeric(segments$Segment_Mean) <= max_age &
         as.numeric(segments$Segment_Mean) >= min_age,
         segments$Segment_Mean, NA)
}

#' Step 2: Filter by plateau variance
#'
#' @param segments Segment data.frame
#' @param filter_1 Result from filter_age_range()
#' @param threshold Variance threshold (default 0.1192)
#' @return Numeric vector (Filter_2): age value or NA
#' @keywords internal
filter_variance <- function(segments, filter_1, threshold = 0.1192) {
  ifelse(as.numeric(segments$Variance) <= threshold, filter_1, NA)
}

#' Step 3: Filter by minimum plateau resolution
#'
#' Removes plateaus shorter than the minimum resolution (default >= 5s).
#'
#' @param segments Segment data.frame
#' @param filter_2 Result from filter_variance()
#' @param min_res Minimum plateau duration in seconds (NULL or < 5 defaults to 5)
#' @return Numeric vector (Filter_3): age value or NA
#' @keywords internal
filter_resolution <- function(segments, filter_2, min_res = NULL) {
  if (is.null(min_res) || is.na(min_res) || min_res < 5) min_res <- 5
  result <- filter_2
  result[-1] <- ifelse(as.numeric(segments$Time_step[-1]) >= min_res,
                       filter_2[-1], NA)
  result
}

#' Step 4: Directional plateau selection (Forward / Reverse)
#'
#' Forward: keep all ascending plateaus; if not ascending, keep the first
#'   (youngest) plateau. If volatile, additionally keep the most stable.
#' Reverse: keep all descending plateaus; if not descending, keep the first
#'   (youngest) plateau. If volatile, additionally keep the most stable.
#'
#' @param segments Segment data.frame
#' @param filter_3 Result from filter_resolution()
#' @param direction "Forward" (keep ascending) or "Reverse" (keep descending)
#' @return Numeric vector (Filter_4): final age value or NA
#' @keywords internal
filter_direction <- function(segments, filter_3, direction = c("Forward", "Reverse")) {
  direction <- match.arg(direction)
  NO_NULL_Age <- filter_3[!is.na(filter_3)]
  filter_4 <- rep(NA, nrow(segments))

  if (length(NO_NULL_Age) == 0) return(filter_4)

  is_ascending  <- all(diff(NO_NULL_Age, na.rm = TRUE) >= 0)
  is_descending <- all(diff(NO_NULL_Age, na.rm = TRUE) < 0)

  if (direction == "Forward") {
    if (is_ascending) {
      filter_4 <- filter_3
    } else {
      first_idx <- which(filter_3 == NO_NULL_Age[1])[1]
      filter_4[first_idx] <- NO_NULL_Age[1]
      if (!is_descending) {
        remaining <- which(!is.na(filter_3) & is.na(filter_4))
        if (length(remaining) > 0) {
          best_idx <- remaining[which.min(segments$Variance[remaining])]
          filter_4[best_idx] <- filter_3[best_idx]
        }
      }
    }
  } else {  # Reverse
    if (is_descending) {
      filter_4 <- filter_3
    } else {
      first_idx <- which(filter_3 == NO_NULL_Age[1])[1]
      filter_4[first_idx] <- NO_NULL_Age[1]
      if (!is_ascending) {
        remaining <- which(!is.na(filter_3) & is.na(filter_4))
        if (length(remaining) > 0) {
          best_idx <- remaining[which.min(segments$Variance[remaining])]
          filter_4[best_idx] <- filter_3[best_idx]
        }
      }
    }
  }

  return(filter_4)
}

#' Apply all four filtering steps
#'
#' @param segments Segment data.frame
#' @param min_age, max_age Age limits (Ma)
#' @param var_threshold Variance threshold
#' @param min_res Minimum plateau resolution (seconds)
#' @param direction "Forward" or "Reverse"
#' @return Updated segments data.frame with Filter_1..Filter_4 and derived columns
#' @keywords internal
apply_filters <- function(segments, min_age = 0, max_age = 4540,
                          var_threshold = 0.1192, min_res = NULL,
                          direction = c("Forward", "Reverse")) {
  direction <- match.arg(direction)

  segments$Filter_1 <- filter_age_range(segments, min_age, max_age)
  segments$Filter_2 <- filter_variance(segments, segments$Filter_1, var_threshold)
  segments$Filter_3 <- filter_resolution(segments, segments$Filter_2, min_res)
  segments$Filter_4 <- filter_direction(segments, segments$Filter_3, direction)

  segments$Final_total_uncertainty <- ifelse(!is.na(segments$Filter_4),
                                             segments$Total_uncertainty, NA)
  segments$Plateau_Numbers <- length(segments$Filter_4)
  segments$Final_Serial_Number <- NA
  non_na <- which(!is.na(segments$Filter_4))
  segments$Final_Serial_Number[non_na] <- seq_along(non_na)
  segments$Final_Age <- segments$Filter_4
  segments$Final_step_Numbers <- sum(!is.na(segments$Filter_4))

  return(segments)
}

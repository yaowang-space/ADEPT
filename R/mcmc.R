# ADEPT: Automated Depth Profiling Technique
# mcmc.R — Bayesian MCMC posterior analysis (optional)
#
# Uses the `mcp` package to fit a changepoint model and extract
# posterior estimates for each plateau segment.

#' Run Bayesian MCMC changepoint analysis
#'
#' Fits an intercept-only changepoint model using the mcp package
#' and returns posterior summaries merged with segment data.
#'
#' @param df data.frame with Time and standardized_loess columns
#' @param segments Segment definitions data.frame
#' @return Updated segments data.frame with MCMC mean, lower, upper, sigma, Rhat, n.eff, step columns
#' @keywords internal
run_mcmc <- function(df, segments) {
  if (!requireNamespace("mcp", quietly = TRUE)) {
    warning("Package 'mcp' is not installed. Skipping MCMC analysis.")
    return(segments)
  }

  n_segments <- nrow(segments)
  model <- list(standardized_loess ~ 1)
  if (n_segments > 1) {
    model <- c(model, replicate(n_segments - 1, ~ 1))
  }

  mcp_fit <- mcp::mcp(model, df, par_x = "Time")
  sum_mcp <- summary(mcp_fit)

  intercepts_info <- sum_mcp[grep("int", sum_mcp$name), ]
  cp_info         <- sum_mcp[grep("cp", sum_mcp$name), ]
  sigma_info      <- sum_mcp[grep("sigma", sum_mcp$name), ]

  intercepts_info <- data.frame(intercepts_info)
  intercepts_info$sigma   <- as.numeric(sigma_info$upper)
  intercepts_info$sd_mean <- as.numeric(intercepts_info$mean)
  intercepts_info$mean    <- as.numeric(intercepts_info$mean)
  intercepts_info$lower   <- as.numeric(intercepts_info$lower)
  intercepts_info$upper   <- as.numeric(intercepts_info$upper)

  segments <- cbind(segments, intercepts_info)

  mean_vals <- cp_info$mean
  mcp_step  <- numeric(length(mean_vals) + 1)

  if (sum(!is.na(mcp_step)) >= 3) {
    mcp_step[1] <- mean_vals[1] - df$Time[1]
    mcp_step[2:length(mean_vals)] <- diff(mean_vals)
    mcp_step[length(mcp_step)] <- df$Time[length(df$Time)] - mean_vals[length(mean_vals)]
  } else {
    mcp_step[1] <- mean_vals[1] - df$Time[1]
    mcp_step[length(mcp_step)] <- df$Time[length(df$Time)] - mean_vals[length(mean_vals)]
  }

  segments <- cbind(segments, data.frame(mcp_step = mcp_step))
  segments$Filter_mcp <- ifelse(!is.na(segments$Filter_4), segments$mean, NA)

  return(segments)
}

# ADEPT: Automated Depth Profiling Technique

An R package for automated identification and extraction of age plateaus from
LA-ICP-MS zircon U-Pb depth profiling data.

**Reference:** Wang et al. (2025) — *Automated Depth Profiling Technique for
quantitative identification of depth profiling age plateaus.*

## Installation

```r
# Install from GitHub
if (!require("remotes")) install.packages("remotes")
remotes::install_github("user/ADEPT")
```

### Dependencies

The package requires:

- **Required:** `stats`, `ggplot2`, `zoo`, `changepoint`, `readxl`, `forecast`, `openxlsx`, `writexl`, `dplyr`, `MASS`, `gridExtra`
- **Optional:** `IsoplotR` (for isotope ratio → age conversion), `mcp` (for Bayesian MCMC)

## Quick Start

```r
library(ADEPT)

# Basic usage with defaults (Forward direction, MCMC on)
result <- adept("Input.xlsx")

# Reverse direction, no MCMC (faster)
result <- adept("Input.xlsx", filter_direction = "Reverse", mcmc = FALSE)

# Access results
head(result$summary)   # Simplified summary (confirmed plateaus only)
head(result$full)      # Full results (all plateaus, all columns)
result$plots[[1]]      # First depth profile plot
```

## Main Function

```r
adept(
  file_path,                   # Path to input Excel file
  chunk_size             = 411,
  lower_ablation_time    = 29,     # seconds
  upper_ablation_time    = 58,     # seconds
  max_age_limit          = 4540,   # Ma
  min_age_limit          = 0,      # Ma
  min_plateau_resolution = NULL,   # NULL = auto (>= 5 s)
  variance_threshold     = 0.1192,
  filter_direction       = "Forward",  # "Forward" or "Reverse"
  mcmc                   = TRUE,       # Bayesian MCMC?
  plot                   = TRUE,       # Generate PDF plots?
  output_path            = NULL,       # NULL = auto-named
  plot_dir               = NULL        # NULL = same as input
)
```

## Parameter Details

| Parameter | Type | Default | Description |
|---|---|---|---|
| `file_path` | character | (required) | Path to input `.xlsx` file |
| `chunk_size` | numeric | 411 | Rows per processing chunk |
| `lower_ablation_time` | numeric | 29 | Min ablation time in seconds |
| `upper_ablation_time` | numeric | 58 | Max ablation time in seconds |
| `max_age_limit` | numeric | 4540 | Maximum valid age (Ma) |
| `min_age_limit` | numeric | 0 | Minimum valid age (Ma) |
| `min_plateau_resolution` | numeric / NULL | NULL | Min plateau duration (s). NULL defaults to 5 |
| `variance_threshold` | numeric | 0.1192 | Maximum plateau variance |
| `filter_direction` | character | "Forward" | "Forward" = keep ascending; "Reverse" = keep descending |
| `mcmc` | logical | TRUE | Run MCMC posterior analysis? |
| `plot` | logical | TRUE | Save depth profile PDFs? |
| `output_path` | character / NULL | NULL | Output Excel path. NULL = auto |
| `plot_dir` | character / NULL | NULL | Plot PDF directory. NULL = input dir |

## Input Data Format

The input Excel file must contain sheets with at least one of the following:

### Format 1: Direct ages
Columns: `Analysis`, `Time`, `Age68`, `Age75`, `Age76`

### Format 2: Raw isotope counts
Columns: `Analysis`, `Time`, `Pb206`, `Pb207`, `U238`

### Format 3: Isotopic ratios
Columns: `Analysis`, `Time`, `Pb206_U238`, `Pb207_U235`, `Pb207_Pb206`

Any additional numeric columns (e.g., trace elements) will be automatically
detected and their plateau means will be included in the output.

## Output

### Excel file (two sheets)

**Sheet 1 — Summary:** Confirmed age plateaus only.

| Column | Description |
|---|---|
| Analysis | Sample analysis ID |
| Group | Chunk group number |
| Points | Number of data points in chunk |
| Final serial number | Plateau rank (1, 2, ...) |
| Integration time | Plateau duration (s) |
| Final age (Ma) | Plateau age |
| Final total uncertainty (Ma) | Combined uncertainty |
| Concordance (%) | Age68/Age75 × 100 |
| Pb206/U238 age mean (Ma) | Mean 206Pb/238U age |
| Pb206/U238 total uncertainty (Ma) | Its uncertainty |
| Pb207/U235 age mean (Ma) | Mean 207Pb/235U age |
| Pb207/U235 total uncertainty (Ma) | Its uncertainty |
| Pb207/Pb206 age mean (Ma) | Mean 207Pb/206Pb age |
| Pb207/Pb206 total uncertainty (Ma) | Its uncertainty |
| `*_Mean` | Means of extra input columns |

**Sheet 2 — Full Results:** All plateaus with complete statistics including
MCMC posterior estimates (if `mcmc = TRUE`), slope/intercept, and filter flags.

### Return Value

```r
list(
  summary = data.frame,   # Sheet 1 content
  full    = data.frame,   # Sheet 2 content
  plots   = list()        # ggplot objects (if plot = TRUE)
)
```

## Processing Pipeline

```
Input Excel → Format Detection → ARIMA Outliers → Discordance Filter
→ Mean Fill → LOESS Smoothing → Standardization → PELT Segmentation
→ Plateau Statistics → 4-Step Filtering → [MCMC] → Output
```

## Citation

Wang et al. (2025). Automated Depth Profiling Technique for quantitative
identification of depth profiling age plateaus. *Journal of Geophysical
Research: Solid Earth.*

## License

MIT

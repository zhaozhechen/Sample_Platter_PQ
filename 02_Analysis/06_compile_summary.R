# Build the requested cross-site table and a GitHub-readable complete report.
# This script summarizes existing TE tables; it does not recalculate TE.

.libPaths(c(normalizePath("R_library"), .libPaths()))
source("01_Functions/TE_implementation.R")

site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))
audit <- data.table::fread("04_Results/Data_Audit/Input_data_summary.csv")
configured <- site_info[run == TRUE]
rows <- vector("list", nrow(configured))

for (i in seq_len(nrow(configured))) {
  spec <- configured[i]
  result_file <- file.path("04_Results", "Tables", paste0("TE_df_", spec$Analysis_ID, ".csv"))
  if (!file.exists(result_file)) stop("Missing result table: ", result_file)
  te_df <- data.table::fread(result_file)
  audit_row <- audit[Analysis_ID == spec$Analysis_ID]
  if (nrow(audit_row) != 1L) stop("Missing or duplicate audit row: ", spec$Analysis_ID)

  peak_index <- which.max(te_df$TE)
  significant <- is.finite(te_df$TE) & is.finite(te_df$TEcrit) & te_df$TE > te_df$TEcrit
  significant_count <- sum(significant)
  # User definition: total significant TE divided by the configured lag count.
  # There are 73 evaluated positions (0:72), but the denominator is 72 lag steps.
  cumulative_te <- sum(te_df$TE[significant], na.rm = TRUE) / as.integer(spec$max_lag_steps)
  notes <- spec$provenance_note
  if (!significant_count) notes <- paste(notes, "No TE value exceeded the shuffled threshold.")

  rows[[i]] <- data.table::data.table(
    Analysis_ID = spec$Analysis_ID,
    Site = spec$Site_Name,
    Variable_pair = paste0(spec$source_label, " -> ", spec$sink_label),
    Native_resolution_hours = audit_row$step_hours,
    Native_resolution_minutes = audit_row$step_hours * 60,
    Data_length_observations = audit_row$rows,
    Data_duration_days = audit_row$duration_days,
    Start_date = audit_row$start,
    End_date = audit_row$end,
    Peak_TE = te_df$TE[peak_index],
    Peak_lag_steps = te_df$Lag[peak_index],
    Peak_lag_hours = te_df$Lag_time[peak_index],
    Peak_TE_significant = significant[peak_index],
    Significant_lag_count = significant_count,
    Cumulative_TE = cumulative_te,
    Max_lag_steps = spec$max_lag_steps,
    Max_lag_hours = max(te_df$Lag_time),
    Notes = notes
  )
}

summary <- data.table::rbindlist(rows, fill = TRUE)
summary_file <- "04_Results/Tables/All_sites_TE_summary.csv"
data.table::fwrite(summary, summary_file)
# Keep the earlier filename as a compatibility copy with the expanded fields.
data.table::fwrite(summary, "04_Results/Tables/All_sites_peak_TE.csv")

fmt <- function(x, digits = 5L) format(x, digits = digits, trim = TRUE, scientific = FALSE)
md_escape <- function(x) gsub("\\|", "\\\\|", x)

summary_rows <- vapply(seq_len(nrow(summary)), function(i) {
  x <- summary[i]
  paste0(
    "| ", md_escape(x$Site), " | ", md_escape(x$Variable_pair), " | ",
    fmt(x$Native_resolution_minutes, 5), " min | ", x$Data_length_observations, " | ",
    substr(x$Start_date, 1, 10), " | ", substr(x$End_date, 1, 10), " | ",
    fmt(x$Peak_TE, 5), " | ", x$Peak_lag_steps, " (", fmt(x$Peak_lag_hours, 5), " h) | ",
    ifelse(x$Peak_TE_significant, "yes", "no"), " | ", fmt(x$Cumulative_TE, 5), " | ",
    md_escape(x$Notes), " |"
  )
}, character(1))

detail_sections <- unlist(lapply(seq_len(nrow(summary)), function(i) {
  x <- summary[i]
  c(
    paste0("### ", x$Analysis_ID),
    "",
    paste0("- Site and pair: ", x$Site, "; ", x$Variable_pair),
    paste0("- Record: ", x$Data_length_observations, " complete observations from ", x$Start_date, " through ", x$End_date, "."),
    paste0("- Native resolution: ", fmt(x$Native_resolution_minutes), " minutes; tested lags: 0-", x$Max_lag_steps, " steps (0-", fmt(x$Max_lag_hours), " hours)."),
    paste0("- Peak TE: ", fmt(x$Peak_TE), " at ", x$Peak_lag_steps, " steps (", fmt(x$Peak_lag_hours), " hours); above shuffled threshold: ", ifelse(x$Peak_TE_significant, "yes", "no"), "."),
    paste0("- Significant lag positions: ", x$Significant_lag_count, "; cumulative TE: ", fmt(x$Cumulative_TE), "."),
    paste0("- Notes: ", x$Notes),
    paste0("- [Full lag table](../04_Results/Tables/TE_df_", x$Analysis_ID, ".csv) | [PNG figure](../04_Results/Figures/TE_lag_", x$Analysis_ID, ".png) | [PDF figure](../04_Results/Figures/TE_lag_", x$Analysis_ID, ".pdf)"),
    ""
  )
}), use.names = FALSE)

report <- c(
  "# Sampler Platter transfer-entropy results",
  "",
  "## Methods",
  "",
  "Each analysis tests lags 0 through 72 in the dataset's native observation steps. `Lag_time` and every plot x-axis convert those native steps to hours using the modal timestamp interval. Thus the physical horizon varies with native resolution.",
  "",
  "The generalized workflow preserves the project histogram estimator, quantile-folded extremes, optional structural-zero bins, and 300 shuffled null realizations. Lag windows are retained only when timestamps are contiguous at the detected native resolution.",
  "",
  "Temperature analyses use five-day same-clock anomalies, consistent with the prior project workflow. Hubbard Brook growing-season runs retain April through October. Casper is omitted as requested. No stream-temperature-to-discharge direction is analyzed.",
  "",
  "Cumulative TE is defined here exactly as requested: the sum of TE values that exceed their shuffled TE threshold, divided by 72 configured lag steps. The lag-zero position can contribute to the numerator, but the fixed denominator remains 72.",
  "",
  "## Results summary",
  "",
  "| Site | Variable pair | Native resolution | N | Start | End | Peak TE | Peak lag (steps and hours) | Peak significant | Cumulative TE | Notes |",
  "|---|---|---:|---:|---|---|---:|---:|---:|---:|---|",
  summary_rows,
  "",
  "## Analysis details",
  "",
  detail_sections,
  "## Interpretation notes",
  "",
  "A peak at lag zero describes same-timestamp predictive information under this estimator and should not be interpreted as a delayed causal response. A peak marked as not significant did not exceed the shuffled threshold. Cumulative TE is a comparative summary across the common 72-step lag range; because physical horizons differ by native resolution, cross-resolution comparisons should retain that distinction.",
  "",
  "## Reproducibility and verification",
  "",
  "Configuration and dataset decisions are in `00_Data/Site_info.csv`. The complete input audit is in `04_Results/Data_Audit/Input_data_summary.csv`, and the verification manifest is in `04_Results/verification_manifest.csv`."
)
writeLines(report, "03_Reports/TE_analysis_report.md")
message("Wrote expanded summary and report for ", nrow(summary), " analyses.")

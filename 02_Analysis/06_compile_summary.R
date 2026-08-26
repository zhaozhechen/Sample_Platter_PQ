# Build the cross-site summary and a complete HTML report from TE outputs.

.libPaths(c(normalizePath("R_library"), .libPaths()))
source("01_Functions/TE_implementation.R")
args <- commandArgs(trailingOnly = TRUE)
report_only <- "--report-only" %in% args
output_arg <- grep("^--output=", args, value = TRUE)
report_file <- if (length(output_arg)) {
  sub("^--output=", "", output_arg[[1L]])
} else {
  "03_Reports/TE_analysis_report.html"
}

.html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

.html_table <- function(df, digits = 6L) {
  display <- as.data.frame(df, stringsAsFactors = FALSE)
  for (name in names(display)) {
    if (is.numeric(display[[name]])) {
      display[[name]] <- format(
        display[[name]], digits = digits, trim = TRUE, scientific = FALSE
      )
    } else if (is.logical(display[[name]])) {
      display[[name]] <- ifelse(display[[name]], "Yes", "No")
    }
  }
  header <- paste0("<th>", .html_escape(names(display)), "</th>", collapse = "")
  body <- vapply(seq_len(nrow(display)), function(i) {
    paste0(
      "<tr>",
      paste0("<td>", .html_escape(display[i, ]), "</td>", collapse = ""),
      "</tr>"
    )
  }, character(1))
  paste0(
    "<div class=\"table-wrap\"><table><thead><tr>", header,
    "</tr></thead><tbody>", paste0(body, collapse = ""),
    "</tbody></table></div>"
  )
}

.base64_encode <- function(path) {
  bytes <- as.integer(readBin(path, what = "raw", n = file.info(path)$size))
  alphabet <- strsplit(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/",
    "",
    fixed = TRUE
  )[[1]]
  padding <- (3L - length(bytes) %% 3L) %% 3L
  padded <- c(bytes, rep.int(0L, padding))
  triplets <- matrix(padded, ncol = 3L, byrow = TRUE)
  indices <- cbind(
    bitwShiftR(triplets[, 1L], 2L),
    bitwOr(bitwShiftL(bitwAnd(triplets[, 1L], 3L), 4L), bitwShiftR(triplets[, 2L], 4L)),
    bitwOr(bitwShiftL(bitwAnd(triplets[, 2L], 15L), 2L), bitwShiftR(triplets[, 3L], 6L)),
    bitwAnd(triplets[, 3L], 63L)
  )
  encoded <- alphabet[as.vector(t(indices)) + 1L]
  if (padding > 0L) encoded[(length(encoded) - padding + 1L):length(encoded)] <- "="
  paste0(encoded, collapse = "")
}

.image_data_uri <- function(path) {
  if (!file.exists(path)) stop("Missing report figure: ", path)
  paste0("data:image/png;base64,", .base64_encode(path))
}

.format_metric <- function(value, suffix = "", digits = 6L) {
  if (length(value) != 1L || !is.finite(value)) return("NA")
  paste0(format(value, digits = digits, trim = TRUE), suffix)
}

site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))
audit <- data.table::fread("04_Results/Data_Audit/Input_data_summary.csv")
configured <- site_info[run == TRUE]
rows <- vector("list", nrow(configured))

for (i in seq_len(nrow(configured))) {
  spec <- configured[i]
  result_file <- file.path("04_Results", "Tables", paste0("TE_df_", spec$Analysis_ID, ".csv"))
  if (!file.exists(result_file)) stop("Missing result table: ", result_file)
  te_df <- data.table::fread(result_file)
  data.table::setorder(te_df, Lag)
  audit_row <- audit[Analysis_ID == spec$Analysis_ID]
  if (nrow(audit_row) != 1L) stop("Missing or duplicate audit row: ", spec$Analysis_ID)

  significant <- is.finite(te_df$TEnorm) & is.finite(te_df$TEnormcrit) &
    te_df$TEnorm > te_df$TEnormcrit
  significant_count <- sum(significant)
  if (significant_count > 0L) {
    significant_indices <- which(significant)
    peak_index <- significant_indices[which.max(te_df$TEnorm[significant_indices])]
    peak_te <- te_df$TEnorm[peak_index]
    peak_lag <- te_df$Lag_time[peak_index]
    cumulative_te <- mean(te_df$TEnorm[significant_indices])
  } else {
    peak_te <- NA_real_
    peak_lag <- NA_real_
    cumulative_te <- NA_real_
  }
  memory_hours <- NA_real_
  if (length(significant) > 0L && isTRUE(significant[[1L]])) {
    first_insignificant <- which(!significant)
    if (length(first_insignificant) > 0L) {
      memory_hours <- te_df$Lag_time[first_insignificant[[1L]]]
    }
  }
  excluded_zero <- if ("excluded_Q_zero_observations" %in% names(audit_row)) {
    audit_row$excluded_Q_zero_observations
  } else {
    0L
  }
  notes <- spec$provenance_note
  if (excluded_zero > 0) {
    notes <- paste0(
      notes, "; excluded ", excluded_zero,
      " Q=0 observations from TE and retained their timestamps as gaps"
    )
  }
  if (!significant_count) notes <- paste(notes, "No TE value exceeded the shuffled threshold.")

  rows[[i]] <- data.table::data.table(
    Analysis_ID = spec$Analysis_ID,
    Site = spec$Site_Name,
    Variable_pair = paste0(spec$source_label, " -> ", spec$sink_label),
    Native_resolution_hours = audit_row$step_hours,
    Native_resolution_minutes = audit_row$step_hours * 60,
    Data_length_observations = audit_row$rows,
    Excluded_Q_zero_observations = excluded_zero,
    Data_duration_days = audit_row$duration_days,
    Start_date = audit_row$start,
    End_date = audit_row$end,
    Peak_normalized_TE_percent = peak_te,
    Peak_lag_hours = peak_lag,
    Significant_lag_count = significant_count,
    Cumulative_normalized_TE_percent = cumulative_te,
    Memory_hours = memory_hours,
    Max_lag_steps = spec$max_lag_steps,
    Max_lag_hours = max(te_df$Lag_time),
    Notes = notes
  )
}

summary <- data.table::rbindlist(rows, fill = TRUE)
if (!report_only) {
  data.table::fwrite(summary, "04_Results/Tables/All_sites_TE_summary.csv", na = "NA")
  data.table::fwrite(summary, "04_Results/Tables/All_sites_peak_TE.csv", na = "NA")
}

summary_display <- summary[, .(
  Site,
  Variable_pair,
  Resolution_minutes = Native_resolution_minutes,
  Observations = Data_length_observations,
  Q_zero_excluded = Excluded_Q_zero_observations,
  Start = substr(Start_date, 1, 10),
  End = substr(End_date, 1, 10),
  Peak_normalized_TE_percent,
  Peak_lag_hours,
  Cumulative_normalized_TE_percent,
  Memory_hours,
  Notes
)]

css <- paste0(
  "body{margin:0;background:#f4f7f8;color:#24303a;font-family:Segoe UI,Arial,sans-serif;line-height:1.6}",
  "main{width:min(1320px,calc(100% - 32px));margin:32px auto;padding:48px 56px;background:#fff;border:1px solid #d9e1e5;border-radius:12px;box-shadow:0 8px 28px rgba(36,48,58,.08)}",
  "h1,h2,h3{line-height:1.25}h1{margin:0 0 8px;font-size:2.25rem}h2{margin-top:42px;padding-bottom:8px;border-bottom:2px solid #e8f4f0}h3{margin-top:34px;color:#176f61}",
  ".subtitle{color:#61717d}.notice{margin:22px 0;padding:14px 18px;background:#e8f4f0;border-left:4px solid #2a8f78;border-radius:4px}",
  ".table-wrap{overflow-x:auto}table{width:100%;margin:16px 0 24px;border-collapse:collapse;font-variant-numeric:tabular-nums;font-size:.92rem}",
  "th,td{padding:8px 10px;border:1px solid #d9e1e5;text-align:right;vertical-align:top}th{background:#f0f5f6;white-space:nowrap}th:first-child,td:first-child{text-align:left}",
  "code{padding:2px 5px;background:#f0f3f4;border-radius:4px;font-family:Consolas,monospace}img{display:block;width:100%;height:auto;margin:22px auto;border:1px solid #d9e1e5}",
  "a{color:#176f61}.analysis{margin-top:38px;padding-top:6px}.meta{padding:14px 18px;background:#f8fafb;border:1px solid #d9e1e5;border-radius:6px}.outputs{columns:3}.outputs li{margin-bottom:7px;break-inside:avoid}",
  "@media(max-width:760px){main{padding:28px 22px}.outputs{columns:1}}@media print{body{background:white}main{width:100%;margin:0;padding:0;border:0;box-shadow:none}}"
)

metrics <- data.frame(
  Metric = c("TE", "TE threshold", "Normalized TE", "MI", "Correlation", "Cumulative TE", "Memory"),
  Meaning = c(
    "Directional information from lagged source to sink, in bits.",
    "Shuffled critical value used to flag TE above the null reference.",
    "TE divided by sink entropy and multiplied by 100.",
    "Non-directional mutual information between lagged source and current sink.",
    "Pearson correlation between lagged source and current sink.",
    "Mean of the significant normalized TE values (sum divided by the number of significant lag times); reported as percent uncertainty reduction.",
    "First lag, in hours, at which an initially significant TE sequence becomes insignificant. Memory is NA if TE is not significant at lag zero or remains significant through the evaluated lags."
  ),
  stringsAsFactors = FALSE
)

html <- c(
  "<!doctype html>",
  "<html lang=\"en\"><head><meta charset=\"utf-8\">",
  "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
  "<title>Sampler Platter Transfer-Entropy Results</title>",
  paste0("<style>", css, "</style></head><body><main>"),
  "<h1>Sampler Platter Transfer-Entropy Results</h1>",
  paste0("<p class=\"subtitle\">Generated ", format(Sys.Date(), "%B %d, %Y"), "</p>"),
  "<div class=\"notice\">This report summarizes all configured site and variable-pair analyses, their native temporal resolution, record coverage, and TE results.</div>",
  "<h2>Methods</h2>",
  "<p>Each analysis evaluates lags 0 through 72 in native observation steps. Tables and figure x-axes convert those steps to hours using the native timestamp interval, so the temporal coverage varies by dataset resolution.</p>",
  "<p>TE algorithm uses histogram discretization, quantile-folded extremes, optional structural-zero bins, and 300 shuffled null realizations. </p>",
  "<div class=\"notice\"><strong>Q to stream temperature:</strong> observations with missing Q or Q equal to zero are excluded from TE. Their timestamps are not collapsed, so lag windows cannot connect values across those gaps. Q=0 periods are shaded in both time-series panels.</div>",
  "<p>Temperature analyses use five-day same-clock anomalies. Hubbard Brook growing-season runs retain April through October. </p>",
  "<h2>Metrics</h2>",
  .html_table(metrics),
  "<h2>Results summary</h2>",
  .html_table(summary_display)
)

for (i in seq_len(nrow(summary))) {
  x <- summary[i]
  safe_id <- x$Analysis_ID
  figure_file <- file.path("04_Results", "Figures", paste0("TE_lag_", safe_id, ".png"))
  result_table <- data.table::fread(file.path(
    "04_Results", "Tables", paste0("TE_df_", safe_id, ".csv")
  ))
  result_display <- data.frame(
    Statistic = c(
      "Native resolution", "Complete TE observations", "Q=0 observations excluded",
      "Start", "End", "Maximum lag", "Peak normalized TE", "Peak lag",
      "Significant lag positions", "Cumulative normalized TE", "Memory"
    ),
    Value = c(
      paste0(format(x$Native_resolution_minutes, trim = TRUE), " minutes"),
      x$Data_length_observations,
      x$Excluded_Q_zero_observations,
      x$Start_date,
      x$End_date,
      paste0(x$Max_lag_steps, " steps (", format(x$Max_lag_hours, trim = TRUE), " hours)"),
      .format_metric(x$Peak_normalized_TE_percent, " %"),
      .format_metric(x$Peak_lag_hours, " hours"),
      x$Significant_lag_count,
      .format_metric(x$Cumulative_normalized_TE_percent, " %"),
      .format_metric(x$Memory_hours, " hours")
    ),
    stringsAsFactors = FALSE
  )
  html <- c(
    html,
    paste0("<section class=\"analysis\"><h3>", i, ". ", .html_escape(x$Site), ": ", .html_escape(x$Variable_pair), "</h3>"),
    paste0("<div class=\"meta\"><strong>Analysis ID:</strong> ", .html_escape(safe_id), "<br><strong>Data notes:</strong> ", .html_escape(x$Notes), "</div>"),
    .html_table(result_display),
    paste0("<img src=\"", .image_data_uri(figure_file), "\" alt=\"TE diagnostic figure for ", .html_escape(safe_id), "\">"),
    "</section>"
  )
}

html <- c(
  html,
  "<p></p>",
  "</main></body></html>"
)

writeLines(html, report_file, useBytes = TRUE)
message("Wrote expanded summary and HTML report for ", nrow(summary), " analyses to ", report_file, ".")

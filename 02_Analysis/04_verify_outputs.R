# Verify completeness and physical lag axes for every configured analysis.

.libPaths(c(normalizePath("R_library"), .libPaths()))
site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))[run == TRUE]
input_summary <- data.table::fread("04_Results/Data_Audit/Input_data_summary.csv")
status <- data.table::fread("02_Analysis/logs/run_status.csv")

if (nrow(status[status != "success"]) > 0L) stop("At least one analysis is not successful.")
if (!setequal(status$Analysis_ID, site_info$Analysis_ID)) stop("Run-status IDs do not match configuration.")

checks <- lapply(seq_len(nrow(site_info)), function(i) {
  spec <- site_info[i]
  audit <- input_summary[Analysis_ID == spec$Analysis_ID]
  table_file <- file.path("04_Results/Tables", paste0("TE_df_", spec$Analysis_ID, ".csv"))
  png_file <- file.path("04_Results/Figures", paste0("TE_lag_", spec$Analysis_ID, ".png"))
  pdf_file <- file.path("04_Results/Figures", paste0("TE_lag_", spec$Analysis_ID, ".pdf"))
  if (!file.exists(table_file)) stop("Missing table: ", table_file)
  if (!file.exists(png_file) || file.info(png_file)$size < 10000) stop("Missing/empty PNG: ", png_file)
  if (!file.exists(pdf_file) || file.info(pdf_file)$size < 10000) stop("Missing/empty PDF: ", pdf_file)
  result <- data.table::fread(table_file)
  required <- c(
    "Lag", "Lag_time", "Lag_unit", "N_complete_windows", "TE", "TEcrit",
    "TEnorm", "TEnormcrit", "MI", "MIcrit", "Corr", "Corrcrit"
  )
  if (!all(required %in% names(result))) stop("Missing result columns for ", spec$Analysis_ID)
  expected_rows <- audit$max_lag_steps + 1L
  if (nrow(result) != expected_rows) stop("Wrong lag-row count for ", spec$Analysis_ID)
  if (!all(result$Lag_unit == "hours")) stop("Non-hour lag unit for ", spec$Analysis_ID)
  if (abs(min(result$Lag_time) - 0) > 1e-8) stop("Lag axis does not start at zero for ", spec$Analysis_ID)
  expected_max_hours <- spec$max_lag_steps * audit$step_hours
  if (abs(max(result$Lag_time) - expected_max_hours) > 1e-4) stop("Lag axis does not match 72 native steps for ", spec$Analysis_ID)
  if (any(result$N_complete_windows <= 0)) stop("Empty lag window for ", spec$Analysis_ID)
  if (any(!is.finite(result$TE)) || any(!is.finite(result$TEcrit))) stop("Non-finite TE result for ", spec$Analysis_ID)
  data.frame(
    Analysis_ID = spec$Analysis_ID,
    rows = nrow(result),
    lag_step_hours = if (nrow(result) > 1L) result$Lag_time[2] - result$Lag_time[1] else NA_real_,
    max_lag_steps = spec$max_lag_steps,
    max_lag_hours = max(result$Lag_time),
    min_complete_windows = min(result$N_complete_windows),
    table_bytes = file.info(table_file)$size,
    png_bytes = file.info(png_file)$size,
    pdf_bytes = file.info(pdf_file)$size,
    table_md5 = unname(tools::md5sum(table_file)),
    png_md5 = unname(tools::md5sum(png_file)),
    pdf_md5 = unname(tools::md5sum(pdf_file)),
    stringsAsFactors = FALSE
  )
})

manifest <- data.table::rbindlist(checks)
data.table::fwrite(manifest, "04_Results/verification_manifest.csv")
report_file <- "03_Reports/TE_analysis_report.html"
if (!file.exists(report_file) || file.info(report_file)$size < 10000) {
  stop("Missing or empty HTML report: ", report_file)
}
if (file.exists("03_Reports/TE_analysis_report.md")) {
  stop("Obsolete Markdown report is still present.")
}
message("VERIFICATION PASSED for ", nrow(manifest), " analyses")

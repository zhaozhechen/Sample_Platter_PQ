# Run every unambiguous site/variable-pair analysis in 00_Data/Site_info.csv.

.libPaths(c(normalizePath("R_library"), .libPaths()))
source("01_Functions/Data_loading.R")
source("01_Functions/TE_implementation.R")

workers <- max(1L, min(4L, future::availableCores() - 1L))
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))
site_info <- site_info[run == TRUE]
dir.create("02_Analysis/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Results/Tables", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Results/Figures", recursive = TRUE, showWarnings = FALSE)

peak_rows <- list()
status_rows <- list()

for (i in seq_len(nrow(site_info))) {
  spec <- site_info[i]
  started <- Sys.time()
  message(sprintf("[%d/%d] %s", i, nrow(site_info), spec$Analysis_ID))
  result <- tryCatch({
    site_df <- load_site_analysis(spec)
    te_df <- run_TE(
      site_df = site_df,
      Site_ID = spec$Analysis_ID,
      Site_Name = spec$Site_Name,
      source_varname = spec$source_label,
      sink_varname = spec$sink_label,
      max_lag_steps = spec$max_lag_steps,
      nshuffle = spec$nshuffle,
      ZFlagSource = spec$zflag_source,
      ZFlagSink = spec$zflag_sink,
      seed = 111L + i,
      parallel = workers > 1L
    )
    best <- .TE_best_metric(te_df, "TE")
    peak_rows[[length(peak_rows) + 1L]] <- data.frame(
      Analysis_ID = spec$Analysis_ID,
      Site_Name = spec$Site_Name,
      source_label = spec$source_label,
      sink_label = spec$sink_label,
      step_hours = attr(site_df, "step_hours"),
      rows = nrow(site_df),
      detected_TE_lag_hours = best$lag,
      peak_TE = te_df$TE[best$index],
      TE_threshold = te_df$TEcrit[best$index],
      significant = best$significant,
      stringsAsFactors = FALSE
    )
    list(status = "success", error = "")
  }, error = function(e) list(status = "failed", error = conditionMessage(e)))

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    Analysis_ID = spec$Analysis_ID,
    status = result$status,
    error = result$error,
    elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(
    data.table::rbindlist(status_rows, fill = TRUE),
    "02_Analysis/logs/run_status.csv"
  )
  if (result$status == "failed") message("FAILED: ", result$error)
  gc()
}

if (length(peak_rows)) {
  data.table::fwrite(
    data.table::rbindlist(peak_rows, fill = TRUE),
    "04_Results/Tables/All_sites_peak_TE.csv"
  )
}
message("All configured analyses attempted. See 02_Analysis/logs/run_status.csv")

# Run or rerun one configured analysis.
# Usage: Rscript 02_Analysis/03_run_one_TE.R Baltimore_P_Q

.libPaths(c(normalizePath("R_library"), .libPaths()))
source("01_Functions/Data_loading.R")
source("01_Functions/TE_implementation.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply exactly one Analysis_ID.")
analysis_id <- args[1]

site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))
spec <- site_info[Analysis_ID == analysis_id]
if (nrow(spec) != 1L) stop("Analysis_ID not found or not unique: ", analysis_id)

workers <- max(1L, min(4L, future::availableCores() - 1L))
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)
started <- Sys.time()

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
  seed = 111L + which(site_info$Analysis_ID == analysis_id),
  parallel = workers > 1L
)

best <- .TE_best_metric(te_df, "TE")
peak <- data.table::data.table(
  Analysis_ID = spec$Analysis_ID,
  Site_Name = spec$Site_Name,
  source_label = spec$source_label,
  sink_label = spec$sink_label,
  step_hours = attr(site_df, "step_hours"),
  rows = nrow(site_df),
  detected_TE_lag_hours = best$lag,
  peak_TE = te_df$TE[best$index],
  TE_threshold = te_df$TEcrit[best$index],
  significant = best$significant
)
peak_file <- "04_Results/Tables/All_sites_peak_TE.csv"
peaks <- if (file.exists(peak_file)) data.table::fread(peak_file) else peak[0]
peaks <- peaks[Analysis_ID != analysis_id]
data.table::fwrite(data.table::rbindlist(list(peaks, peak), fill = TRUE), peak_file)

status_file <- "02_Analysis/logs/run_status.csv"
status <- if (file.exists(status_file)) data.table::fread(status_file) else data.table::data.table()
status <- status[Analysis_ID != analysis_id]
new_status <- data.table::data.table(
  Analysis_ID = analysis_id,
  status = "success",
  error = "",
  elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))
)
data.table::fwrite(data.table::rbindlist(list(status, new_status), fill = TRUE), status_file)
message("Completed ", analysis_id)

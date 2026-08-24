# Author: Zhaozhe Chen
# Generalized Sampler Platter data loading and preprocessing

.SP_parse_time <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "UTC"))
  x <- as.character(x)
  x <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", x)
  x <- sub("Z$", "+0000", x)
  formats <- c(
    "%Y-%m-%d %H:%M:%OS%z", "%Y-%m-%dT%H:%M:%OS%z",
    "%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS",
    "%m/%d/%Y %H:%M:%OS", "%m/%d/%Y %H:%M",
    "%m/%d/%Y", "%Y-%m-%d"
  )
  out_num <- rep(NA_real_, length(x))
  for (fmt in formats) {
    remaining <- which(is.na(out_num) & !is.na(x) & nzchar(x))
    if (!length(remaining)) break
    parsed <- suppressWarnings(strptime(x[remaining], format = fmt, tz = "UTC"))
    good <- !is.na(parsed)
    if (any(good)) out_num[remaining[good]] <- as.numeric(parsed[good])
  }
  # Some USGS downloads switch from formatted UTC timestamps to Unix seconds
  # within the same column. Parse those remaining values explicitly.
  remaining <- which(is.na(out_num) & !is.na(x) & nzchar(x))
  if (length(remaining)) {
    epoch <- suppressWarnings(as.numeric(x[remaining]))
    plausible <- is.finite(epoch) & epoch > 1e8 & epoch < 4e9
    out_num[remaining[plausible]] <- epoch[plausible]
  }
  as.POSIXct(out_num, origin = "1970-01-01", tz = "UTC")
}

.SP_mode_step_hours <- function(time) {
  delta <- diff(as.numeric(time)) / 3600
  delta <- delta[is.finite(delta) & delta > 0]
  if (!length(delta)) stop("No positive time step could be detected.")
  rounded <- round(delta, 8)
  counts <- sort(table(rounded), decreasing = TRUE)
  step <- as.numeric(names(counts)[1])
  list(
    step_hours = step,
    regular_fraction = mean(abs(delta - step) <= max(1e-7, step * 1e-6)),
    positive_differences = length(delta)
  )
}

.SP_apply_filter <- function(dt, column, value) {
  if (is.na(column) || !nzchar(column)) return(dt)
  if (!column %in% names(dt)) stop("Filter column not found: ", column)
  dt[as.character(get(column)) == as.character(value)]
}

.SP_clock_anomaly <- function(time, value, width_days = 5L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  work <- data.table::data.table(
    row_id = seq_along(value),
    time = time,
    value = value,
    clock = format(time, "%H:%M:%S", tz = "UTC")
  )
  data.table::setorder(work, clock, time)
  work[, moving_mean := data.table::frollmean(
    value, n = width_days, align = "center", fill = NA_real_, na.rm = TRUE
  ), by = clock]
  work[, anomaly := value - moving_mean]
  data.table::setorder(work, row_id)
  work$anomaly
}

.SP_transform <- function(time, value, method) {
  if (method == "raw") return(value)
  if (method == "anomaly_5day_clock") {
    return(.SP_clock_anomaly(time, value, width_days = 5L))
  }
  stop("Unknown transform: ", method)
}

load_site_analysis <- function(spec, data_root = "../../Data") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  path <- file.path(data_root, spec$data_file)
  if (!file.exists(path)) stop("Input file not found: ", path)
  has_sink_file <- !is.na(spec$sink_data_file) && nzchar(spec$sink_data_file)

  needed <- unique(c(
    spec$time_col, spec$source_col,
    if (!has_sink_file) spec$sink_col else NA_character_,
    spec$filter_col_1, spec$filter_col_2
  ))
  needed <- needed[!is.na(needed) & nzchar(needed)]
  dt <- data.table::fread(path, select = needed, showProgress = FALSE)
  dt <- .SP_apply_filter(dt, spec$filter_col_1, spec$filter_value_1)
  dt <- .SP_apply_filter(dt, spec$filter_col_2, spec$filter_value_2)

  source_data <- data.table::data.table(
    time = .SP_parse_time(dt[[spec$time_col]]),
    source_raw = suppressWarnings(as.numeric(dt[[spec$source_col]]))
  )

  input_files <- normalizePath(path, winslash = "/")
  if (has_sink_file) {
    sink_path <- file.path(data_root, spec$sink_data_file)
    if (!file.exists(sink_path)) stop("Sink input file not found: ", sink_path)
    sink_needed <- unique(c(spec$sink_time_col, spec$sink_col, spec$sink_filter_col))
    sink_needed <- sink_needed[!is.na(sink_needed) & nzchar(sink_needed)]
    sink_dt <- data.table::fread(sink_path, select = sink_needed, showProgress = FALSE)
    sink_dt <- .SP_apply_filter(sink_dt, spec$sink_filter_col, spec$sink_filter_value)
    sink_data <- data.table::data.table(
      time = .SP_parse_time(sink_dt[[spec$sink_time_col]]),
      sink_raw = suppressWarnings(as.numeric(sink_dt[[spec$sink_col]]))
    )
    source_data <- source_data[!is.na(time)]
    sink_data <- sink_data[!is.na(time)]
    data.table::setorder(source_data, time)
    data.table::setorder(sink_data, time)
    source_data <- source_data[!duplicated(time)]
    sink_data <- sink_data[!duplicated(time)]
    out <- merge(source_data, sink_data, by = "time", all = FALSE)
    input_files <- paste(
      input_files, normalizePath(sink_path, winslash = "/"), sep = " | "
    )
  } else {
    out <- source_data
    out[, sink_raw := suppressWarnings(as.numeric(dt[[spec$sink_col]]))]
  }
  out <- out[!is.na(time)]
  data.table::setorder(out, time)
  out <- out[!duplicated(time)]

  if (spec$season == "growing") {
    out <- out[as.integer(format(time, "%m", tz = "UTC")) %in% 4:10]
  } else if (spec$season != "all") {
    stop("Unknown season value: ", spec$season)
  }

  out[, source := .SP_transform(time, source_raw, spec$source_transform)]
  out[, sink := .SP_transform(time, sink_raw, spec$sink_transform)]
  plot_data <- out[, .(time, source, sink)]
  # W3 contains discharge outside the period with precipitation coverage.
  # Keep those rows out of the diagnostic time-series panels while retaining
  # all internal gaps within the overlapping P-Q period.
  if (grepl("^HB_w3_", as.character(spec$Analysis_ID))) {
    overlap <- out[is.finite(source) & is.finite(sink), time]
    if (length(overlap)) {
      plot_data <- plot_data[time >= min(overlap) & time <= max(overlap)]
    }
  }
  exclude_source_zero <- identical(as.character(spec$source_label), "Discharge") &&
    grepl("Stream temperature", as.character(spec$sink_label), fixed = TRUE)
  zero_excluded_rows <- if (exclude_source_zero) {
    sum(is.finite(out$source) & out$source == 0 & is.finite(out$sink))
  } else {
    0L
  }
  out <- out[is.finite(source) & is.finite(sink)]
  if (exclude_source_zero) out <- out[source != 0]
  data.table::setorder(out, time)

  if (nrow(out) < 10) stop("Fewer than 10 complete rows for ", spec$Analysis_ID)
  step <- .SP_mode_step_hours(out$time)
  attr(out, "input_file") <- input_files
  attr(out, "step_hours") <- step$step_hours
  attr(out, "regular_fraction") <- step$regular_fraction
  attr(out, "plot_data") <- plot_data
  attr(out, "exclude_source_zero") <- exclude_source_zero
  attr(out, "zero_excluded_rows") <- zero_excluded_rows
  out
}

summarize_site_input <- function(site_df, spec) {
  step_hours <- attr(site_df, "step_hours")
  data.frame(
    Analysis_ID = spec$Analysis_ID,
    Site_Name = spec$Site_Name,
    input_file = attr(site_df, "input_file"),
    source_column = spec$source_col,
    sink_column = spec$sink_col,
    source_transform = spec$source_transform,
    sink_transform = spec$sink_transform,
    season = spec$season,
    rows = nrow(site_df),
    duration_days = as.numeric(difftime(max(site_df$time), min(site_df$time), units = "days")),
    start = format(min(site_df$time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    end = format(max(site_df$time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    step_hours = step_hours,
    regular_transition_fraction = attr(site_df, "regular_fraction"),
    excluded_Q_zero_observations = attr(site_df, "zero_excluded_rows"),
    max_lag_steps = as.integer(spec$max_lag_steps),
    max_lag_hours = as.integer(spec$max_lag_steps) * step_hours,
    notes = spec$provenance_note,
    stringsAsFactors = FALSE
  )
}

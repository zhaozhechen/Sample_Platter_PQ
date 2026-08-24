# Author: Zhaozhe Chen
# Generalized transfer-entropy implementation for Sampler Platter v2
#
# Preserves the project estimator:
# TE(X -> Y) = H(Yt,Yt-1) + H(Yt-1,Xt-lag) - H(Yt-1)
#              - H(Yt,Yt-1,Xt-lag)
# Continuous values are discretized with quantile-folded extremes and optional
# structural-zero bins. A configured lag is a native observation step; plotted
# and tabulated lag values are converted to hours using the detected resolution.

.TE_lapply <- function(x, fun, parallel = FALSE) {
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required when parallel = TRUE.")
    }
    return(future.apply::future_lapply(x, fun, future.seed = TRUE))
  }
  lapply(x, fun)
}

.TE_entropy <- function(counts) {
  counts <- counts[counts > 0]
  if (!length(counts)) return(NA_real_)
  probabilities <- counts / sum(counts)
  -sum(probabilities * log2(probabilities))
}

.TE_metrics <- function(counts, nbins) {
  array_3d <- array(counts, dim = c(nbins, nbins, nbins))
  # joint_id is encoded with Y(t-1) varying fastest, then Y(t), then X.
  # These margins intentionally match the original Sampler Platter estimator.
  m_x <- apply(array_3d, 3, sum)
  m_y <- apply(array_3d, 2, sum)
  m_y1 <- apply(array_3d, 1, sum)
  m_xy <- apply(array_3d, c(2, 3), sum)
  m_yy1 <- apply(array_3d, c(1, 2), sum)
  m_y1x <- apply(array_3d, c(1, 3), sum)
  h_x <- .TE_entropy(m_x)
  h_y <- .TE_entropy(m_y)
  h_y1 <- .TE_entropy(m_y1)
  h_xy <- .TE_entropy(m_xy)
  h_yy1 <- .TE_entropy(m_yy1)
  h_y1x <- .TE_entropy(m_y1x)
  h_xyy1 <- .TE_entropy(counts)
  c(
    Hx = h_x,
    Hy = h_y,
    MI = h_x + h_y - h_xy,
    TE = h_yy1 + h_y1x - h_y1 - h_xyy1
  )
}

.TE_make_bins <- function(x, nbins, lower_qt, upper_qt, zero_adjust) {
  nonzero <- x[is.finite(x) & x != 0]
  if (length(unique(nonzero)) < 3) {
    stop("Each variable needs at least three distinct non-zero values.")
  }
  bounds <- as.numeric(stats::quantile(
    nonzero, probs = c(lower_qt, upper_qt), na.rm = TRUE, names = FALSE
  ))
  target <- if (zero_adjust) nonzero else x[is.finite(x)]
  interior <- target[target > bounds[1] & target < bounds[2]]
  if (length(unique(interior)) < 2) {
    stop("Insufficient interior variation after quantile folding.")
  }
  working_bins <- if (zero_adjust) nbins - 1L else nbins
  breaks <- seq(min(interior), max(interior), length.out = working_bins + 1L)
  breaks[1] <- breaks[1] - 1e-3
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-3

  out <- rep(NA_integer_, length(x))
  finite <- which(is.finite(x))
  if (zero_adjust) {
    zero <- finite[x[finite] == 0]
    nonzero_index <- finite[x[finite] != 0]
    out[zero] <- 1L
    cut_bins <- cut(
      x[nonzero_index], breaks = breaks, include.lowest = TRUE,
      right = FALSE, labels = FALSE
    )
    cut_bins[x[nonzero_index] <= bounds[1]] <- 1L
    cut_bins[x[nonzero_index] >= bounds[2]] <- working_bins
    out[nonzero_index] <- as.integer(cut_bins) + 1L
  } else {
    cut_bins <- cut(
      x[finite], breaks = breaks, include.lowest = TRUE,
      right = FALSE, labels = FALSE
    )
    cut_bins[x[finite] <= bounds[1]] <- 1L
    cut_bins[x[finite] >= bounds[2]] <- working_bins
    out[finite] <- as.integer(cut_bins)
  }
  out
}

.TE_window_indices <- function(time_numeric, lag, step_seconds) {
  n <- length(time_numeric)
  if (lag + 2L > n) return(list(x = integer(), yt = integer(), y1 = integer()))
  x_index <- 2L:(n - lag)
  yt_index <- (lag + 2L):n
  y1_index <- (lag + 1L):(n - 1L)
  step_tolerance <- max(1e-6, step_seconds * 1e-6)
  lag_tolerance <- max(step_tolerance, max(1L, lag) * step_seconds * 1e-6)
  regular <- abs(time_numeric[yt_index] - time_numeric[y1_index] - step_seconds) <= step_tolerance
  if (lag > 0L) {
    regular <- regular &
      abs(time_numeric[yt_index] - time_numeric[x_index] - lag * step_seconds) <= lag_tolerance
  } else {
    regular <- regular & abs(time_numeric[yt_index] - time_numeric[x_index]) <= step_tolerance
  }
  list(x = x_index[regular], yt = yt_index[regular], y1 = y1_index[regular])
}

.TE_counts <- function(x_bin, yt_bin, y1_bin, nbins) {
  good <- is.finite(x_bin) & is.finite(yt_bin) & is.finite(y1_bin)
  joint_id <- (x_bin[good] - 1L) * nbins^2L +
    (yt_bin[good] - 1L) * nbins + y1_bin[good]
  list(counts = tabulate(joint_id, nbins^3L), good = good)
}

.TE_shuffle_column <- function(x) {
  out <- x
  index <- which(!is.na(x))
  out[index] <- sample(x[index], length(index), replace = FALSE)
  out
}

.TE_critical <- function(
    source_bin, sink_bin, source_raw, sink_raw, indices, nbins,
    nshuffle, alpha, parallel) {
  x_bin <- source_bin[indices$x]
  yt_bin <- sink_bin[indices$yt]
  y1_bin <- sink_bin[indices$y1]
  x_raw <- source_raw[indices$x]
  yt_raw <- sink_raw[indices$yt]

  shuffled <- .TE_lapply(seq_len(nshuffle), function(i) {
    bx <- .TE_shuffle_column(x_bin)
    byt <- .TE_shuffle_column(yt_bin)
    by1 <- .TE_shuffle_column(y1_bin)
    counts <- .TE_counts(bx, byt, by1, nbins)$counts
    metrics <- .TE_metrics(counts, nbins)
    rx <- .TE_shuffle_column(x_raw)
    ryt <- .TE_shuffle_column(yt_raw)
    c(
      MI = unname(metrics["MI"]),
      TE = unname(metrics["TE"]),
      Corr = stats::cor(rx, ryt, use = "complete.obs")
    )
  }, parallel = parallel)
  null <- do.call(rbind, shuffled)
  threshold <- stats::qt(1 - alpha, df = 100)
  c(
    MIcrit = mean(null[, "MI"], na.rm = TRUE) + threshold * stats::sd(null[, "MI"], na.rm = TRUE),
    TEcrit = mean(null[, "TE"], na.rm = TRUE) + threshold * stats::sd(null[, "TE"], na.rm = TRUE),
    Corrcrit = mean(null[, "Corr"], na.rm = TRUE) + threshold * stats::sd(null[, "Corr"], na.rm = TRUE)
  )
}

Cal_TE_MI_main <- function(
    Source, Sink, Time, nbins = 11L, nshuffle = 300L, alpha = 0.05,
    Maxlag, ZFlagSink = TRUE, ZFlagSource = TRUE,
    Lag_Dependent_Crit = FALSE, lower_qt = 0.001,
    upper_qt = 1 - lower_qt, step_hours, parallel = FALSE) {
  if (Lag_Dependent_Crit) {
    stop("Lag-dependent thresholds are not supported in the optimized v2 engine.")
  }
  source_bin <- .TE_make_bins(Source, nbins, lower_qt, upper_qt, ZFlagSource)
  sink_bin <- .TE_make_bins(Sink, nbins, lower_qt, upper_qt, ZFlagSink)
  time_numeric <- as.numeric(Time)
  step_seconds <- step_hours * 3600
  max_indices <- .TE_window_indices(time_numeric, Maxlag, step_seconds)
  critical <- .TE_critical(
    source_bin, sink_bin, Source, Sink, max_indices, nbins,
    nshuffle, alpha, parallel
  )

  results <- lapply(0:Maxlag, function(lag) {
    index <- .TE_window_indices(time_numeric, lag, step_seconds)
    counted <- .TE_counts(
      source_bin[index$x], sink_bin[index$yt], sink_bin[index$y1], nbins
    )
    metrics <- .TE_metrics(counted$counts, nbins)
    good <- counted$good
    correlation <- stats::cor(
      Source[index$x][good], Sink[index$yt][good], use = "complete.obs"
    )
    data.frame(
      Lag = lag,
      Lag_time = lag * step_hours,
      Lag_unit = "hours",
      N_complete_windows = sum(good),
      MI = unname(metrics["MI"]),
      MIcrit = unname(critical["MIcrit"]),
      TE = unname(metrics["TE"]),
      TEcrit = unname(critical["TEcrit"]),
      Corr = correlation,
      Corrcrit = unname(critical["Corrcrit"]),
      Hx = unname(metrics["Hx"]),
      Hy = unname(metrics["Hy"])
    )
  })
  do.call(rbind, results)
}

.TE_theme <- function(base_size = 18) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.title = ggplot2::element_text(size = base_size, face = "bold"),
      axis.text = ggplot2::element_text(size = base_size * 0.82, color = "black"),
      plot.title = ggplot2::element_text(size = base_size * 1.08, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = base_size * 0.86),
      legend.title = ggplot2::element_text(size = base_size * 0.9, face = "bold")
    )
}

.TE_best_metric <- function(te_df, metric) {
  threshold <- paste0(metric, "crit")
  significant <- which(
    is.finite(te_df[[metric]]) & is.finite(te_df[[threshold]]) &
      te_df[[metric]] > te_df[[threshold]]
  )
  if (!length(significant)) {
    return(list(index = which.max(te_df[[metric]]), lag = NA_real_, significant = FALSE))
  }
  index <- significant[which.max(te_df[[metric]][significant])]
  list(index = index, lag = te_df$Lag_time[index], significant = TRUE)
}

TE_lag_plot <- function(te_df, metric, color = "#66C2A5", show_title = TRUE) {
  y_label <- c(
    TE = "TE (bits)", TEnorm = "Uncertainty reduction (%)",
    MI = "MI (bits)", Corr = "Correlation"
  )[[metric]]
  threshold <- paste0(metric, "crit")
  best <- .TE_best_metric(te_df, metric)
  plot_df <- data.frame(
    Lag_time = te_df$Lag_time,
    Metric = te_df[[metric]],
    Critical = te_df[[threshold]]
  )
  g <- ggplot2::ggplot(plot_df, ggplot2::aes(Lag_time, Metric)) +
    ggplot2::geom_line(color = color, linewidth = 1.05) +
    ggplot2::geom_line(
      ggplot2::aes(y = Critical), color = "#8DA0CB",
      linewidth = 0.85, linetype = "dashed"
    ) +
    ggplot2::labs(x = "Lag (hours)", y = y_label) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.2))) +
    .TE_theme()
  if (best$significant) {
    label <- paste0("Detected lag = ", format(best$lag, digits = 4, trim = TRUE), " h")
    x_range <- diff(range(plot_df$Lag_time, na.rm = TRUE))
    y_range <- diff(range(plot_df$Metric, na.rm = TRUE))
    if (!is.finite(y_range) || y_range == 0) y_range <- max(abs(plot_df$Metric), na.rm = TRUE)
    x_offset <- max(x_range * 0.055, .Machine$double.eps)
    raw_label_x <- if (best$lag > max(plot_df$Lag_time) * 0.62) {
      best$lag - x_offset
    } else {
      best$lag + x_offset
    }
    x_limits <- range(plot_df$Lag_time, na.rm = TRUE)
    label_x <- min(
      max(raw_label_x, x_limits[1] + x_range * 0.25),
      x_limits[2] - x_range * 0.25
    )
    label_hjust <- 0.5
    label_y <- plot_df$Metric[best$index] + max(y_range * 0.12, .Machine$double.eps)
    g <- g +
      ggplot2::geom_vline(
        xintercept = best$lag, color = "#FC8D62", linewidth = 0.9,
        linetype = "dotted"
      ) +
      ggplot2::geom_point(
        data = plot_df[best$index, ], color = "#FC8D62",
        fill = "white", shape = 21, size = 3.5, stroke = 1.2
      ) +
      ggplot2::annotate(
        "label", x = label_x, y = label_y, label = label,
        hjust = label_hjust, vjust = 0, size = 5.4,
        label.padding = grid::unit(0.22, "lines"), linewidth = 0.3,
        fill = "white", color = "black"
      )
  } else {
    g <- g + ggplot2::annotate(
      "label", x = min(plot_df$Lag_time), y = max(plot_df$Metric, na.rm = TRUE),
      label = "Detected lag = NA", hjust = 0, vjust = 1,
      size = 5.4, linewidth = 0.3, fill = "white"
    )
  }
  if (show_title) g <- g + ggplot2::ggtitle(y_label)
  g
}

.TE_zero_periods <- function(plot_data, step_hours) {
  zero <- plot_data[is.finite(plot_data$source) & plot_data$source == 0, ]
  if (!nrow(zero)) return(data.frame())
  step_seconds <- step_hours * 3600
  tolerance <- max(1e-6, step_seconds * 1e-6)
  breaks <- c(TRUE, abs(diff(as.numeric(zero$time)) - step_seconds) > tolerance)
  zero$period <- cumsum(breaks)
  periods <- lapply(split(zero$time, zero$period), function(x) {
    data.frame(
      xmin = min(x) - step_seconds / 2,
      xmax = max(x) + step_seconds / 2
    )
  })
  out <- do.call(rbind, periods)
  out$xmin <- as.POSIXct(out$xmin, origin = "1970-01-01", tz = "UTC")
  out$xmax <- as.POSIXct(out$xmax, origin = "1970-01-01", tz = "UTC")
  out
}

.TE_plot_series <- function(
    time, value, label, color, shade_periods = NULL,
    show_shade_legend = FALSE) {
  n <- length(value)
  delta <- diff(as.numeric(time))
  positive <- delta[is.finite(delta) & delta > 0]
  rounded <- round(positive, 6)
  expected <- as.numeric(names(sort(table(rounded), decreasing = TRUE)[1]))
  tolerance <- max(1e-6, expected * 1e-6)
  segment <- cumsum(c(TRUE, abs(delta - expected) > tolerance))
  index <- if (n > 60000L) unique(round(seq(1, n, length.out = 60000L))) else seq_len(n)
  plot_df <- data.frame(time = time[index], value = value[index], segment = segment[index])
  g <- ggplot2::ggplot(plot_df, ggplot2::aes(time, value, group = segment))
  if (!is.null(shade_periods) && nrow(shade_periods)) {
    g <- g + ggplot2::geom_rect(
      data = shade_periods,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                   fill = "Q = 0; excluded from TE"),
      inherit.aes = FALSE, alpha = 0.28, color = NA,
      show.legend = show_shade_legend
    ) + ggplot2::scale_fill_manual(
      values = c("Q = 0; excluded from TE" = "#F2D6A2"), name = NULL
    )
  }
  g +
    ggplot2::geom_line(color = color, linewidth = 0.45, na.rm = TRUE) +
    ggplot2::labs(x = NULL, y = label, title = label) +
    .TE_theme(16)
}

.TE_plot_hist <- function(value, label, color, remove_zero = FALSE, nbins = 11L) {
  value <- value[is.finite(value)]
  title <- "All values"
  if (remove_zero) {
    value <- value[value != 0]
    title <- "Zeros removed"
  }
  ggplot2::ggplot(data.frame(value = value), ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = if (remove_zero) nbins - 1L else nbins,
                            fill = color, color = "black", linewidth = 0.25) +
    ggplot2::labs(x = label, y = "Count", title = title) +
    .TE_theme(16)
}

plot_TE_result <- function(
    te_df, site_df, site_name, source_label, sink_label, output_file,
    nbins = 11L, dpi = 300) {
  colors <- c(source = "#8DA0CB", sink = "#FC8D62", metric = "#66C2A5")
  plot_data <- attr(site_df, "plot_data")
  if (is.null(plot_data)) plot_data <- site_df[, c("time", "source", "sink")]
  shade_periods <- NULL
  if (isTRUE(attr(site_df, "exclude_source_zero"))) {
    shade_periods <- .TE_zero_periods(plot_data, attr(site_df, "step_hours"))
  }
  source_row <- cowplot::plot_grid(
    .TE_plot_series(
      plot_data$time, plot_data$source, source_label, colors["source"],
      shade_periods, show_shade_legend = TRUE
    ),
    .TE_plot_hist(site_df$source, source_label, colors["source"], FALSE, nbins),
    .TE_plot_hist(site_df$source, source_label, colors["source"], TRUE, nbins),
    nrow = 1, rel_widths = c(2, 1, 1), align = "hv", axis = "tblr"
  )
  sink_row <- cowplot::plot_grid(
    .TE_plot_series(
      plot_data$time, plot_data$sink, sink_label, colors["sink"],
      shade_periods, show_shade_legend = FALSE
    ),
    .TE_plot_hist(site_df$sink, sink_label, colors["sink"], FALSE, nbins),
    .TE_plot_hist(site_df$sink, sink_label, colors["sink"], TRUE, nbins),
    nrow = 1, rel_widths = c(2, 1, 1), align = "hv", axis = "tblr"
  )
  metric_row <- cowplot::plot_grid(
    TE_lag_plot(te_df, "TE", colors["metric"], FALSE),
    TE_lag_plot(te_df, "TEnorm", colors["metric"], FALSE),
    TE_lag_plot(te_df, "MI", colors["metric"], FALSE),
    TE_lag_plot(te_df, "Corr", colors["metric"], FALSE),
    nrow = 1, align = "hv", axis = "tblr"
  )
  title <- cowplot::ggdraw() + cowplot::draw_label(
    paste0(site_name, ": ", source_label, " -> ", sink_label),
    x = 0.01, hjust = 0, fontface = "bold", size = 25
  )
  note <- cowplot::ggdraw() + cowplot::draw_label(
    paste0(
      "Solid = metric; dashed = shuffled threshold; dotted/point = strongest significant lag",
      if (isTRUE(attr(site_df, "exclude_source_zero")))
        "; light shading = Q = 0 observations excluded from TE" else ""
    ),
    x = 0.5, hjust = 0.5, size = 14
  )
  combined <- cowplot::plot_grid(
    title, source_row, sink_row, metric_row, note,
    ncol = 1, rel_heights = c(0.13, 0.85, 0.85, 0.9, 0.08)
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(output_file)
  ggplot2::ggsave(paste0(base, ".png"), combined, width = 18, height = 12,
                  units = "in", dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(base, ".pdf"), combined, width = 18, height = 12,
                  units = "in", bg = "white")
  invisible(c(png = paste0(base, ".png"), pdf = paste0(base, ".pdf")))
}

run_TE <- function(
    site_df, Site_ID, Site_Name, source_varname, sink_varname,
    n_bin = 11L, max_lag_steps = 72L, nshuffle = 300L, alpha = 0.05,
    ZFlagSource = TRUE, ZFlagSink = TRUE, lower_qt = 0.001,
    upper_qt = 1 - lower_qt, seed = 111L, parallel = FALSE,
    Output_path = "04_Results") {
  required <- c("time", "source", "sink")
  if (!all(required %in% names(site_df))) stop("site_df must contain time, source, and sink.")
  step_hours <- attr(site_df, "step_hours")
  if (!is.finite(step_hours) || step_hours <= 0) stop("site_df is missing a valid step_hours attribute.")
  max_lag <- as.integer(max_lag_steps)
  if (!is.finite(max_lag) || max_lag < 1L) stop("max_lag_steps must be a positive integer.")
  set.seed(seed)
  te_df <- Cal_TE_MI_main(
    Source = site_df$source, Sink = site_df$sink, Time = site_df$time,
    nbins = n_bin, nshuffle = nshuffle, alpha = alpha, Maxlag = max_lag,
    ZFlagSink = ZFlagSink, ZFlagSource = ZFlagSource,
    Lag_Dependent_Crit = FALSE, lower_qt = lower_qt, upper_qt = upper_qt,
    step_hours = step_hours, parallel = parallel
  )
  te_df$TEnorm <- te_df$TE / te_df$Hy * 100
  te_df$TEnormcrit <- te_df$TEcrit / te_df$Hy * 100
  te_df$TE_significant <- te_df$TE > te_df$TEcrit
  te_df$MI_significant <- te_df$MI > te_df$MIcrit

  safe_id <- gsub("[^A-Za-z0-9._-]+", "_", Site_ID)
  table_dir <- file.path(Output_path, "Tables")
  figure_dir <- file.path(Output_path, "Figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(te_df, file.path(table_dir, paste0("TE_df_", safe_id, ".csv")), row.names = FALSE)
  plot_TE_result(
    te_df, site_df, Site_Name, source_varname, sink_varname,
    file.path(figure_dir, paste0("TE_lag_", safe_id, ".png")), nbins = n_bin
  )
  te_df
}

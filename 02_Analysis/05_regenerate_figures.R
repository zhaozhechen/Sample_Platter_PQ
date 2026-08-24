# Regenerate figures from existing verified TE tables without recalculating TE.

.libPaths(c(normalizePath("R_library"), .libPaths()))
source("01_Functions/Data_loading.R")
source("01_Functions/TE_implementation.R")

site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))[run == TRUE]
for (i in seq_len(nrow(site_info))) {
  spec <- site_info[i]
  message(sprintf("[%d/%d] %s", i, nrow(site_info), spec$Analysis_ID))
  site_df <- load_site_analysis(spec)
  te_df <- data.table::fread(file.path(
    "04_Results/Tables", paste0("TE_df_", spec$Analysis_ID, ".csv")
  ))
  plot_TE_result(
    te_df = te_df,
    site_df = site_df,
    site_name = spec$Site_Name,
    source_label = spec$source_label,
    sink_label = spec$sink_label,
    output_file = file.path(
      "04_Results/Figures", paste0("TE_lag_", spec$Analysis_ID, ".png")
    )
  )
  gc()
}
message("Regenerated ", nrow(site_info), " figure pairs")


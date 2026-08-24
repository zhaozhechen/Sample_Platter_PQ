# Audit and verify every configured TE input without calculating TE.

.libPaths(c(normalizePath("R_library"), .libPaths()))
source("01_Functions/Data_loading.R")

site_info <- data.table::fread("00_Data/Site_info.csv", na.strings = c("", "NA"))
site_info <- site_info[run == TRUE]

summaries <- lapply(seq_len(nrow(site_info)), function(i) {
  spec <- site_info[i]
  message("Loading ", spec$Analysis_ID)
  site_df <- load_site_analysis(spec)
  summarize_site_input(site_df, spec)
})

summary_df <- data.table::rbindlist(summaries, fill = TRUE)
dir.create("04_Results/Data_Audit", recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(summary_df, "04_Results/Data_Audit/Input_data_summary.csv")
message("Wrote 04_Results/Data_Audit/Input_data_summary.csv")


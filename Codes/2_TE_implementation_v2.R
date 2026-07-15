# Author: Zhaozhe Chen
# Update Date: 2026.7.15

# This codes conducts TE from P -> Q for selected hourly sites

# ------- Global -------
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
# Source functions
source("Functions.R")

# Decide which site to process
arrayid <- 7

# Parameters for TE implementation ===============
n_bin <- 11 # Number of bins for TE discritization of continuous data (e.g., SM)
# Consider 3 days
max_lag <- 72 # Maximum lag to consider (This should be adjusted according to the processes and the temporal resolution of data)
Lag_Dependent_Crit <- FALSE # Determine if critical TE is lag-dependent
nshuffle <- 300 # Number of shuffles (bootstrap) for critical TE for statistical inference
alpha <- 0.05 # Confidence level for critical TE
# Set parallel session
plan(multisession,workers = availableCores()-1)
# Ensure reproducibility
set.seed(111)
# Determines if zero should be adjusted for the Sink and Source variables
ZFlagSink <- TRUE
ZFlagSource <- TRUE
# These are folding parameters to deal with extreme values (outliers) in the time series
# i.e., extreme values will be binned into the first or last bin
lower_qt <- 0.001
upper_qt <- 1-lower_qt

my_color <- brewer.pal(3,"Set2")

# Data path
Input_path <- "../../../Data/"

# Site info
site_info <- read.csv("../Site_info.csv")

# Output path
Output_path <- "../Results/TE_results_sites/"

# Get site info data ======
Site_ID <- site_info$Site_ID[arrayid]
# Data file name
filename <- site_info$filename[arrayid]
# Time variable name
time_name <- site_info$time_name[arrayid]
source_name <- site_info$source_name[arrayid]
sink_name <- site_info$sink_name[arrayid]
site_NO <- site_info$site_NO[arrayid]
site_filter <- site_info$filter[arrayid]

# Data processing =========
# Read in dataset
Site_df <- read.csv(paste0(Input_path,filename))

# Filter for sites are too long
if(!is.na(site_filter)){
  Site_df <- Site_df %>%
    filter(site_no == site_NO) %>%
    # This is for testing only. Cut some time for faster running
    slice_head(n = 100000)
}

# Filter site_ID for Atlanta sites
if(!is.na(site_NO)){
  # Select the target site
  Site_df <- Site_df %>%
    filter(site_no == site_NO)
}

Site_df <- Site_df %>%
  # Only select required variables
  select(Time = all_of(time_name),
         var1 = all_of(source_name),
         var2 = all_of(sink_name)) %>%
  # Format time 
  mutate(
    Time = parse_date_time(
      Time,
      orders = c(
        "Y-m-d H:M:S", "Y-m-d",
        "m/d/Y H:M:S", "m/d/Y H:M",
        "m/d/Y"
      )
    )
  ) %>%
  na.omit()

# Detrend hourly TS using 5-day windows =========
# Get 5-day moving window average
# Detect temporal resolution first
time_step <- diff(Site_df$Time)
time_step_hr <- as.numeric(
  names(sort(table(as.numeric(time_step, units = "hours")), decreasing = TRUE)[1])
)

# Number of observations in a 5-day window
window_size <- round(5 * 24 / time_step_hr)

Site_df <- Site_df %>%
  arrange(Time) %>%
  mutate(
    # Identify the same clock time across days
    clock_time = format(Time, "%H:%M:%S")
  ) %>%
  group_by(clock_time) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(
    # Five-day mean at the same clock time
    var1_same_time_mean = rollapply(
      var1,
      width = 5,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE,
      na.rm = TRUE
    ),
    var2_same_time_mean = rollapply(
      var2,
      width = 5,
      FUN = mean,
      fill = NA,
      align = "center",
      partial = TRUE,
      na.rm = TRUE
    ),
    
    # Anomalies
    var1_anom = var1 - var1_same_time_mean,
    var2_anom = var2 - var2_same_time_mean
  ) %>%
  ungroup() %>%
  select(-clock_time)

gc()

# Run TE using anomaly ======
run_TE(Site_ID,Site_df,source_name = "var1_anom",sink_name = "var2_anom")




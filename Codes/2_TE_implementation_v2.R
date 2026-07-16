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

# 6 G maximum allowed size
options(future.globals.maxSize = 4 * 1024^3) 

# Decide which site to process
for(arrayid in 16){
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
  Output_path <- "../Results/TE_results_sites_v3/"
  
  # Get site info data ======
  # Name for this TE implementation
  Site_ID <- site_info$Site_ID[arrayid]
  # Name for the site
  Site_Name <- site_info$Site_Name[arrayid]
  # Data file name
  filename <- site_info$filename[arrayid]
  # Time variable name
  time_name <- site_info$time_name[arrayid]
  # variable names
  source_name <- site_info$source_name[arrayid]
  sink_name <- site_info$sink_name[arrayid]
  # site no, specific for Atlanta sites
  site_NO <- site_info$site_NO[arrayid]
  # Filter data for growing season
  site_filter <- site_info$filter[arrayid]
  # Whether should use anomaly for source variable
  source_anomaly <- site_info$source_anomaly[arrayid]
  sink_anomaly <- site_info$sink_anomaly[arrayid]
  # Name to for the labels
  source_varname <- site_info$source_varname[arrayid]
  sink_varname <- site_info$sink_varname[arrayid]
  
  # Data processing =========
  # Read in dataset
  Site_df <- read.csv(paste0(Input_path,filename))
  
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
  
  # Apply optional seasonal filtering
  if (!is.na(site_filter)) {
    
    if (site_filter == "GS") {
      # Keep April through October
      Site_df <- Site_df %>%
        filter(month(Time) %in% 4:10)
      
    } else {
      warning(
        "Unknown site_filter value: ",
        site_filter,
        ". No filtering was applied."
      )
    }
  }
  
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
    select(-clock_time) %>%
    arrange(Time)
  
  gc()
  
  # Run TE ======
  if(!is.na(source_anomaly) & !is.na(sink_anomaly)){
    run_TE(Site_ID,Site_df,source_name = "var1_anom",sink_name = "var2_anom",
           source_varname = source_varname,sink_varname = sink_varname)
  }else if(!is.na(source_anomaly) & is.na(sink_anomaly)){
    run_TE(Site_ID,Site_df,source_name = "var1_anom",sink_name = "var2",
           source_varname = source_varname,sink_varname = sink_varname)
  }else if(is.na(source_anomaly) & is.na(sink_anomaly)){
    run_TE(Site_ID,Site_df,source_name = "var1",sink_name = "var2",
           source_varname = source_varname,sink_varname = sink_varname)
  }else if(is.na(source_anomaly) & !is.na(sink_anomaly)){
    run_TE(Site_ID,Site_df,source_name = "var1",sink_name = "var2_anom",
           source_varname = source_varname,sink_varname = sink_varname)
  }  
  
  message(arrayid)
}







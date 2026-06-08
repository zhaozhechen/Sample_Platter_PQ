# Author: Zhaozhe Chen
# Update Date: 2026.6.8

# This codes conducts TE from P -> Q for selected hourly sites

# ------- Global -------
library(dplyr)
library(tidyr)
library(lubridate)

# Data path ==========
Input_path <- "../../../Data/"
# Output path
Output_path <- "../Results/"
# This Site_ID is just for output
Site_ID <- "HB_w3"
# Data file name
filename <- "Hubbard_Brook_Alix/precip_discharge_w3.csv"

# Variable names in the dataset ====
# Time variable name
time_name <- "DATETIME"
# Precipitation variable name
P_name <- "precip"
# Q variable name
Q_name <- "avg_discharge"

# Source functions
source("Functions.R")

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

# ------- Main ---------
# Data processing ==========
# Read in dataset
Site_df <- read.csv(paste0(Input_path,filename)) %>%
  # Only select required variables
  select(Time = all_of(time_name),
         P = all_of(P_name),
         Q = all_of(Q_name)) %>%
  # Format time 
  mutate(
    Time = parse_date_time(
      Time,
      orders = c("Y-m-d H:M:S", "Y-m-d")
    )
  ) %>%
  na.omit()

# Implement hourly TE calculation ============
# Timing the TE calculation
start_time <- Sys.time()
# Run TE
TE_df <- Cal_TE_MI_main(Source = Site_df$P,
                        Sink = Site_df$Q,
                        nbins = n_bin,
                        Maxlag = max_lag,
                        alpha = alpha,
                        nshuffle = nshuffle,
                        upper_qt = upper_qt,
                        lower_qt = lower_qt,
                        ZFlagSource = ZFlagSource,
                        ZFlagSink = ZFlagSink,
                        Lag_Dependent_Crit = Lag_Dependent_Crit)
end_time <- Sys.time()
run_time <- end_time - start_time
message(run_time)

# Make TE vs lag plot
g_TE <- lag_plots_all(TE_df,Site_ID)

# Output TE df
write.csv(TE_df,paste0(Output_path,"/TE_results_sites/TE_df_",Site_ID,".csv"))
# Output TE plot
print_g(g_TE,
        paste0("/TE_results_sites/TE_lag_",Site_ID),
        14,13)



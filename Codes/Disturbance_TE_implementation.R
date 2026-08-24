# Author: Zhaozhe Chen
#Date: 2026.7.16

# This code it to run TE for the synthesis data for Disturbance project

# ------- Global -------
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
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

# Data path
Input_path <- "../../../Data/Synthesis_data/"

# Output path
Output_path <- "../Results/Disturbance_TE/"

# Desync =========
# Read in df
Site_df <- read.csv(paste0(Input_path,"synthetic_desync.csv")) %>%
  mutate(Time = time)
# Run pre-disturb
Site_df_pre <- Site_df %>%
  filter(phase == "pre")
Site_Name <- "Desync_pre"
run_TE(Site_ID = "Desync_pre",Site_df_pre,source_name = "climate",sink_name = "response",
       source_varname = "climate",sink_varname = "response")

# Run post-disturb
Site_df_pre <- Site_df %>%
  filter(phase == "post")
Site_Name <- "Desync_post"
run_TE(Site_ID = "Desync_post",Site_df_pre,source_name = "climate",sink_name = "response",
       source_varname = "climate",sink_varname = "response")

# Return =========
# Read in df
Site_df <- read.csv(paste0(Input_path,"synthetic_return.csv")) %>%
  mutate(Time = time)
# Run pre-disturb
Site_df_pre <- Site_df %>%
  filter(phase == "pre")
Site_Name <- "Return_pre"
run_TE(Site_ID = "Return_pre",Site_df_pre,source_name = "climate",sink_name = "response",
       source_varname = "climate",sink_varname = "response")

# Run post-disturb
Site_df_pre <- Site_df %>%
  filter(phase == "post")
Site_Name <- "Return_post"
run_TE(Site_ID = "Return_post",Site_df_pre,source_name = "climate",sink_name = "response",
       source_varname = "climate",sink_varname = "response")

# Shift =========
# Read in df
Site_df <- read.csv(paste0(Input_path,"synthetic_shift.csv")) %>%
  mutate(Time = time)
# Run pre-disturb
Site_df_pre <- Site_df %>%
  filter(phase == "pre")
Site_Name <- "Shift_pre"
run_TE(Site_ID = "Shift_pre",Site_df_pre,source_name = "climate",sink_name = "response",
       source_varname = "climate",sink_varname = "response")

# Run post-disturb
Site_df_pre <- Site_df %>%
  filter(phase == "post")
Site_Name <- "Shift_post"
run_TE(Site_ID = "Shift_post",Site_df_pre,source_name = "climate",sink_name = "response",
       source_varname = "climate",sink_varname = "response")






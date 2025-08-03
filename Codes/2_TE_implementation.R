# Author: Zhaozhe Chen
# Update Date: 2025.8.3

# This codes conducts TE from P -> Q for selected Macroshed sites

# ------- Global -------
library(dplyr)
library(here)
library(tidyr)

# Data path
Input_path <- here("../../../Data/")
# Output path
Output_path <- here("../Results/")
# Source functions
source(here("Functions.R"))

# Parameters for TE implementation
n_bin <- 11 # Number of bins for TE discritization of continuous data (e.g., SM)
max_lag <- 365 # Maximum lag to consider (This should be adjusted according to the processes and the temporal resolution of data)
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
# Data processing ---------------------
# Read in Macroshed dataset
MS_data <- read.csv(here(Input_path,"p_q_et_ms_data.csv"))
# Site info for tested sites
Site_info <- read.csv(here(Input_path,"refSites_for_PQanalysis.csv"))

# Get site names to be analyzed
Site_ID_ls <- Site_info$site_code
# Remove unwanted sites
Site_ID_ls <- Site_ID_ls[!Site_ID_ls %in% c(
  "WS-4", # No P
  "allequash_creek" # No P
)]

# Only keep data for these sites, and only keep required variables (P and Q)
MS_df <- MS_data %>%
  filter(site_code %in% Site_ID_ls,
         var == "discharge"|var == "precipitation")
# Remove problematic data
MS_df$val[!is.na(MS_df$ms_status) & MS_df$ms_status == 1] <- NA

# Initialize a list to store all TE results
TE_ls <- list()
# Initialize a list to store all TE figures
g_ls <- list()

# Loop over the sites
for(i in 1:length(Site_ID_ls)){
  Site_ID <- Site_ID_ls[i]
  # Get P and Q TS for this site
  P_df <- MS_df %>%
    filter(site_code == Site_ID,
           var == "precipitation") %>%
    select(date,P = val)
  Q_df <- MS_df %>%
    filter(site_code == Site_ID,
           var == "discharge") %>%
    select(date,Q = val)
  # Combine the two df to match their time
  Site_df <- merge(P_df,Q_df,by="date")
  Site_df$date <- as.Date(Site_df$date)
  # Sort the dataset
  Site_df <- Site_df[order(Site_df$date),]
  # Output this site df
  write.csv(Site_df,paste0(Output_path,"/Processed_data/PQ_df_",Site_ID,".csv"))
  # Implement TE -------------------------
  # Timing the TE calculation
  start_time <- Sys.time()
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
  run_time <- as.character(end_time - start_time)
  message(run_time)
  TE_ls[[i]] <- TE_df
  
  # Make TE vs lag plot ------------
  g <- lag_plots_all(TE_df,Site_ID)
  g_ls[[i]] <- g
}

# Output TE df list
saveRDS(TE_ls,paste0(Output_path,"/TE_results/TE_ls.rds"))
# Combine all TE lag figures
g_all <- plot_grid(plotlist = g_ls,
                   ncol=1)
print_g(g_all,paste0("/TE_results/TE_lag_all_sites"),
        14,4*13)



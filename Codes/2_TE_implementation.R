# Author: Zhaozhe Chen
# Update Date: 2026.7.6

# This codes conducts TE from P -> Q for selected hourly sites

# ------- Global -------
library(dplyr)
library(tidyr)
library(lubridate)
# Source functions
source("Functions.R")

# Data path ==========
Input_path <- "../../../Data/"
# Output path
Output_path <- "../Results/TE_results_sites/"

# --- Below variables need to be revised for each site --------
# This Site_ID is just for output file name
#Site_ID <- "HB_w3"
#Site_ID <- "HB_w8"
#Site_ID <- "Atlanta"
#Site_ID <- "Baltimore"
#Site_ID <- "Konza"
#Site_ID <- "Angelo_DRY_Q_Ts"
#Site_ID <- "Angelo_ELDER_Q_Ts"
#Site_ID <- "Angelo_DRY_Ta_Ts"
Site_ID <- "Angelo_ELDER_Ta_Ts"

# File specific var names
if(Site_ID == "HB_w3"){
  # Data file name
  filename <- "Hubbard_Brook_Alix/precip_discharge_w3.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DATETIME"
  # Precipitation variable name
  var1_name <- "precip"
  # Q variable name
  var2_name <- "avg_discharge"  
}else if (Site_ID == "HB_w8"){
  # Data file name
  filename <- "Hubbard_Brook_Alix/precip_discharge_w8.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DATETIME"
  # Precipitation variable name
  var1_name <- "precip"
  # Q variable name
  var2_name <- "avg_discharge"  
}else if (Site_ID == "Atlanta"){
  # Data file name
  filename <- "Atlanta/Atlanta_USGS_QP_0925.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "dateTime"
  # Precipitation variable name
  var1_name <- "Precip_Inst"
  # Q variable name
  var2_name <- "Flow_Inst"  
} else if (Site_ID == "Baltimore"){
  # Data file name
  filename <- "Baltimore/BES_USGS_QP_0917.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DateTime_EST_rnd"
  # Precipitation variable name
  var1_name <- "Precip_mm"
  # Q variable name
  var2_name <- "Q_cfs"  
} else if (Site_ID == "Konza"){
  # Data file name
  filename <- "Konza_Erin/Konza_df.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DateTime"
  # Precipitation variable name
  var1_name <- "P"
  # Q variable name
  var2_name <- "Q"  
} else if (Site_ID == "Angelo_DRY_Q_Ts"){
  # Data file name
  filename <- "Angelo_Laurel/CompiledAngeloSamplerPlatter_Hourly.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DateTime"
  # Precipitation variable name
  var1_name <- "DRY_discharge_cms"
  # Q variable name
  var2_name <- "DRY_temp_C"  
} else if (Site_ID == "Angelo_ELDER_Q_Ts"){
  # Data file name
  filename <- "Angelo_Laurel/CompiledAngeloSamplerPlatter_Hourly.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DateTime"
  # Precipitation variable name
  var1_name <- "ELDER_discharge_cms"
  # Q variable name
  var2_name <- "ELDER_temp_C"  
} else if (Site_ID == "Angelo_DRY_Ta_Ts"){
  # Data file name
  filename <- "Angelo_Laurel/CompiledAngeloSamplerPlatter_Hourly.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DateTime"
  # Precipitation variable name
  var1_name <- "DRY_air_temp"
  # Q variable name
  var2_name <- "DRY_temp_C"  
} else if (Site_ID == "Angelo_ELDER_Ta_Ts"){
  # Data file name
  filename <- "Angelo_Laurel/CompiledAngeloSamplerPlatter_Hourly.csv"
  # Variable names in the dataset ====
  # Time variable name
  time_name <- "DateTime"
  # Precipitation variable name
  var1_name <- "ELDER_air_temp"
  # Q variable name
  var2_name <- "ELDER_temp_C"  
}

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
         var1 = all_of(var1_name),
         var2 = all_of(var2_name)) %>%
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

# Make exploratory plots of the data ==========
# Time series of P and Q
g_P_TS <- TS_plot("var1",Site_df,Site_ID,my_color[3])+
  labs(y = var1_name)
g_Q_TS <- TS_plot("var2",Site_df,Site_ID,my_color[2])+
  labs(y = var2_name)
# Histogram of P and Q
g_P_hist <- plot_hist(Site_df,"var1",n_bin=11,zero_remove = FALSE,my_color[3])+
  labs(y = var1_name)
g_P_hist_no0 <- plot_hist(Site_df,"var1",n_bin=11,zero_remove = TRUE,my_color[3])+
  labs(y = var1_name)
g_Q_hist <- plot_hist(Site_df,"var2",n_bin = 11,zero_remove = FALSE,my_color[2])+
  labs(y = var2_name)
g_Q_hist_no0 <- plot_hist(Site_df,"var2",n_bin = 11,zero_remove = TRUE,my_color[2])+
  labs(y = var2_name)

# Combine these plots
g_data <- plot_grid(g_P_TS,g_P_hist,g_P_hist_no0,
                    g_Q_TS,g_Q_hist,g_Q_hist_no0,
                    nrow=2,align="hv",axis="btlr",
                    rel_widths = c(1,0.5,0.5))

# Implement hourly TE calculation ============
# Timing the TE calculation
start_time <- Sys.time()
# Run TE
TE_df <- Cal_TE_MI_main(Source = Site_df$var1,
                        Sink = Site_df$var2,
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
write.csv(TE_df,paste0(Output_path,"/TE_df_",Site_ID,".csv"))

# Combine all figures
g_all <- plot_grid(g_data,g_TE,nrow=2,
                   align="hv",axis="tblr",
                   rel_heights = c(2,1))

# Output TE plot
print_g(g_all,
        paste0("TE_lag_",Site_ID),
        16,12)



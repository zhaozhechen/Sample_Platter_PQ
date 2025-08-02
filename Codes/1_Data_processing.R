# Author: Zhaozhe Chen
# Date: 2025.8.1

# This codes processes Macroshed dataset for TE analysis

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
# Colors
my_color <- brewer.pal(5,"Dark2")


# ------- Main ---------
# Read in Macroshed dataset
MS_data <- read.csv(here(Input_path,"p_q_et_ms_data.csv"))
# Site info for tested sites
Site_info <- read.csv(here(Input_path,"refSites_for_PQanalysis.csv"))

# Get site names to be analyzed
Site_ID_ls <- Site_info$site_code
# Only keep data for these sites, and only keep required variables (P and Q)
MS_df <- MS_data %>%
  filter(site_code %in% Site_ID_ls,
         var == "discharge"|var == "precipitation")

# Get summary statistics of the dataset
summary_stat <- MS_df %>%
  group_by(site_code,var) %>%
  summarise(
    n_total = n(),
    n_status1 = sum(ms_status == 1,na.rm=TRUE),
    n_interp1 = sum(ms_interp == 1,na.rm=TRUE)
  ) %>%
  mutate(var = recode(var,discharge ="Q",precipitation = "P")) %>%
  pivot_wider(id_cols = site_code,
              names_from = var,
              values_from = c(n_total,n_status1,n_interp1),
              names_glue = "{var}_{.value}") %>%
  mutate(Q_stats1_pt = round(Q_n_status1/Q_n_total,3)*100,
         Q_interp1_pt = round(Q_n_interp1/Q_n_total,3)*100,
         P_stats1_pt = round(P_n_status1/P_n_total,3)*100,
         P_interp1_pt = round(P_n_interp1/P_n_total,3)*100) %>%
  select(site_code,Q_n_total,Q_stats1_pt,Q_interp1_pt,
         P_n_total,P_stats1_pt,P_interp1_pt)

# Initialize vectors to store interpolation info
interp_n <- c()
interp_interval_mean <- c()
interp_interval_median <- c()
# Initiliaze a list to store output figures
g_all <- list()
# Loop over all sites
for(i in 1:length(Site_ID_ls)){
  Site_ID <- Site_ID_ls[i]
  # Get data for this site
  Site_df <- MS_df %>%
    filter(site_code == Site_ID)
  # Remove questionable data
  Site_df$val[!is.na(Site_df$ms_status) & Site_df$ms_status == 1] <- NA
  
  # Get length distribution of interpolated data
  # Note: since no interpolated Q, so only did this for precipitation
  interp_lengths <- interpolated_interval("precipitation",Site_df)
  # Get the number of interpolated segments
  if(length(interp_lengths)==1){
    if(is.na(interp_lengths)){
      interp_n <- c(interp_n,0) 
    }
  }else{
    interp_n <- c(interp_n,length(interp_lengths))
  }
  # Get the mean interpolation interval
  interp_interval_mean <- c(interp_interval_mean,mean(interp_lengths))
  # Get the median interpolation interval
  interp_interval_median <- c(interp_interval_median,median(interp_lengths))
  
  # Make a TS plot of Q and P
  g_P <- TS_plot("precipitation",Site_df,Site_ID)
  g_Q <- TS_plot("discharge",Site_df,Site_ID)
  
  # Get the distribution of P and Q
  g_hist_P <- Hist_plot("precipitation",Site_df)
  g_hist_Q <- Hist_plot("discharge",Site_df)

  # Distribution after removing 0
  
  # Make the distribution of precipitation interpolation length
  
  # Combine them together
  
  g <- plot_grid(g_Q,g_hist_Q,
                 g_P,g_hist_P,
                 nrow=2,align="hv",axis="tblr")
  g_all[[i]] <- g
}

# Put interpolation info into a df
# Note: since no interpolated Q, so only did this for precipitation
interp_summary_df <- data.frame(
  site_code = Site_ID_ls,
  interp_n = interp_n,
  interp_interval_mean = interp_interval_mean,
  interp_interval_median = interp_interval_median
)
# Combine this with the summary_stat df
summary_stat <- summary_stat %>%
  right_join(interp_summary_df,by="site_code")
# Output this summary table
write.csv(summary_stat,here(Output_path,"Site_summary.csv"))



# Author: Zhaozhe Chen
# Date: 2025.10.15

# This code is to explore site data


# ------- Global ----------

library(dplyr)

Bal_df <- read.csv("D:/OneDrive - UW-Madison/Research/ET Synchrony/CZ_Sync/Sampler platter/Data/Atlanta/Altanta_USGS_QP_0925.csv")


test <- Bal_df %>%
  filter(site_no == "2335350")
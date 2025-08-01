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

# Get interpolated value intervals



# Output this summary table
write.csv(summary_stat,here(Output_path,"Site_summary.csv"))

# Loop over all sites
for(i in length(Site_ID_ls)){
  Site_ID <- Site_ID_ls[i]
  # Get data for this site
  Site_df <- MS_df %>%
    filter(site_code == Site_ID)


}








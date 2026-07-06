# Author: Zhaozhe Chen
# Update Date: 2026.7.6

# This codes pre-process df for TE implementation for each site

# ------- Global -------
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(lubridate)

# --- Main --------
# Konza ============
P_folder <- "D:/OneDrive - UW-Madison/Research/ET Synchrony/CZ_Sync/Sampler platter/Data/Konza_Erin/precip_data/precip_data/"

files <- list.files(P_folder,
                    pattern = "\\.txt$",
                    full.names = TRUE)

P_df <- map_dfr(files, function(f){
  
  dat <- read.table(
    f,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Keep only rows with 3 columns
  dat <- dat[,1:3]
  names(dat) <- c("ID","Time","P")
  
  # Remove bad rows
  dat <- dat %>%
    filter(!is.na(ID),
           !is.na(Time),
           !is.na(P))
  
  # Force precipitation to numeric
  dat$P <- as.numeric(dat$P)
  
  # Remove rows that failed conversion
  dat <- filter(dat, !is.na(P))
  
  date <- as.Date(substr(dat$ID,
                         nchar(dat$ID)-5,
                         nchar(dat$ID)),
                  "%y%m%d")
  
  time <- sprintf("%04d", dat$Time)
  
  tibble(
    DateTime = as.POSIXct(
      paste(
        format(date,"%m/%d/%Y"),
        paste0(substr(time,1,2),":",substr(time,3,4))
      ),
      format="%m/%d/%Y %H:%M"
    ),
    P = dat$P
  )
})

# Q_df
Q_df <- read.csv("D:/OneDrive - UW-Madison/Research/ET Synchrony/CZ_Sync/Sampler platter/Data/Konza_Erin/discharge data/discharge data/ASD021.csv")
Q_df <- Q_df %>%
  select(RecYear,RecMonth,Recday,RecHour,Discharge) %>%
  mutate(
    RecHour = sprintf("%04d", RecHour),
    DateTime = as.POSIXct(
      paste(
        sprintf("%02d/%02d/%04d", RecMonth, Recday, RecYear),
        paste0(substr(RecHour, 1, 2), ":", substr(RecHour, 3, 4))
      ),
      format = "%m/%d/%Y %H:%M"
    )
  ) %>%
  transmute(
    DateTime,
    Q = Discharge
  )

# Combine P and Q
Konza_df <- inner_join(P_df, Q_df, by = "DateTime") %>%
  select(DateTime, P, Q) %>%
  arrange(DateTime)
write.csv(Konza_df,"D:/OneDrive - UW-Madison/Research/ET Synchrony/CZ_Sync/Sampler platter/Data/Konza_Erin/Konza_df.csv")


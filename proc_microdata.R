#################################################################################
# proc_microdata.r
# Microclimate data input and processing for Athens Mosquito Microclimate project
# written by Mike Wimberly, University of Oklahoma
#################################################################################

library(tidyverse)
library(lubridate)

# Read in the raw data from the Google Drive archive
raw_data_path <- "c:/Users/wimb0002/Google Drive/DataLoggerData/2018_rawdata/"

# Function to read in the microclimate logger files
my_read_csv <- function(x) {
  #out <- read_csv(x, col_types = cols(.default = "c"))
  out <- read_csv(x)
  if ( nrow(out) >= 1 ) {
    site <- substr(x, str_length(x) - 7, str_length(x) - 4)
    cbind(Site=site, out)
  }
}

# Read all logger files into a single R object
micro_in <-
  list.files(path = raw_data_path,
             pattern = "*.csv", 
             recursive = TRUE,
             full.names = T) %>% 
  map_df(~my_read_csv(.)) %>%
  separate(Time, into = c("Date", "Time"), sep = " ") %>%
  mutate(Date = mdy(Date))

names(micro_in)[c(4, 7)] <- c("temp", "rh")

# Remote bad sites
# R3A2, S1A4, and U34 have malfunctioning loggers
# other values are bad site names
allsites <- unique(micro_in$Site)
badsites <- c("R3A2", "S1A4", "U3A4", "2A12", "6UNK", "0820")
goodsites <- allsites[!allsites %in% badsites]

# Generate daily summaries of mean, min, and max RH and temp values
micro_daily <- micro_in %>%
  group_by(Date, Site) %>%
  filter(Site %in% goodsites) %>%
  summarise(num_temp = sum(!is.na(temp)),
            num_rh = sum(!is.na(rh)),
            mean_temp = mean(temp),
            min_temp = mean(temp[temp <= quantile(temp, 0.04)]),
            max_temp = mean(temp[temp >= quantile(temp, 0.96)]),
            mean_rh = mean(rh),
            min_rh = mean(rh[rh <= quantile(rh, 0.04)]),
            max_rh = mean(rh[rh >= quantile(rh, 0.96)])) %>%
  ungroup() %>%
  arrange(Date, Site)

# Crosstabulate sites by date
micro_xtab <- micro_daily %>%
  dplyr::select(Date, Site, num_temp) %>%
  spread(Site, num_temp)

# Save data
# Raw data
write_csv(micro_in, "./outdata/athens_2018_rawdata.csv")
# Summarized daily data
write_csv(micro_daily, "./outdata/micro_2018_daily.csv")
# Cross-tabulation of daily data by logger
write_csv(micro_xtab, "./outdata/micro_2018_xtab.csv")

# Remote temporary objects
rm(list = ls())

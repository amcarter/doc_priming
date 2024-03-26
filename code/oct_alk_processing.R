#Fall alkalinity processing

library(tidyverse)
library(dplyr)

dat <- read_csv('data/alkalinity/fl_ybc_oct_alkalinity.csv')

dat_filt <- dat[!grepl("\\S", dat$notes),]

avgs <- dat_filt %>%
    group_by(site) %>%
    summarize(avg_alk_meqL = mean(alk_meqL, na.rm = TRUE))
print(avgs)

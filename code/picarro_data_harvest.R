# Picarro data handling attempt.
# Stealing from https://github.com/bpbond/R-data-picarro/blob/master/scripts/picarro.R

# Updated summer 2023 by A. Carter
library(plyr)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(dygraphs)

# Set your working directory if needed
# setwd("C:/Users/julia/Dropbox/Kirby_priming/Data_Processing/Data_Harvest_By_Day/Kirby_trial_run_2023-20230711T182315Z-001/Kirby_trial_run_2023")

###########
## IMPORT DATA
###########
# Put in paths to data files and folders here.
# This should be the only place in the script that needs to be changed when used.
# Data folder where all of the raw data is located including run files from the
# Picarro, the txt files you generated while running it, and any associated metadata:
dat_folder <- "data/Kirby_NyackW/NyackW_08_04_2023/"
run_file <- "NyackW_08_04_2023.txt"
# meta_file <- "" # optional additional Metadata
output_file <- "nw_2_CO2_processed.csv" # provide a name for your output file
batch <- 2

## Import meta data for each day:  In this case metadata consists of a sample
## identifier and time stamp, where time stamp is when the sample syringe hits
## about 10 mL and it looks like the trace is leveling off for that sample.
## This metadata can and should have other stuff in it, such as sample location,
## collection time, pH, volume of sample, volume of equlibration gas (typically 70mL),
## temp of equilibration

meta <- read.table(paste0(dat_folder, run_file),
                   sep=",",header=T)
meta$DATETIME <- lubridate::ymd_hms(as.character(meta$ResetTime))
meta$samplenum<-seq(1:nrow(meta))


## Import the raw data
filelist <-list.files(dat_folder, pattern = "\\.dat", full.names = TRUE)  ##makes a vector of that day's files
rawdata <- data.frame()  #empty dataframe
#below loops over all the files making a list containg a tibble for each file and stuffing it into a list
for(f in filelist) {
  cat("Reading", f, "\n") # prints which file is being read
  dd <- read_table(f) %>%
    select(DATE, TIME, ALARM_STATUS, MPVPosition, #add whichever things we want here
           `12CO2_dry`, #dry mole fraction of 12CO2, corrected for water vapor
           `13CO2`, Delta_Raw_iCO2,
           HP_12CH4, `HP_Delta_iCH4_Raw`) # methane data, if needed.
  rawdata <- bind_rows(rawdata, dd)
}

## Fix datetime
rawdata <-
    mutate(rawdata,
           DATETIME = lubridate::ymd_hms(paste(DATE, TIME))) %>%
    select(-DATE, -TIME) # Get rid of unneeded fields

## average data for each second where we have a measure.  Not optional
rawdata <- rawdata %>%
    group_by(DATETIME) %>%
    summarise_all(mean)

# select the integration window for each sample (in seconds):
## 40 seconds appears correct, but always check this time.
## To short and we lost data.  Too long and we average off the plateau
integration_window = 40

#Plot it. Look ok?

sample_times <- as.POSIXct(c(), tz = 'UTC')
for(i in 1:nrow(meta)){
    dd <- seq(meta$DATETIME[i], meta$DATETIME[i] + integration_window, by = 'sec')
    sample_times <- c(sample_times, dd)
}
sample_times <- data.frame(DATETIME = lubridate::ymd_hms(sample_times)) %>%
    left_join(rawdata[,c(1,4)]) %>%
    dplyr::rename(sample_window = `12CO2_dry`)


dygraphs::dygraph(full_join(rawdata[,c(1,4,6)], sample_times)) %>%
    dySeries("sample_window", drawPoints = TRUE, color = 'red') %>%
    dyRangeSelector()

# plot(rawdata$DATETIME,rawdata$`12CO2_dry`)
# plot(rawdata$DATETIME,rawdata$Delta_Raw_iCO2)
# plot(rawdata$DATETIME,rawdata$HP_12CH4)
# plot(rawdata$DATETIME,rawdata$`HP_Delta_iCH4_Raw`)

##get rid of all the early in the day data to make a smaller file.  This step is optional
#rawdata <- subset(rawdata, rawdata$DATETIME > "2022-07-08 17:30:00")

# Check that there weren't any wierd alarms:
unique(rawdata$ALARM_STATUS) # Don't know what these mean, maybe we should look into them but not rn.

##if ok, then name the rawdat file for that day here
dat <- rawdata %>%
    mutate(CO2_ppm = `12CO2_dry` + `13CO2`) %>%
    rename(Delta_iCO2 = Delta_Raw_iCO2, CH4 = HP_12CH4,
           Delta_iCH4 = HP_Delta_iCH4_Raw) %>%
    select(-ALARM_STATUS, -MPVPosition, -`12CO2_dry`, -`13CO2`)



#############
## SUBSET
#############

## Subset the main datafile based on sample time + integration window,
meta$ENDTIME <- meta$DATETIME + integration_window

subdat<-list()
for(i in 1:length(meta$samplenum)){
    subdat[[i]] <- subset(dat, DATETIME>meta$DATETIME[i]  & DATETIME< meta$ENDTIME[i])
    subdat[[i]]$samplenum <- i
}
subdat <- bind_rows(subdat)

## Plot
# do the sample intervals look relatively flat? If not, adjust your sampling window.
subdat_melt <- subdat %>%
    melt(id.vars=c("DATETIME","samplenum"))
ggplot(subdat_melt, aes(DATETIME, value, color = factor(samplenum)))+
  geom_point()+
  facet_wrap(~variable, ncol=2, scales = "free")


## Summarize per sample
subdatsum <- subdat %>%
    group_by(samplenum) %>%
    summarise(CO2_ppm = mean(CO2_ppm), delCO2=mean(Delta_iCO2),
              CH4=mean(CH4), delCH4=mean(Delta_iCH4))


##merge with meta, and give it a new name

# The metafile should contain: water_temp, equil_temp, pH, ANC,
# air_vol (it won't always be 70), baro_press (absolute, mm Hg)
# not the  kH_CO2 because that is a derived value that we will calc each time

sum_file <- full_join(meta, subdatsum, by = 'samplenum')

sum_file$batch <- batch
# name your file something informative and save it in the same folder.
write_csv(sum_file, paste0(dat_folder, output_file))


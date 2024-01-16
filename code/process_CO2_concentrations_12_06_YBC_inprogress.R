# Process CO2 data to:
# 1. Correct concentrations based on standards if desired
# 2. Calculate the concentration in H2O before equilibration
# Run this code after you have run the picarro_data_harvest.R
# code on all picarro runs in a batch

# setup ####

library(tidyverse)
library(lubridate)
# setwd("/Users/niko/Dropbox/leachate_linder")

# data ####
# read in sample metadata file for the entire duration of experiment
# This will include columns the following columns at a minimum:
# Syringe: The syringe number, matching what was entered when you ran the Picarro
# Treatment: What teatment is each sample from?
# rundate: the date that samples were run on the picarro (in picarro dates!)
# sample_date: the date that samples were taken
# sample_time: the actual times that samples were taken
# equilpressure: the air pressure in the lab when the samples were equilibrated (in units of atm)
dat_folder <- "data/Kirby_YBC/" # the folder that contains subfolders of processed Picarro run data
mdat <- read_csv(paste0(dat_folder, "YBC_meta.csv")) %>%
    mutate(sample_date = as.Date(sample_date, format = '%m/%d/%Y'))
expt_name <- 'ybc_expt' # used for generating figures and saving output.

# functions ####
## CALCULATE CO2 IN H20 (in the post-processing files now)
#I (Bob) think this is correct
CO2calc2<-function(x) {
    x$kH_CO2<-29* exp( 2400 * ((1/(x$equiltemp+273.15)) - 1/(298.15)) ) #temp correction for Henry's law
    x$CO2_mol_water<- x$CO2_conc_air  / x$kH_CO2
    x$CO2_mol_air<- x$CO2_conc_air/ ((x$equiltemp+273.15)*0.08026/x$equilpressure)

    x$CO2_conc_umolL<-(x$CO2_mol_water*x$volwater + x$CO2_mol_air*x$volair)/x$volwater # µmol/L or mmol/m3
    return(x)

}

# convert del value to atomic fraction and back:
deltoAF<- function (del) {
    Rst<-0.0112372
    R<-(del/1000+1)*Rst
    AF<-R/(1+R)

    return(AF)
}

AFtodel <- function (AF) {
    Rst<-0.0112372
    R<- AF/(1-AF)
    del<- (R/Rst-1)*1000

    return(del)
}

# standard_correct <- function(std_CO2, std_d13C, dat){
#
#     stds <- filter(dat, Syringe == 'Std_1')
#     CO2 <- mean(stds$CO2, na.rm = T)
#     dat$CO2_norm <- dat$CO2 * std_CO2/CO2
#
#     # to correct for the del13C value, convert to AF and conc of 13C
#     del13 <- mean(stds$delCO2, na.rm = T)
#     CO2_13 <- CO2 * deltoAF(del13)
#     std_CO2_13 <- std_CO2 * deltoAF(std_d13C)
#     dd <- tibble(AF = deltoAF(dat$delCO2),
#                  CO2_13 = dat$CO2 * AF)
#     dd$CO2_13_norm <- dd$CO2_13 * std_CO2_13/CO2_13
#     dd$del13_norm <- AFtodel(dd$CO2_13_norm/dat$CO2_norm)
#
#     dat$delCO2_norm = dd$del13_norm
#
#     return(dat)
# }


## Compile data files ####
# correct each file based on the standard concentration
# caclucate CO2 concentration in the H2O before equilibration

# standard 1 known concentrations:
# std_CO2 = 435.5 # ppm
# std_d13C = -2.25 # del13C

dat <- tibble() # create empty data frame to compile samples
filelist <- list.files(dat_folder, pattern = 'processed\\.csv',
                       recursive = TRUE)

for(f in filelist){ # compile all processed data from experiment into one table

    dd <- read_csv(paste0(dat_folder, f)) %>%
        mutate(rundate = as.Date(DATETIME))

    # dd <- standard_correct(std_CO2, std_d13C, dd)

    dd <- dd %>%
      mutate(CO2_conc_air = CO2_ppm)  %>%
        dplyr::rename(equiltemp = eq_temp,
                      volair = air_mL,
                      volwater = H2O_mL)


    dat <- bind_rows(dat, dd)
}

# dat <- filter(dat, Syringe != '70')
# mdat <- filter(mdat, Syringe != '70')

# add metadata to the compiled list of data and clean up:
dat <- left_join(mdat, dat, by = c('batch', 'Syringe')) %>%
    mutate(sample_datetime = lubridate::ymd_hms(paste(sample_date, sample_time, sep = " ")),
           Treatment = as.factor(Treatment)) %>%
    select(-ResetTime, -ENDTIME, -sample_date, -sample_time) %>%
    rename(run_datetime = DATETIME)

## Run function on datafile
sum_file <- CO2calc2(dat)

sum_file %>% group_by(Treatment, batch) %>%
    summarize(n = n()) %>% arrange(batch)

ggplot(sum_file, aes(batch, CO2_conc_umolL, col = Treatment))+
    geom_point(size=3)+
    theme(axis.text.x=element_text(angle=45, hjust = 1))

dat <- sum_file %>%
    mutate(AF = deltoAF(delCO2),
           CO2_13C_umolL = AF * CO2_conc_umolL)

## plot the data

ggplot(dat, aes(x = sample_datetime, y = CO2_conc_umolL, col = Treatment))+
    geom_point() + geom_line()

ggplot(dat, aes(x = sample_datetime, y = CO2_13C_umolL, col = Treatment))+
  geom_point()+
    geom_line()

# dat %>% filter(batch == 3) %>%
#     ggplot(aes(x = sample_datetime, y = CO2_conc_umolL, col = Treatment))+
#     geom_
#wasn't working (JDK 8/1)

dat %>% filter(Treatment==7, batch==2, replicate=="C")

dat %>%
    group_by(batch, Treatment) %>%
    summarize(across(.cols = c('CO2_conc_umolL', 'delCO2'),
                     .fns = list(mean, sd)))

# ggplot( aes(x = batch, y = delCO2, col = Treatment))+
#   geom_boxplot()
#wasn't working (JDK 8/1)

# example code to save a figure:
# png(filename = paste0('figures/', expt_name, '_13CO2.png')) # naming it based on expt name will prevent us from overwriting it
# ggplot(dat, aes(x = sample_datetime, y = CO2_13C_umolL, col = Treatment))+
#     geom_point()
# dev.off()

write_csv(dat, paste0(dat_folder, expt_name, '.csv'))


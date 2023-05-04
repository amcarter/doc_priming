# Process CO2 data to:
# 1. Correct concentrations based on standards
# 2. Calculate the concentration in H2O before equilibration

# setup ####
library(tidyverse)
library(lubridate)
setwd("/Users/niko/Dropbox/leachate_linder")

# data ####
# read in sample metadata file for the entire duration of experiment
mdat <- read_csv("Sampletimes2.csv")
# functions ####
## CALCULATE CO2 IN H20 (in the post-processing files now)
#I (Bob) think this is correct
CO2calc2<-function(x) {
    x$kH_CO2<-29* exp( 2400 * ((1/(x$equiltemp+273.15)) - 1/(298.15)) ) #temp correction for Henry's law
    x$CO2_mol_water<- x$CO2_vol_air  / x$kH_CO2
    x$CO2_mol_air<- x$CO2_vol_air/ ((x$equiltemp+273.15)*0.08026/x$pressure)

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
filelist <- list.files('processed_CO2/expone_7_12_22')
t = 0
for(f in filelist){ # compile all processed data from experiment into one table

    dd <- read_csv(paste0("processed_CO2/expone_7_12_22/", f)) %>%
        mutate(rundate = as.Date(DATETIME),
               pressure = 0.9,
               sample_interval = paste0("T", t))
    t = t + 1

    # dd <- standard_correct(std_CO2, std_d13C, dd)


    dd <- dd %>%
      mutate(CO2_vol_air = CO2)  %>%
        dplyr::rename(SampleID = ...1,
               equiltemp = eq_temp,
               volair = air_mL,
               volwater = H2O_mL)

    ## Run function on datafile
    sum_file <- CO2calc2(dd)

    ggplot(sum_file, aes(SampleID, CO2_conc_umolL))+
        geom_point(size=3)+
        theme(axis.text.x=element_text(angle=45, hjust = 1))

    dat <- bind_rows(dat, sum_file)
}

## summarize data so far
# calculate atomic fractions:
dat <- dat %>%
    dplyr::rename(time_interval = sample_interval,
           syringe_number = Syringe)

dat <-left_join(mdat,dat, by = c('rundate', 'time_interval', 'syringe_number'))

ggplot(dat, aes(x = sample_datetime, y = CO2_conc_umolL, col = treatment))+
    geom_point()

ggplot(dat, aes(x = sample_datetime, y = delCO2, col = treatment))+
  geom_point()

dat <- dat %>%
    mutate(AF = deltoAF(delCO2),
           CO2_13C_umolL = AF * CO2_conc_umolL)

ggplot(dat, aes(x = sample_datetime, y = CO2_13C_umolL, col = treatment))+
    geom_point()


write_csv(dat,"processed_CO2/experiment1_summary2.csv")


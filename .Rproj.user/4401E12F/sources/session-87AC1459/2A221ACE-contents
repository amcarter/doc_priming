library(tidyverse)
library(lubridate)
library(dplyr)


# setup ####
setwd("/Users/niko/Dropbox/leachate_linder/processed_CO2")
# dat <- read_csv("~/Dropbox/leachate_linder/processed_CO2/experiment1_summary2.csv")
dat <- read_csv("processed_CO2/experiment1_summary2.csv")
glimpse(dat)

deltoAF<- function (del) {
  Rst<-0.0112372
  R<-(del/1000+1)*Rst
  AF<-R/(1+R)
  AF

}


AFtodel <- function (AF) {
  Rst<-0.0112372
  R<- AF/(1-AF)
  del<- (R/Rst-1)*1000

  del
}

# calculate DIC in the water ####
###testing DIC from pCO2-alk
###partial pressures of CO2 pCO2 in uatm
alk <- 1.588 *10 ^(-3) # mol/L
pCO2 <- 30 * 10^(-6) # mol/L
temp <- 273+21
K1 <- exp(290.9097-14554.21/(temp)-45.0575*log(temp))
K2 <- exp(207.6548-11843.79/(temp)-33.6485*log(temp))


DICfrom_C_A <- function(K1, K2, C, A){
  H_minus <- (((-K1*C))-sqrt(((K1*C)^2)-(4*-1*A*2*K1*K2*C)))/(2*-1*A)
  pH <- -1*log10(H_minus)
  B <- (K1*C)/H_minus
  Ca <- (K2*B)/H_minus
  D <- C + B + Ca
  A_check <- B + 2*Ca

  l <- list(H_minus, pH, C, B, Ca, D, A, A_check)
  names(l) <- c("H", "pH", "C", "B", "Ca", "D", "A","A check")
  DIC = l$D
  return(DIC)
}

out <-DICfrom_C_A(K1=K1, K2=K2,C=pCO2,A=alk)


dat$DIC_molL <- DICfrom_C_A(K1, K2, dat$CO2_conc_umolL * 10^(-6), alk)
dat$DIC_mgL <- dat$DIC_molL*12.011*1000
dat$AF = deltoAF(dat$delCO2)
dat$DIC_13C_umolL = dat$DIC_molL * dat$AF *10^6

dat %>% dplyr::filter(treatment != 'LW_nolechate') %>%
    filter(!is.na(treatment)) %>%
ggplot(aes(x = sample_datetime, y = DIC_molL, col = treatment))+
  geom_point()
View(dat)

DIC_dataframe <- data.frame(
  time_interval = dat$time_interval,
  treatment = dat$treatment,
  DIC_13C_umolL = dat$DIC_13C_umolL
)

DIC_dataframe <- dat %>%
    select(-volwater, -volair, -equiltemp, -ResetTime, -DATETIME, -samplenum,
           -ENDTIME, -CH4, -delCH4, -pressure, -CO2_mol_water, -kH_CO2, -CO2_vol_air,
           -CO2_mol_air) %>%
    rename(raw_CO2 = CO2)


# dat <- dat %>%
#   mutate(add_column(dat, ifelse(treatment = C,CNP)))
write_csv(DIC_dataframe, "processed_CO2/experiment1_processed.csv")
DIC <- read_csv('processed_CO2/experiment1_processed.csv')

# calculate excess 13CDIC in each treatment ####
sumDIC <- filter(DIC, time_interval %in% c('T0', 'T5')) %>%
    select(time_interval, treatment, sample_datetime, DIC_molL, DIC_13C_umolL) %>%
    group_by(time_interval, treatment) %>%
    # summarize(across(.fns = list(mean = ~mean(.), sd = sd(.))))
    summarize(sample_datetime = mean(sample_datetime),
              DIC_molL = mean(DIC_molL),
              DIC_13C_sd = sd(DIC_13C_umolL),
              DIC_13C_umolL = mean(DIC_13C_umolL)) %>%
    ungroup()

DIC_table <- sumDIC %>%
    mutate(excess_13C_DIC_umolL = DIC_13C_umolL - sumDIC$DIC_13C_umolL[sumDIC$time_interval == 'T0'],
           excess_DIC_umolL = (DIC_molL - sumDIC$DIC_molL[sumDIC$time_interval == 'T0'])*10^6,
           excess_13C_DIC_umolL_sd = DIC_13C_sd + sumDIC$DIC_13C_sd[sumDIC$time_interval == 'T0'])


#creating excess accumulated 13C CO2 (from DIC) plot  ####
# setwd('/Users/niko/Dropbox/leachate_linder/processed_CO2')
# DIC_table <- read_csv('~/Dropbox/leachate_linder/processed_CO2/DIC_table.csv')
DIC_table$nutrients <- c(0,0,1,0,1)
DIC_table$sugar <- c(0,1,1,0,0)
# DIC_table$DICXS <- DIC_table$DICXS*1000
# DIC_table$STDev_T5 <- DIC_table$STDev_T5*1000
# DIC_table$STDev_T0 <- DIC_table$STDev_T0*1000
# DIC_table$STDev <- DIC_table$STDev_T0 + DIC_table$STDev_T5

DIC_table <- slice(DIC_table, -1)

plot(DIC_table$sugar[DIC_table$nutrients == 1],
     DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1], type = 'b',
     xlim = c(-0.25,1.25), ylim = c(0,1.5), pch = 15, col = 'blue',
     xaxt = 'n', xlab = 'No Sugar                                         Sugar',
     ylab = 'Accumulated Excess 13C CO2 (µg)')


arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 1],
       y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1],
       x1 = DIC_table$sugar[DIC_table$nutrients == 1],
       y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1]+
         DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==1], length = 0.1, angle = 90)
arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 1],
      y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1],
      x1 = DIC_table$sugar[DIC_table$nutrients == 1],
      y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1]-
        DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==1], angle = 90, length = 0.1)

lines(DIC_table$sugar[DIC_table$nutrients == 0],
      DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0],pch = 1, col = 'forestgreen')
points(DIC_table$sugar[DIC_table$nutrients == 0],
       DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0], pch = 15, col = 'forestgreen')
arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 0],
       y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0],
       x1 = DIC_table$sugar[DIC_table$nutrients == 0],
       y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0]+
         DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==0], length = 0.1, angle = 90)
arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 0],
       y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0],
       x1 = DIC_table$sugar[DIC_table$nutrients == 0],
       y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0]-
         DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==0], angle = 90, length = 0.1)

legend('topleft', legend = c('Nutrients', 'No Nutrients'),
       col = c('green', 'blue'), lty = c(1,1), pch = c(15, 15), border = NA)

# plot(DIC_table$sample_datetime, dat$DIC_13C, type = 'b', pch = 15)
lines(dat$sample_datetime)

arrows()

library(tidyverse)
library(lubridate)
library(dplyr)


# setup ####
dat <- read_csv("data/Kirby_YBC/ybc_expt.csv")
#alk <- read_csv('data/alkalinity/nw_alkalinity.csv') %>%
 #   filter(!is.na(trt)) %>%
  #  mutate(date = as.Date(date, format = '%m/%d/%Y')) %>%
   # select(-samplenum, -treatment, -vol_acid)
glimpse(dat)


# Functions ####
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

# metadata ####
trts <- data.frame(Treatment = 1:8,
                   trt = c('W',  'NP', 'C', 'CNP','W',  'NP', 'C', 'CNP'),
                   leachate = c(0,0,0,0,1,1,1,1))

# add actual data for this part once we have DOC from Adam
c_source <- read_csv('data/DOC_stock_concentrations.csv')

# There is no clear pattern in the alkalinity data. Take an average:
#alk <- alk %>% filter(final_pH > 4.3) %>%
 #   summarize(alk_meqL = mean(alk_meqL))

dat <- left_join(dat, trts, by = 'Treatment') %>%
    mutate(date = as.Date(sample_datetime))
    #left_join(alk, by = c('batch', 'trt', 'leachate'))

# calculate DIC in the water ####
###testing DIC from pCO2-alk
###partial pressures of CO2 pCO2 in uatm
alk <- 1.51 # meq/L
# pCO2 <- 30 * 10^(-6) # mol/L
temp <- 273+21
K1 <- exp(290.9097-14554.21/(temp)-45.0575*log(temp))
K2 <- exp(207.6548-11843.79/(temp)-33.6485*log(temp))

# out <-DICfrom_C_A(K1=K1, K2=K2,C=pCO2,A=alk)
dat$DIC_molL <- DICfrom_C_A(K1, K2, dat$CO2_conc_umolL * 10^(-6),
                            alk* 10^-3)
# dat$DIC_13molL <- DICfrom_C_A(K1, K2, dat$CO2_13C_umolL * 10^(-6), alk)
dat$DIC_mgL <- dat$DIC_molL*12.011*1000
# dat$AF = deltoAF(dat$delCO2)
dat$DIC_13C_umolL = dat$DIC_molL * dat$AF *10^6
dat$DIC_13C_ugL = dat$DIC_13C_umolL * 13


dat %>% dplyr::filter(!is.na(Treatment)) %>%
ggplot(aes(x = sample_datetime, y = DIC_13C_ugL))+
  geom_point() + geom_line()

# # Clean up data ####

#dat <- dat %>%
 #   filter(!(batch == 2 & Syringe == '28'),
  #         !is.na(Treatment))

# # Remove high time zero point:
# dat %>% filter(batch == 1, Treatment == 6)
# # It looks like the high one is syringe 8G. Check this in the notes.
#
# dat <- dat %>%
#     filter(!(batch == 1 & Syringe == '8G'),
#            !is.na(Treatment))

hist(dat$delCO2)


# START OF FINAL FIGURES CODE####
# filter to use only the start and end time points
dat_sum <- dat %>%
    group_by(batch, trt, leachate, replicate) %>%
    summarize(across(any_of(c('sample_datetime', 'delCO2', 'DIC_mgL',
                              'DIC_13C_ugL', 'alk_meqL')),
                     .fns = ~mean(., na.rm = T))) %>%
    ungroup()

# Calculate the change over time for each bottle measured at time T
bottle_dif <-
    dat_sum %>% select(-starts_with(c('delCO2', 'alk_meqL'))) %>%
    mutate(batch = case_when(batch == 1 ~ 'T0',
                             batch == 2 ~ 'T1'))%>%
    pivot_wider(names_from = c('batch', 'replicate'),
                values_from = c('sample_datetime', 'DIC_mgL', 'DIC_13C_ugL')) %>%
    rename_with(~sub('_NA$', '', .)) %>%
    mutate(across(starts_with('DIC_mgL_T1'), .fns = ~ . - DIC_mgL_T0),
           across(starts_with('DIC_13C_ugL_T1'), .fns = ~ . - DIC_13C_ugL_T0),
           across(starts_with('sample_datetime_T1'), .fns = ~ . - sample_datetime_T0)) %>%
    select(-ends_with('T0')) %>%
    pivot_longer(cols = ends_with(c('_A', '_B', '_C', '_D')),
                 names_to = c('.value', 'rep'),
                 names_sep = '_T1_') %>%
    mutate(sample_datetime = as.numeric(sample_datetime),
           DIC_ugLd = DIC_mgL/sample_datetime*1000,
           DIC_13C_ugLd = DIC_13C_ugL/sample_datetime) %>%
    select(-sample_datetime, -DIC_mgL, -DIC_13C_ugL) %>%
    mutate(site = "ybc")

ggplot(bottle_dif, aes(trt, DIC_13C_ugLd, fill = factor(leachate))) +
    geom_violin()
write_csv(bottle_dif, 'data/DIC_change_by_bottle_ybc.csv')

dat_sum <- dat_sum %>%
    group_by(batch, trt, leachate) %>%
    summarize(across(any_of(c('sample_datetime', 'delCO2', 'DIC_mgL',
                              'DIC_13C_ugL', 'alk_meqL')),
                     .fns = list(mean = ~mean(., na.rm = T),
                                 sd = ~sd(., na.rm = T)))) %>%
    rename(sample_datetime = sample_datetime_mean) %>%
    select(-sample_datetime_sd) %>%
    mutate(across(ends_with('sd'), ~case_when(is.na(.) ~ 0,
                                              TRUE ~ .))) %>%
    mutate(site = "ybc")
write_csv(dat_sum, 'data/DIC_overtime_ybc.csv')

p12 <- ggplot(dat_sum, aes(sample_datetime, DIC_mgL_mean, col = trt, lty = factor(leachate)))+
    geom_point() + geom_line(size = 1) +
    geom_errorbar(aes(ymin = DIC_mgL_mean - DIC_mgL_sd,
                      ymax = DIC_mgL_mean + DIC_mgL_sd), lty = 1)+
    xlab('Date') +theme_bw()
p13 <- ggplot(dat_sum, aes(sample_datetime, DIC_13C_ugL_mean, col = trt, lty = factor(leachate)))+
    geom_point() + geom_line(size = 1) +
    geom_errorbar(aes(ymin = DIC_13C_ugL_mean - DIC_13C_ugL_sd,
                      ymax = DIC_13C_ugL_mean + DIC_13C_ugL_sd), lty = 1)+
    xlab('Date') +theme_bw()

png(filename = 'figures/ybc_expt_13CO2overtime.png',
    width = 8, height = 5, res = 300, units = 'in')
ggpubr::ggarrange(p12, p13, common.legend = T)
dev.off()

#removing baseline ####

dat_dif <- dat_sum %>%
    mutate(AF = deltoAF(delCO2_mean)) %>%
    select(-starts_with(c('delCO2', 'alk_meqL'))) %>%
    mutate(batch = case_when(batch == 1 ~ 'T0',
                             batch == 2 ~ 'T1'))%>%
    rename(DIC_mgL = DIC_mgL_mean, DIC_13C_ugL = DIC_13C_ugL_mean) %>%
    pivot_wider(names_from = 'batch',
                values_from = c('sample_datetime', 'AF', 'DIC_mgL', 'DIC_13C_ugL',
                                'DIC_mgL_sd', 'DIC_13C_ugL_sd')) %>%
    mutate(AF_exc = AF_T1 - AF_T0,
           DIC_pool = DIC_mgL_T1,
           DIC_mgL_T1 = DIC_mgL_T1 - DIC_mgL_T0,
           DIC_mgL_T0 = DIC_mgL_T0 - DIC_mgL_T0,
           DIC_13C_ugL_T1 = DIC_13C_ugL_T1 - DIC_13C_ugL_T0,
           DIC_13C_ugL_T0 = DIC_13C_ugL_T0 - DIC_13C_ugL_T0) %>%
    pivot_longer(cols = ends_with(c('T0', 'T1')),
                 names_to = c('.value', 'time'),
                 names_sep = '_T')

#code to save figure(change expt name)
dat_dif %>% filter(DIC_mgL != 0) %>%
ggplot( aes(trt, AF_exc * DIC_pool, col = factor(leachate), group = leachate))+
    geom_point()
geom_bar(stat = 'identity')
p12 <- ggplot(dat_dif, aes(sample_datetime, DIC_mgL, col = trt, lty = factor(leachate)))+
    geom_point() + geom_line(size = 1) +
    geom_errorbar(aes(ymin = DIC_mgL - DIC_mgL_sd,
                      ymax = DIC_mgL + DIC_mgL_sd), lty = 1)+
    xlab('Date') +theme_bw()
p13 <- ggplot(dat_dif, aes(sample_datetime, DIC_13C_ugL, col = trt, lty = factor(leachate)))+
    geom_point() + geom_line(size = 1) +
    geom_errorbar(aes(ymin = DIC_13C_ugL - DIC_13C_ugL_sd,
                      ymax = DIC_13C_ugL + DIC_13C_ugL_sd), lty = 1)+
    xlab('Date') +theme_bw()
png(filename = 'figures/ybc_expt_13CO2.png',
    width = 8, height = 5, res = 300, units = 'in')
ggpubr::ggarrange(p12, p13, common.legend = T)
dev.off()


#tables ####

dd <- dat_sum %>% select(-starts_with(c('delCO2', 'alk_meqL'))) %>%
    mutate(batch = case_when(batch == 1 ~ 'T0',
                             batch == 2 ~ 'T'))%>%
    rename(DIC_mgL = DIC_mgL_mean, DIC_13C_ugL = DIC_13C_ugL_mean) %>%
    pivot_wider(names_from = 'batch',
                values_from = c('sample_datetime', 'DIC_mgL', 'DIC_13C_ugL',
                                'DIC_mgL_sd', 'DIC_13C_ugL_sd')) %>%
    mutate(inc_time = sample_datetime_T - sample_datetime_T0,
           Delta_T = as.numeric(inc_time),
           Delta_DIC_ugLd = (DIC_mgL_T - DIC_mgL_T0)*1000/Delta_T,
           Delta_DIC_ugLd_sd = (DIC_mgL_sd_T + DIC_mgL_sd_T0)*1000/Delta_T,
           Delta_DIC_13_ugLd = (DIC_13C_ugL_T - DIC_13C_ugL_T0)/Delta_T,
           Delta_DIC_13_ugLd_sd = (DIC_13C_ugL_sd_T + DIC_13C_ugL_sd_T0)/Delta_T) %>%
    select(-ends_with(c('_T', '_T0')))

mutate(dd, AF = Delta_DIC_13_ugLd/(Delta_DIC_13_ugLd + Delta_DIC_ugLd),
       del = AFtodel(AF))

excess13CDIC <- dd %>% select( -inc_time) %>%
    pivot_wider(id_cols = trt, names_from = 'leachate',
                values_from = c('Delta_DIC_ugLd', 'Delta_DIC_13_ugLd',
                                'Delta_DIC_ugLd_sd', 'Delta_DIC_13_ugLd_sd')) %>%
    mutate(excess_13C_DIC = Delta_DIC_13_ugLd_1 - Delta_DIC_13_ugLd_0,
           excess_13C_sd = Delta_DIC_13_ugLd_sd_1 + Delta_DIC_13_ugLd_sd_0,
           excess_DIC = Delta_DIC_ugLd_1 - Delta_DIC_ugLd_0,
           excess_DIC_sd = Delta_DIC_ugLd_sd_1 + Delta_DIC_ugLd_sd_0) %>%
    select(-ends_with(c('_1','_0'))) %>%
    mutate(site = 'ybc', .before = trt)

excess13CDIC$carbon = c('Glucose', 'Glucose', 'No Glucose', 'No Glucose')
excess13CDIC$carbon <- factor(excess13CDIC$carbon, levels = c('No Glucose','Glucose'))
excess13CDIC$nuts = c('No Nutrients', 'Nutrients', 'No Nutrients', 'Nutrients')
png(filename = 'figures/ybc_expt_13CDIC_pres.png',
    width = 5, height = 5, res = 300, units = 'in')
ggplot(excess13CDIC, aes(carbon, excess_13C_DIC, col = nuts,fill = nuts, group = nuts)) +
    geom_line() + geom_point() +
    # geom_ribbon(aes(ymin = excess_13C_DIC - excess_13C_sd,
    #                   ymax = excess_13C_DIC + excess_13C_sd),col = 'transparent', alpha = 0.2,
    #             outline.type = 'full')+
    ggtitle('Yellow Bay Creek')+
    xlab('Glucose Treatment') +theme_bw()+
    ylab('Excess 13C DIC')+
    guides(color = guide_legend(title = "Nutrient Treatment"))+
    guides(fill = FALSE)
dev.off()

excess13CDIC$carbon = c('Glucose', 'Glucose', 'No Glucose', 'No Glucose')
excess13CDIC$carbon <- factor(excess13CDIC$carbon, levels = c('No Glucose','Glucose'))
excess13CDIC$nuts = c('No Nutrients', 'Nutrients', 'No Nutrients', 'Nutrients')
png(filename = 'figures/ybc_expt_13CDIC.png',
    width = 5, height = 5, res = 300, units = 'in')
ggplot(excess13CDIC, aes(carbon, excess_13C_DIC, col = nuts, group = nuts)) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = excess_13C_DIC - excess_13C_sd,
                      ymax = excess_13C_DIC + excess_13C_sd))+
    ggtitle('Yellow Bay Creek')
dev.off()

write_csv(excess13CDIC, 'data/excessDIC_ybc.csv')

png(filename = 'figures/ybc_expt_13CDIC_pres.png',
    width = 5, height = 5, res = 300, units = 'in')
ggplot(excess13CDIC, aes(carbon, excess_13C_DIC, col = nuts,fill = nuts, group = nuts)) +
    geom_line() + geom_point() +
    # geom_ribbon(aes(ymin = excess_13C_DIC - excess_13C_sd,
    #                   ymax = excess_13C_DIC + excess_13C_sd),col = 'transparent', alpha = 0.2,
    #             outline.type = 'full')+
    ggtitle('Yellow Bay Creek')+
    xlab('Glucose Treatment') +theme_bw()+
    ylab('Excess 13C DIC')+
    guides(color = guide_legend(title = "Nutrient Treatment"))+
    guides(fill = FALSE)
dev.off()

# write_csv(DIC_dataframe, "processed_CO2/experiment1_processed.csv")
# DIC <- read_csv('processed_CO2/experiment1_processed.csv')
#
# # calculate excess 13CDIC in each treatment ####
# sumDIC <- filter(DIC, time_interval %in% c('T0', 'T5')) %>%
#     select(time_interval, treatment, sample_datetime, DIC_molL, DIC_13C_umolL) %>%
#     group_by(time_interval, treatment) %>%
#     # summarize(across(.fns = list(mean = ~mean(.), sd = sd(.))))
#     summarize(sample_datetime = mean(sample_datetime),
#               DIC_molL = mean(DIC_molL),
#               DIC_13C_sd = sd(DIC_13C_umolL),
#               DIC_13C_umolL = mean(DIC_13C_umolL)) %>%
#     ungroup()
#
# DIC_table <- sumDIC %>%
#     mutate(excess_13C_DIC_umolL = DIC_13C_umolL - sumDIC$DIC_13C_umolL[sumDIC$time_interval == 'T0'],
#            excess_DIC_umolL = (DIC_molL - sumDIC$DIC_molL[sumDIC$time_interval == 'T0'])*10^6,
#            excess_13C_DIC_umolL_sd = DIC_13C_sd + sumDIC$DIC_13C_sd[sumDIC$time_interval == 'T0'])
#
#
# #creating excess accumulated 13C CO2 (from DIC) plot  ####
# # setwd('/Users/niko/Dropbox/leachate_linder/processed_CO2')
# # DIC_table <- read_csv('~/Dropbox/leachate_linder/processed_CO2/DIC_table.csv')
# DIC_table$nutrients <- c(0,0,1,0,1)
# DIC_table$sugar <- c(0,1,1,0,0)
# # DIC_table$DICXS <- DIC_table$DICXS*1000
# # DIC_table$STDev_T5 <- DIC_table$STDev_T5*1000
# # DIC_table$STDev_T0 <- DIC_table$STDev_T0*1000
# # DIC_table$STDev <- DIC_table$STDev_T0 + DIC_table$STDev_T5
#
# DIC_table <- slice(DIC_table, -1)
#
# plot(DIC_table$sugar[DIC_table$nutrients == 1],
#      DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1], type = 'b',
#      xlim = c(-0.25,1.25), ylim = c(0,1.5), pch = 15, col = 'blue',
#      xaxt = 'n', xlab = 'No Sugar                                         Sugar',
#      ylab = 'Accumulated Excess 13C CO2 (µg)')
#
#
# arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 1],
#        y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1],
#        x1 = DIC_table$sugar[DIC_table$nutrients == 1],
#        y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1]+
#          DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==1], length = 0.1, angle = 90)
# arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 1],
#       y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1],
#       x1 = DIC_table$sugar[DIC_table$nutrients == 1],
#       y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 1]-
#         DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==1], angle = 90, length = 0.1)
#
# lines(DIC_table$sugar[DIC_table$nutrients == 0],
#       DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0],pch = 1, col = 'forestgreen')
# points(DIC_table$sugar[DIC_table$nutrients == 0],
#        DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0], pch = 15, col = 'forestgreen')
# arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 0],
#        y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0],
#        x1 = DIC_table$sugar[DIC_table$nutrients == 0],
#        y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0]+
#          DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==0], length = 0.1, angle = 90)
# arrows(x0 = DIC_table$sugar[DIC_table$nutrients == 0],
#        y0 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0],
#        x1 = DIC_table$sugar[DIC_table$nutrients == 0],
#        y1 = DIC_table$excess_13C_DIC_umolL[DIC_table$nutrients == 0]-
#          DIC_table$excess_13C_DIC_umolL_sd[DIC_table$nutrients ==0], angle = 90, length = 0.1)
#
# legend('topleft', legend = c('Nutrients', 'No Nutrients'),
#        col = c('green', 'blue'), lty = c(1,1), pch = c(15, 15), border = NA)
#
# # plot(DIC_table$sample_datetime, dat$DIC_13C, type = 'b', pch = 15)
# lines(dat$sample_datetime)
#
# arrows()

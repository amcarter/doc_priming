
library(tidyverse)
library(lme4)

fl <- read_csv('data/DIC_change_by_bottle_fl.csv')
fl2 <- read_csv('data/DIC_change_by_bottle_fl2.csv')
ybc <- read_csv('data/DIC_change_by_bottle_ybc.csv')
nsc <- read_csv('data/DIC_change_by_bottle_nsc.csv')
nmc <- read_csv('data/DIC_change_by_bottle_nmc.csv')
tntp <- read_csv('data/nutrients/TNTP_2024.csv')

tntp <- tntp %>%
    filter(!is.na(Sample_ID)) %>%
    mutate(`P_ug/L` = case_when(`P_ug/L`== '< 1.5' ~ 0.75,
              TRUE ~ as.numeric(`P_ug/L`))) %>%
    select(-5, -6)

dat <- bind_rows(fl2, ybc, fl, nsc, nmc) %>%
    mutate(carbon = case_when(trt %in% c('C', 'CNP') ~ 1,
                              TRUE ~ 0),
           nutrients = case_when(trt %in% c('NP', 'CNP') ~ 1,
                                 # leachate == 1 ~ 1,
                                 TRUE ~ 0))

# calculate the increase in 13 C in each treatment that can be attributed to
# the breakdown of glucose or of lake carbon:
# for now, we will assume that the carbon in the nyack is similar in del value
# to that in the lake because I don't believe the measured values.

c_source <- read_csv('data/DOC_stock_concentrations.csv')
lake_AF <- c_source$AF[c_source$Sample_ID == 'lake water']
glu_AF <- c_source$AF[c_source$Sample_ID == 'glucose']

lch_N_ugL <- as.numeric(tntp[9, 3])
nutoct_N_ugL <- 4 * (1/1e6) * 14.0067 * 1e6
nutjul_N_ugl <- 400 * (1/1e6) * 14.0067 * 1e6
fl_N_ugl <- as.numeric(tntp[1, 3])
nmc_N_ugl <- as.numeric(tntp[2, 3])
nsc_N_ugl <- as.numeric(tntp[3, 3])
fl2_N_ugl <- as.numeric(tntp[5, 3])
ybc_N_ugl <- as.numeric(tntp[6, 3])

lch_P_ugL <- as.numeric(tntp[9, 4])
nutoct_P_ugL <- 0.25 * (1/1e6) * 30.9738 * 1e6
nutjul_P_ugl <- 50 * (1/1e6) * 30.9738 * 1e6
fl_P_ugl <- as.numeric(tntp[1, 4])
nmc_P_ugl <- as.numeric(tntp[2, 4])
nsc_P_ugl <- as.numeric(tntp[3, 4])
fl2_P_ugl <- as.numeric(tntp[5, 4])
ybc_P_ugl <- as.numeric(tntp[6, 4])

dat <- dat %>%
    mutate(site = case_when(site == 'fl' ~ 'flathead lake jul',
                            site == 'fl2' ~ 'flathead lake oct',
                            TRUE ~ site)) %>%
    mutate(lake_13C = case_when(DIC_ugLd < 0 ~ 0,
                                TRUE ~ DIC_ugLd * lake_AF/(1 - lake_AF)),
           glucose_13C = case_when(DIC_ugLd < 0 ~ 0,
                                   trt %in% c('C', 'CNP') ~ DIC_ugLd * glu_AF/(1 - glu_AF),
                                   TRUE ~ NA)) %>%
    rowwise() %>%
    mutate(background_13DIC_ugLd = mean(c_across(c('glucose_13C', 'lake_13C')),
                                        na.rm = T)) %>%
    mutate(range_bgrd_13DIC_ugLd = abs((glucose_13C - lake_13C)/2),
           excess_13C_DIC = case_when(DIC_13C_ugLd < 0 ~ 0,
                                      TRUE ~ DIC_13C_ugLd - background_13DIC_ugLd),
           DIC_12_ugLd = DIC_ugLd - DIC_13C_ugLd) %>%
    select(-glucose_13C, -lake_13C)

dat <- dat %>%
    mutate(experiment = ifelse(grepl("flathead lake oct|ybc", site), "fall", "summer")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), 1, 0)) %>%
    mutate(N_lch_ugL = case_when(leachate == 0 ~ 0,
                             TRUE ~ ((5 * lch_N_ugL)/80)),
           N_nuts_july_ugL = case_when(nutrients == 1 & experiment == "summer" ~ nutjul_N_ugl,
                                 TRUE ~ 0),
           N_nuts_october_ugL = case_when(nutrients == 1 & experiment == "fall" ~ nutoct_N_ugL,
                                          TRUE ~ 0),
           N_fl_ugL = case_when(site == "flathead lake oct" ~ fl_N_ugl,
                                TRUE ~ 0),
           N_nmc_ugL = case_when(site == "nyack_mc" ~ nmc_N_ugl,
                                 TRUE ~ 0),
           N_nsc_ugL = case_when(site == "nyack_sc" ~ nsc_N_ugl,
                                 TRUE ~ 0),
           N_fl2_ugL = case_when(site == "flathead lake jul" ~ fl2_N_ugl,
                                 TRUE ~ 0),
           N_ybc_ugL = case_when(site == "ybc" ~ ybc_N_ugl,
                                 TRUE ~ 0)) %>%
    mutate(totalN_ugL = N_lch_ugL + N_nuts_july_ugL + N_nuts_october_ugL + N_fl_ugL +
               N_nmc_ugL + N_nsc_ugL + N_fl2_ugL + N_ybc_ugL) %>%
    select(-starts_with("N_")) %>%
    mutate(P_lch_ugL = case_when(leachate == 0 ~ 0,
                                 TRUE ~ ((5 * lch_N_ugL)/80)),
           P_nuts_july_ugL = case_when(nutrients == 1 & experiment == "summer" ~ nutjul_P_ugl,
                                       TRUE ~ 0),
           P_nuts_october_ugL = case_when(nutrients == 1 & experiment == "fall" ~ nutoct_P_ugL,
                                          TRUE ~ 0),
           P_fl_ugL = case_when(site == "flathead lake oct" ~ fl_P_ugl,
                                TRUE ~ 0),
           P_nmc_ugL = case_when(site == "nyack_mc" ~ nmc_P_ugl,
                                 TRUE ~ 0),
           P_nsc_ugL = case_when(site == "nyack_sc" ~ nsc_P_ugl,
                                 TRUE ~ 0),
           P_fl2_ugL = case_when(site == "flathead lake jul" ~ fl2_P_ugl,
                                 TRUE ~ 0),
           P_ybc_ugL = case_when(site == "ybc" ~ ybc_P_ugl,
                                 TRUE ~ 0)) %>%
    mutate(totalP_ugL = P_lch_ugL + P_nuts_july_ugL + P_nuts_october_ugL + P_fl_ugL +
               P_nmc_ugL + P_nsc_ugL + P_fl2_ugL + P_ybc_ugL) %>%
    select(-starts_with("P_")) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose"))

## Nitrogen Plots

png(filename = 'figures/Nplot_site.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    filter(leachate == 1) %>%
    # filter(totalN_ugL < 1000) %>%
    ggplot(aes(totalN_ugL, excess_13C_DIC, col = carbon, shape = site)) +
    geom_point() +
    geom_smooth(data = . %>% filter(carbon == "Glucose"), method = "lm", se = FALSE) + # Trendline for glucose
    geom_smooth(data = . %>% filter(carbon == "No Glucose"), method = "lm", se = FALSE) + # Trendline without glucose
    facet_wrap(.~site, scales = "free") +
    ylab('Excess 13C DIC (ugL/d)') + xlab('Total N(ug/L)') +
    theme_bw()

dev.off()

png(filename = 'figures/Nplot_continuous.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    # filter(leachate == 1) %>%
    # filter(totalN_ugL < 1000) %>%
    ggplot(aes(totalN_ugL, excess_13C_DIC, col = carbon, shape = site)) +
    geom_point() +
    # geom_smooth(data = . %>% filter(carbon == "Glucose"), method = "lm", se = FALSE) + # Trendline for glucose
    # geom_smooth(data = . %>% filter(carbon == "No Glucose"), method = "lm", se = FALSE) + # Trendline without glucose
    ylab('Excess 13C DIC (ugL/d)') + xlab('Total N (ug/L)') +
    scale_x_continuous(trans = "log10") +
    theme_bw()

dev.off()

png(filename = 'figures/Nplot_12DIC.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    # filter(leachate == 1) %>%
    # filter(totalN_ugL < 1000) %>%
    ggplot(aes(totalN_ugL, DIC_12_ugLd, col = carbon, shape = site)) +
    geom_point() +
    # geom_smooth(data = . %>% filter(carbon == "Glucose"), method = "lm", se = FALSE) + # Trendline for glucose
    # geom_smooth(data = . %>% filter(carbon == "No Glucose"), method = "lm", se = FALSE) + # Trendline without glucose
    ylab('Excess 13C DIC (ugL/d)') + xlab('Total N (ug/L)') +
    scale_x_continuous(trans = "log10") +
    theme_bw()

dev.off()

## Phosphorous Plots

png(filename = 'figures/Pplot_site.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    filter(leachate == 1) %>%
    ggplot(aes(totalP_ugL, excess_13C_DIC, col = carbon, shape = site)) +
    geom_point() +
    geom_smooth(data = . %>% filter(carbon == "Glucose"), method = "lm", se = FALSE) + # Trendline for glucose
    geom_smooth(data = . %>% filter(carbon == "No Glucose"), method = "lm", se = FALSE) + # Trendline without glucose
    facet_wrap(.~site, scales = "free") +
    ylab('Excess 13C DIC (ugL/d)') + xlab('Total P (ug/L)') +
    theme_bw()

dev.off()

png(filename = 'figures/Pplot_continuous.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    # filter(leachate == 1) %>%
    ggplot(aes(totalP_ugL, excess_13C_DIC, col = carbon, shape = site)) +
    geom_point() +
    # geom_smooth(data = . %>% filter(carbon == "Glucose"), method = "lm", se = FALSE) + # Trendline for glucose
    # geom_smooth(data = . %>% filter(carbon == "No Glucose"), method = "lm", se = FALSE) + # Trendline without glucose
    ylab('Excess 13C DIC (ugL/d)') + xlab('Total P (ug/L)') +
    scale_x_continuous(trans = "log10", ) +
    coord_cartesian(xlim = c(-5, NA)) +
    theme_bw()

dev.off()

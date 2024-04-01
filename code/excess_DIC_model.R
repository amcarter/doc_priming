
library(tidyverse)
library(lme4)

fl <- read_csv('data/DIC_change_by_bottle_fl.csv')
nsc <- read_csv('data/DIC_change_by_bottle_nsc.csv')
nmc <- read_csv('data/DIC_change_by_bottle_nmc.csv')
fl2 <- read_csv('data/DIC_change_by_bottle_fl2.csv')
ybc <- read_csv('data/DIC_change_by_bottle_ybc.csv')

dat <- bind_rows(fl, nsc, nmc, fl2, ybc) %>%
    mutate(carbon = case_when(trt %in% c('C', 'CNP') ~ 1,
                              TRUE ~ 0),
           nutrients = case_when(trt %in% c('NP', 'CNP') ~ 1,
                                 leachate == 1 ~ 1,
                                 TRUE ~ 0))

# calculate the increase in 13 C in each treatment that can be attributed to
# the breakdown of glucose or of lake carbon:
# for now, we will assume that the carbon in the nyack is similar in del value
# to that in the lake because I don't believe the measured values.

c_source <- read_csv('data/DOC_stock_concentrations.csv')
lake_AF <- c_source$AF[c_source$Sample_ID == 'lake water']
glu_AF <- c_source$AF[c_source$Sample_ID == 'glucose']

dat <- dat %>%
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
                                      TRUE ~ DIC_13C_ugLd - background_13DIC_ugLd)) %>%
    select(-glucose_13C, -lake_13C)

# boxplot

png(filename = 'figures/bybottle/excess13CDIC_boxplot_bybottle.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                     TRUE ~ site)) %>%
    ggplot(aes(trt, excess_13C_DIC, fill = factor(leachate))) +
    geom_boxplot() +
    facet_wrap(~ site, nrow = 2, scales = "free_x", strip.position = "bottom") +
    ylab('excess 13C DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# strip plot

png(filename = 'figures/bybottle/excess13CDIC_stripplot_bybottle.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    ggplot(aes(trt, excess_13C_DIC, color = factor(leachate))) +
    geom_point(alpha = 0.5) +
    facet_wrap(~ site, nrow = 2, scales = "free_x", strip.position = "bottom") +
    ylab('excess 13C DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# excess plots
    # glucose presence

dat_mean <- dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    filter(leachate == 1) %>%
    group_by(site, trt, leachate) %>%
    summarize(mean_excess13CDIC = mean(excess_13C_DIC)) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ungroup

png(filename = 'figures/bybottle/excess13CDIC_bybottle.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    filter(leachate == 1) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ggplot(aes(carbon, excess_13C_DIC, color = factor(nuts))) +
    geom_point(alpha = 0.5) +
    geom_line(data = dat_mean, aes(y = mean_excess13CDIC, group = nuts)) +
    facet_wrap(~ site, nrow = 2, scales = "free_x", strip.position = "bottom") +
    ylab('excess 13C DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c("No Glucose", "Glucose"))

dev.off()

    # leachate presence

dat_mean_lch <- dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    filter(carbon == 0) %>%
    group_by(site, leachate, trt) %>%
    summarize(mean_DIC = mean(DIC_ugLd)) %>%
    mutate(leachate = ifelse(grepl(0, leachate), "No Leachate", "Leachate")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ungroup

png(filename = 'figures/bybottle/excessDIC_bybottle.png',
    width = 10, height = 5, res = 300, units = 'in')

dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    filter(carbon == 0) %>%
    mutate(leachate = ifelse(grepl(0, leachate), "No Leachate", "Leachate")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ggplot(aes(x = leachate, y = DIC_ugLd, color = factor(nuts))) +
    geom_point(alpha = 0.5) +
    geom_line(data = dat_mean_lch, aes(y = mean_DIC, group = nuts)) +
    facet_wrap(~ site, nrow = 2, scales = "free_x", strip.position = "bottom") +
    ylab('total DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c("No Leachate", "Leachate"))

dev.off()

    # presence

dat_mean_3 <- dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    # filter(carbon == 0) %>%
    group_by(site, leachate, trt) %>%
    summarize(mean_DIC = mean(DIC_ugLd)) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(leachate = ifelse(grepl(0, leachate), "No Leachate", "Leachate")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ungroup

# png(filename = 'figures/bybottle/excessDIC_bybottle.png',
#     width = 10, height = 5, res = 300, units = 'in')

dat %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site)) %>%
    # filter(carbon == 0) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(leachate = ifelse(grepl(0, leachate), "No Leachate", "Leachate")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ggplot(aes(x = carbon, y = DIC_ugLd, color = factor(nuts))) +
    geom_point(alpha = 0.5) +
    geom_line(data = dat_mean_3, aes(y = mean_DIC, group = nuts)) +
    facet_grid(leachate ~ site, scales = "free_x") +
    ylab('total DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c("No Glucose", "Glucose"))

dev.off()
    ##
    ##

combined <- dat %>%
    group_by(site, trt, leachate) %>%
    summarize(mean_excess13CDIC = mean(excess_13C_DIC)) %>%
    ungroup

excess13CDIC_bybottle <- combined %>%
    pivot_wider(id_cols = c(trt, site,),
                names_from = 'leachate',
                values_from = 'mean_excess13CDIC') %>%
    mutate(excess_13C_DIC = `1` - `0`) %>%
    select(-(c(`1`,`0`))) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients"))

excess13CDIC_bybottle$carbon <- factor(excess13CDIC_bybottle$carbon, levels = c("No Glucose", "Glucose"))


ggplot(excess13CDIC_bybottle, aes(carbon, excess_13C_DIC, col = nuts, group = nuts)) +
    geom_line() +
    # geom_point(data = combined) +
    # geom_errorbar(aes(ymin = excess_13C_DIC - excess_13C_sd,
    #                   ymax = excess_13C_DIC + excess_13C_sd)) +
    geom_point() +
    # geom_point(data = combined, aes(carbon, mean_excess13CDIC, col = nuts), alpha = 0.5) + #starting trying to add points in, not sure if this is the correct data
    facet_grid(cols = vars(site)) +  # Facet by site
    ggtitle('Excess 13C DIC')  # Title for the entire plot

dev.off()

# model trials below here
dat <-
    bind_rows(fl, nsc, nmc, fl2, ybc) %>%
    mutate(carbon = case_when(trt %in% c('C', 'CNP') ~ 1,
                              TRUE ~ 0),
           nutrients = case_when(trt %in% c('NP', 'CNP') ~ 1,
                                 TRUE ~ 0))

ggplot(dat, aes(trt, DIC_ugLd, fill = factor(leachate))) +
    geom_boxplot()+
    facet_wrap(.~site)
ggplot(dat, aes(trt, DIC_13C_ugLd, fill = factor(leachate))) +
    geom_boxplot()+
    facet_wrap(.~site)


mod_fl <- lm(DIC_13C_ugLd ~  carbon + nutrients + leachate,
          data = fl)
mod <- lmer(DIC_13C_ugLd ~ (1|site)  + carbon + nutrients + leachate,
          data = dat)
summary(mod_fl)
ranef(mod)

#
# mod <- lmer(DIC_13C_ugLd ~ (1 + carbon + nutrients|site) +
#                 carbon + nutrients + leachate,
#           data = dat)
#
# summary(mod)
# ranef(mod)

mod <- lm(DIC_13C_ugLd ~ carbon + nutrients + leachate + carbon*nutrients,
          data = fl)
mod <- lmer(DIC_13C_ugLd ~ (1|site) + carbon + nutrients + leachate +
              carbon*leachate + nutrients*leachate + carbon*nutrients,
          data = dat)

summary(mod)
ranef(mod)

# Data without nutrients
no_nuts <- dat %>% filter(nutrients == 0)
mod <- lmer(DIC_13C_ugLd ~ (1|site) + carbon +leachate +
              carbon*leachate,
          data = no_nuts)

summary(mod)
ranef(mod)




library(tidyverse)
library(lme4)

fl <- read_csv('data/DIC_change_by_bottle_fl.csv')
fl2 <- read_csv('data/DIC_change_by_bottle_fl2.csv')
ybc <- read_csv('data/DIC_change_by_bottle_ybc.csv')
nsc <- read_csv('data/DIC_change_by_bottle_nsc.csv')
nmc <- read_csv('data/DIC_change_by_bottle_nmc.csv')

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

dat <- dat %>%
    mutate(site = case_when(site == 'fl' ~ 'flathead lake jul',
                            site == 'fl2' ~ 'flathead lake oct',
                     TRUE ~ site)) %>%
    # filter(site== 'flathead lake oct') %>%
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


dat_mean <- dat %>%
    # filter(leachate == 1) %>%
    group_by(site, trt, leachate) %>%
    mutate(nutrients = factor(mean(nutrients)),
           carbon = factor(mean(carbon), levels = c('0', '1')),
           excess_13C_DIC = mean(excess_13C_DIC),
           DIC_12_ugLd = mean(DIC_12_ugLd))

dat %>%
    ggplot(aes(trt, DIC_12_ugLd, fill = factor(leachate))) +
    geom_boxplot() +
    facet_wrap(.~site) +
    ylab('12C DIC (ugL/d)') + xlab("") +
    theme_bw()
dat %>%
    ggplot(aes(trt, excess_13C_DIC, fill = factor(leachate))) +
    geom_boxplot() +
    facet_wrap(.~site) +
    ylab('excess 13C DIC (ugL/d)') + xlab("") +
    theme_bw()

dat %>%
    # filter(leachate == 0) %>%
    mutate(nutrients = factor(nutrients),
           carbon = factor(carbon, levels = c('0', '1'))) %>%
    ggplot(aes(carbon, DIC_12_ugLd, col = nutrients, group = nutrients)) +
    geom_point() +
    geom_line(data = dat_mean) +
    facet_grid(leachate~site) +
    ylab('12C Breakdown (ugL/d)') + xlab("") +
    theme_bw()



# model trials below here
dat <-
    bind_rows(fl, nsc, nmc) %>%
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


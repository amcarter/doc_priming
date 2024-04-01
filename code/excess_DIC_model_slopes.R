
library(tidyverse)
library(lme4)

fl <- read_csv('data/DIC_change_by_bottle_fl.csv')
nsc <- read_csv('data/DIC_change_by_bottle_nsc.csv')
nmc <- read_csv('data/DIC_change_by_bottle_nmc.csv')
fl2 <- read_csv('data/DIC_change_by_bottle_fl2.csv')
ybc <- read_csv('data/DIC_change_by_bottle_ybc.csv')

c_source <- read_csv('data/DOC_stock_concentrations.csv')
lake_AF <- c_source$AF[c_source$Sample_ID == 'lake water']
glu_AF <- c_source$AF[c_source$Sample_ID == 'glucose']

dat <- bind_rows(fl, nsc, nmc, fl2, ybc) %>%
    mutate(carbon = case_when(trt %in% c('C', 'CNP') ~ 1,
                              TRUE ~ 0),
           nutrients = case_when(trt %in% c('NP', 'CNP') ~ 1,
                                 # leachate == 1 ~ 1,
                                 TRUE ~ 0))

# calculate the increase in 13 C in each treatment that can be attributed to
# the breakdown of glucose or of lake carbon:
# for now, we will assume that the carbon in the nyack is similar in del value
# to that in the lake because I don't believe the measured values.

dat <- dat %>%
    mutate(lake_13C = case_when(DIC_ugLd < 0 ~ 0,
                                TRUE ~ DIC_ugLd * lake_AF/(1 - lake_AF)),
           glucose_13C = case_when(DIC_ugLd < 0 ~ 0,
                                   trt %in% c('C', 'CNP') ~
                                       DIC_ugLd * glu_AF/(1 - glu_AF),
                                   TRUE ~ NA)) %>%
    rowwise() %>%
    mutate(background_13DIC_ugLd = mean(c_across(c('glucose_13C', 'lake_13C')),
                                        na.rm = T)) %>%
    mutate(range_bgrd_13DIC_ugLd = abs((glucose_13C - lake_13C)/2),
           excess_13C_DIC = case_when(DIC_13C_ugLd < 0 ~ 0,
                                      TRUE ~ DIC_13C_ugLd - background_13DIC_ugLd)) %>%
    select(-glucose_13C, -lake_13C) %>%
    mutate(site = case_when(site == 'fl' ~ 'Flathead Lake, July',
                            site == 'fl2'~ 'Flathead Lake, October',
                            site == 'ybc' ~ 'Yellow Bay Creek',
                            site == 'nyack_mc' ~ 'Nyack Main Channel',
                            site == 'nyack_sc' ~ 'Nyack Side Channel',
                            TRUE ~ site))

# excess plots
    # glucose presence

dat_mean <- dat %>%
    # filter(leachate == 1) %>%
    group_by(site, trt, leachate) %>%
    summarize(mean_excess13CDIC = mean(excess_13C_DIC),
              mean_DIC = mean(DIC_ugLd)) %>%
    mutate(leachate = ifelse(grepl(0, leachate), "No Leachate", "Leachate"),
           carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose"),
           nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ungroup


dat %>%
    filter(leachate == 1) %>%
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ggplot(aes(carbon, excess_13C_DIC, color = factor(nuts))) +
    geom_point(alpha = 0.5) +
    geom_line(data = filter(dat_mean, leachate == "Leachate"),
              aes(y = mean_excess13CDIC, group = nuts)) +
    facet_wrap(~ site, nrow = 2, scales = "free_x", strip.position = "bottom") +
    ylab('excess 13C DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c("No Glucose", "Glucose"))


# calculate slopes for excess 13C:
# for now, we will filter out the nutrient treatments because we cannot differentiate
# a nutrient addition from a leachate addition. We will address the nutrient effect
# in the across sites model with continuous nutrients as a covariate

dat_leachate <- filter(dat, leachate == 1, nutrients == 0)
dat_leachate <- filter(dat, leachate == 1)

mod_excess13C <- lmer(excess_13C_DIC ~ (carbon | site),
                      data = dat_leachate)

mod_excess13C <- lmer(excess_13C_DIC ~ (carbon*nutrients | site),
                      data = dat_leachate)

ranef(mod_excess13C)
dotplot.ranef.mer(ranef(mod_excess13C, whichel = c('carbon')), level = 0.95)
dd <- as.data.frame(ranef(mod_excess13C))

dd %>% filter(term != '(Intercept)') %>%
ggplot(aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0) +
    geom_vline(xintercept = 0)+
    ylab('Site')+
    xlab('Treatment effect') +
    theme_bw()

dd %>% filter(term != '(Intercept)') %>%
    rename(site = grp,
           effect = term,
           model_estimate = condval,
           est_sd = condsd) %>%
    mutate(CI_95 = est_sd * 1.96) %>%
    select(-grpvar) %>%
    write_csv('data/model_fits/carbon_effect_on_excess13C.csv')

extract_fit <- function(fit, site_name){
    data <- filter(dat_leachate, site == site_name)
    fl_j <- data.frame(summary(fit)$coefficients)
    colnames(fl_j) <- c('model_estimate', 'est_se', 't_value', 'p_value')
    fl_j <- mutate(fl_j, effect = rownames(fl_j),
                   est_sd = est_se * sqrt(nrow(data)),
                   CI_95 = 1.96 * est_se,
                   site = site_name)
    fl_j <- filter(fl_j, effect != '(Intercept)')
    rownames(fl_j) <- NULL

    return(fl_j)

}
# compare the results to models on individual sites:
mod_fl_j <- lm(excess_13C_DIC ~ carbon * nutrients,
               data = filter(dat_leachate, site == 'Flathead Lake, July')) #significant
fl_j <- extract_fit(mod_fl_j, 'Flathead Lake, July')

mod_fl_o <- lm(excess_13C_DIC ~ carbon * nutrients,
               data = filter(dat_leachate, site == 'Flathead Lake, October'))
summary(mod_fl_o)
site_fits <- extract_fit(mod_fl_o, 'Flathead Lake, October') %>%
    bind_rows(fl_j)
mod_nmc <- lm(excess_13C_DIC ~ carbon * nutrients,
               data = filter(dat_leachate, site == 'Nyack Main Channel'))
summary(mod_nmc) #significant
site_fits <- extract_fit(mod_nmc, 'Nyack Main Channel') %>%
    bind_rows(site_fits)

mod_nsc <- lm(excess_13C_DIC ~ carbon * nutrients,
               data = filter(dat_leachate, site == 'Nyack Side Channel'))
summary(mod_nsc) #maybe
site_fits <- extract_fit(mod_nsc, 'Nyack Side Channel') %>%
    bind_rows(site_fits)

mod_ybc <- lm(excess_13C_DIC ~ carbon * nutrients,
               data = filter(dat_leachate, site == 'Yellow Bay Creek'))
summary(mod_ybc)
site_fits <- extract_fit(mod_ybc, 'Yellow Bay Creek') %>%
    bind_rows(site_fits)


png(filename = 'figures/treatment_effects_all_sites.png',
    width = 8, height = 4, res = 300, units = 'in')
site_fits %>%
    mutate(sig = case_when(p_value <= 0.05 ~ 1,
                           TRUE ~ 0)) %>%
    # ggplot(aes(y = site, x = model_estimate, col = factor(sig))) +
    ggplot(aes(y = site, x = model_estimate)) +
    geom_point() +
    facet_wrap( ~ effect, scales = "free_x") +
    geom_errorbarh(aes(xmin = model_estimate - CI_95,
                       xmax = model_estimate + CI_95), height=0) +
    geom_vline(xintercept = 0)+
    ylab('Site')+
    xlab('Treatment effect') +
    theme_bw()

dev.off()

# model the rate of change in DIC across treatments
dat %>%
    mutate(leachate = ifelse(grepl(0, leachate), "No Leachate", "Leachate"),
           carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose"),
           nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    ggplot(aes(x = carbon, y = DIC_ugLd, color = factor(nuts))) +
    geom_point(alpha = 0.5) +
    geom_line(data = dat_mean, aes(y = mean_DIC, group = nuts)) +
    facet_grid(leachate ~ site, scales = "free_x") +
    ylab('total DIC (ugL/d)') +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c("No Glucose", "Glucose"))

# calculate slopes for accumulated DIC:

dat2 <- mutate(dat, nutrients = case_when(leachate == 1 ~ 1,
                                     TRUE ~ nutrients))
mod_DIC <- lmer(DIC_ugLd ~ (carbon * nutrients| site),
                      data = dat2)

ranef(mod_DIC)
dotplot.ranef.mer(ranef(mod_DIC))
dd <- as.data.frame(ranef(mod_DIC))

# plot the modeled effects with 95% confidence intervals:
dd %>% filter(term != '(Intercept)') %>%
ggplot(aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0) +
    geom_vline(xintercept = 0)+
    ylab('Site')+
    xlab('Treatment effect') +
    theme_bw()

dd %>% filter(term != '(Intercept)') %>%
    rename(site = grp,
           effect = term,
           model_estimate = condval,
           est_sd = condsd) %>%
    mutate(CI_95 = est_sd * 1.96) %>%
    select(-grpvar) %>%
    write_csv('data/model_fits/treatment_effects_on_DIC.csv')



library(tidyverse)

doc <- read_csv('data/DOC_summer_2023.csv')%>%
    select(Sample_ID, sample_type, Collection_date, site, trt, leachate, DOC_mgL,
           Max_12CO2_ppm, CO2_12_Integral, CO2_13_Integral, Delta_13C)


doc %>%
ggplot( aes(DOC_mgL, CO2_12_Integral, col = sample_type)) + geom_point()

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

    del
}
# develop calibration curve:

doc_cal <- doc %>% filter(sample_type == 'cal')
doc_cal %>% filter(grepl('^TOC', Sample_ID)) %>%
    ggplot(aes( CO2_12_Integral, CO2_13_Integral, col = Delta_13C)) +
    geom_point() + facet_wrap(.~Sample_ID, scales = 'free') + theme_bw()


cal <- lm(DOC_mgL ~ CO2_12_Integral, doc_cal)
ggplot(doc_cal, aes(DOC_mgL, CO2_12_Integral, col = Delta_13C)) +
    geom_point() + geom_smooth(method = 'lm', se = FALSE)

coefs <- cal$coefficients

doc <- doc %>%
    mutate(DOC_12_mgL = CO2_12_Integral * coefs[2] + coefs[1],
           AF = deltoAF(Delta_13C),
           # DOC132 = CO2_13_Integral * coefs[2] + coefs[1],
           DOC_13_mgL = DOC_12_mgL * (AF/(1-AF))) %>%
    # select( -Max_12CO2_ppm, -ends_with(c('Integral', 'Baseline'))) %>%
    filter(!is.na(Sample_ID))

# ggplot(doc, aes(DOC132, DOC_13_mgL, col = trt, pch = factor(leachate))) +
    # geom_point() +
    # geom_abline( slope = 1, intercept = 0)

# Determine the concentration of leachate:

conc <- doc %>% filter(sample_type == 'conc')
gluc <- conc %>%
    select(-trt, -Collection_date, -sample_type) %>%
    filter(grepl('^glucose', Sample_ID)) %>%
    mutate(Sample_ID = substr(Sample_ID, 1, nchar(Sample_ID)-2))
conc <- conc %>%
    mutate(Sample_ID = sub("[ _][AB]$", "", Sample_ID)) %>%
    filter(Sample_ID != 'glucose_leachate') %>%
    select(Sample_ID, delta_13C = Delta_13C, AF, conc_stock_mgL = DOC_mgL) %>%
    group_by(Sample_ID) %>%
    summarize(across(everything(), mean)) %>%
    mutate(conc_sample_mgL = conc_stock_mgL)

conc$conc_stock_mgL[conc$Sample_ID == 'glucose'] <- 300
conc$conc_sample_mgL[conc$Sample_ID == 'glucose'] <- 1

AF_gluc <- mean(gluc$AF[gluc$Sample_ID == 'glucose'])

leachate_dilution <- 0.3/100

leach <- gluc %>%
    mutate(glucose_13C = DOC_12_mgL * (AF_gluc/(1-AF_gluc)),
           extra_13C = DOC_13_mgL - glucose_13C,
           leachate_13C = extra_13C/leachate_dilution) %>%
    filter(leachate ==1)

leachate_conc_mgL = mean(leach$leachate_13C)
leach <- data.frame(Sample_ID = 'leachate',
                    delta_13C = Inf,
                    AF = 1,
                    conc_stock_mgL = leachate_conc_mgL,
                    conc_sample_mgL = leachate_conc_mgL/100)
stocks <- rbind(conc, leach)
write_csv(stocks, 'data/DOC_stock_concentrations.csv')

# save DOC data for samples ####
dd <- doc %>% filter(!is.na(trt)) %>%
    select(-sample_type) %>%
    filter(site != 'nyack_well') %>%
    mutate(site = ifelse(site == "flathead", "fl", site)) %>%
    mutate(date = as.Date(Collection_date, format = '%m/%d/%Y'))
write_csv(dd, 'data/DOC_sample_data.csv')

filter(dd, date > as.Date('2023-07-11')) %>%
    group_by(date, trt, leachate, site) %>%
    summarize(across(everything(), mean)) %>%

ggplot( aes(date, DOC_13_mgL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 2) +
    facet_wrap(.~site, scales = 'free') + theme_classic()

#excess13CDOC$Collection_date <- as.Date(excess13CDOC$Collection_date, format = "%m/%d/%Y")%>%
#how to include as.Date function in pipe?

# DOC processing to match calculate excess CO2 code output ####

#### chunk of excess13CDOC code below = "dat_sum" in calculate_excess_CO2

excess13CDOC <- dd %>% select(-Sample_ID, -Collection_date, -Delta_13C, -AF) %>%
    #mutate(date = as.Date(Collection_date, format = "%m/%d/%Y"))
    mutate(carbon = ifelse(grepl("C", trt), "Glucose", "No Glucose")) %>%
    mutate(nuts = ifelse(grepl("NP", trt), "Nutrients", "No Nutrients")) %>%
    group_by(site, trt) %>%
    mutate (batch = case_when(date %in% as.Date(c("2023-07-13", "2023-07-21"))~1,
                              date %in% as.Date(c("2023-07-18", "2023-07-27"))~2,
                              TRUE~NA_real_)) %>%
        ungroup() %>%
    filter(!is.na(batch)) %>%
    group_by(site, trt, leachate, batch, date, carbon, nuts) %>%
    summarise(
        DOC_mgL_mean = mean(DOC_mgL, na.rm = TRUE),
        DOC_mgL_sd = sd(DOC_mgL, na.rm = TRUE),
        DOC_13_mgL_mean = mean(DOC_13_mgL, na.rm = TRUE),
        DOC_13_mgL_sd = sd(DOC_13_mgL, na.rm = TRUE),
        DOC_12_mgL_mean = mean(DOC_12_mgL, na.rm = TRUE)
        )

dat_dif <- excess13CDOC %>%
    select(-DOC_12_mgL_mean, -DOC_13_mgL_sd, -DOC_mgL_sd) %>%
    mutate(batch = case_when(batch == 1 ~ 'T0',
                             batch == 2 ~ 'T1')) %>%
    rename(DOC_mgL = DOC_mgL_mean, DOC_13C_mgL = DOC_13_mgL_mean) %>%
    pivot_wider(names_from = 'batch',
                values_from = c('date', 'DOC_mgL', 'DOC_13C_mgL')) %>%
    mutate( DOC_mgL_T1 = DOC_mgL_T1 - DOC_mgL_T0,
            DOC_mgL_T0 = DOC_mgL_T0 - DOC_mgL_T0,
            DOC_13C_ugL_T1 = (DOC_13C_mgL_T1 - DOC_13C_mgL_T0)*1000,
            DOC_13C_ugL_T0 = (DOC_13C_mgL_T0 - DOC_13C_mgL_T0)*1000) %>%
    pivot_longer(cols = ends_with(c('T0', 'T1')),
                 names_to = c('.value', 'sitetime'), ######what does this do?
                 names_sep = '_T')

#baseline removed plots ####

# For site "fl"
p12_fl <- ggplot(dat_dif[dat_dif$site == "fl", ], aes(date, DOC_mgL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 1) +
    xlab('Date') + theme_bw() +
    facet_wrap(~ site)

p13_fl <- ggplot(dat_dif[dat_dif$site == "fl", ], aes(date, DOC_13C_ugL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 1) +
    xlab('Date') + theme_bw() +
    facet_wrap(~ site)

# For site "nyack_mc"
p12_nmc <- ggplot(dat_dif[dat_dif$site == "nyack_mc", ], aes(date, DOC_mgL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 1) +
    xlab('Date') + theme_bw() +
    facet_wrap(~ site)

p13_nmc <- ggplot(dat_dif[dat_dif$site == "nyack_mc", ], aes(date, DOC_13C_ugL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 1) +
    xlab('Date') + theme_bw() +
    facet_wrap(~ site)

# For site "nyack_sc"
p12_nsc <- ggplot(dat_dif[dat_dif$site == "nyack_sc", ], aes(date, DOC_mgL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 1) +
    xlab('Date') + theme_bw() +
    facet_wrap(~ site)

p13_nsc <- ggplot(dat_dif[dat_dif$site == "nyack_sc", ], aes(date, DOC_13C_ugL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 1) +
    xlab('Date') + theme_bw() +
    facet_wrap(~ site)

# Saving plots
png(filename = 'figures/fl_DOC_nobase.png', width = 8, height = 5, res = 300, units = 'in')
print(ggpubr::ggarrange(p12_fl, p13_fl, common.legend = TRUE))
dev.off()

png(filename = 'figures/nmc_DOC_nobase.png', width = 8, height = 5, res = 300, units = 'in')
print(ggpubr::ggarrange(p12_nmc, p13_nmc, common.legend = TRUE))
dev.off()

png(filename = 'figures/nsc_DOC_nobase.png', width = 8, height = 5, res = 300, units = 'in')
print(ggpubr::ggarrange(p12_nsc, p13_nsc, common.legend = TRUE))
dev.off()

png(filename = 'figures/12_DOCDIC_nobase.png', width = 15, height = 5, res = 300, units = 'in')
ggpubr::ggarrange(
    cowplot::plot_grid(p12_fl_DOC, p12_nmc_DOC, p12_nsc_DOC, ncol = 3),
    cowplot::plot_grid(p12_fl, p12_nmc, p12_nsc, ncol = 3),
    nrow = 2,
    common.legend = TRUE,
    legend = "none"
)
dev.off()

png(filename = 'figures/13_DOCDIC_nobase.png', width = 15, height = 5, res = 300, units = 'in')
ggpubr::ggarrange(
    cowplot::plot_grid(p13_fl_DOC, p13_nmc_DOC, p13_nsc_DOC, ncol = 3),
    cowplot::plot_grid(p13_fl, p13_nmc, p13_nsc, ncol = 3),
    nrow = 2,
    common.legend = TRUE,
    legend = "none"
)
dev.off()



#####

# p12 <- ggplot(excess13CDOC, aes(batch, DOC_mgL_mean, col = trt, lty = factor(leachate)))+
#     geom_point() + geom_line(size = 1) +
#     geom_errorbar(aes(ymin = DOC_mgL_mean - DOC_mgL_sd,
#                       ymax = DOC_mgL_mean + DOC_mgL_sd), lty = 1)+
#     xlab('Date') +theme_bw()
# p13 <- ggplot(excess13CDOC, aes(batch, DOC_13_mgL_mean, col = trt, lty = factor(leachate)))+
#     geom_point() + geom_line(size = 1) +
#     geom_errorbar(aes(ymin = DOC_13_mgL_mean - DOC_13_mgL_sd,
#                       ymax = DOC_13_mgL_mean + DOC_13_mgL_sd), lty = 1)+
#     xlab('Date') +theme_bw()
#
# png(filename = 'figures/DOC_summer_13CO2overtime.png',
#     width = 8, height = 5, res = 300, units = 'in')
# ggpubr::ggarrange(p12, p13, common.legend = T)
# dev.off()

# Split the dataset by site
# Split the dataset by site
site_data <- split(excess13CDOC, excess13CDOC$site)

# Create a list to store the plots for each site
plots_list <- list()

# Iterate over each site
for (site in names(site_data)) {
    # Filter data for the current site
    filtered_data <- site_data[[site]]

    # Create plot for DOC_mgL_mean
    p12 <- ggplot(filtered_data, aes(date, DOC_mgL_mean, col = trt, lty = factor(leachate))) +
        geom_point() + geom_line(size = 1) +
        xlab('Date') + theme_bw() +
        ggtitle(paste('DOC mg/L Over Time - Site:', site))

    # Create plot for DOC_13_mgL_mean
    p13 <- ggplot(filtered_data, aes(date, DOC_13_mgL_mean, col = trt, lty = factor(leachate))) +
        geom_point() + geom_line(size = 1) +
        xlab('Date') + theme_bw() +
        ggtitle(paste('DOC 13 mg/L Over Time - Site:', site))

    # Add the plots to the list
    plots_list[[site]] <- list(p12, p13)
}

# Plot the list of plots for each site
for (site in names(plots_list)) {
    # Arrange the plots
    arranged_plot <- ggpubr::ggarrange(plotlist = plots_list[[site]], common.legend = TRUE)

    # Save the arranged plot
    filename <- paste0('figures/', site, '_DOC_over_time.png')
    ggsave(filename, arranged_plot, width = 8, height = 5, units = 'in', dpi = 300)
}

# Plot for the entire dataset
png(filename = 'figures/12CDOC_over_time_site.png',
    width = 7, height = 5, res = 300, units = 'in')
ggplot(excess13CDOC, aes(date, DOC_12_mgL_mean, col = site, shape = factor(leachate))) +
    geom_point() +
    scale_x_date(date_breaks = "2 day", date_labels = "%b %d") +
    xlab('Date') + theme_bw() +
    ggtitle('12C DOC Over Time')
dev.off()

#### this point is where I can make the figures with baseline removed.
#### chunk of excess13CDOC code below = "dd" in calculate_excess_CO2

excess13CDOC <-  excess13CDOC %>%
    mutate(batch = case_when(batch == 1 ~ 'T0',
                             batch == 2 ~ 'T')) %>%
    rename(DOC_mgL = DOC_mgL_mean, DOC_13_mgL = DOC_13_mgL_mean, DOC_12_mgL = DOC_12_mgL_mean) %>%
    pivot_wider(names_from = 'batch',
                values_from = c('date', 'DOC_mgL', 'DOC_13_mgL',
                                'DOC_mgL_sd', 'DOC_13_mgL_sd', 'DOC_12_mgL')) %>%
    mutate(inc_time = date_T - date_T0,
           Delta_T = as.numeric(inc_time),
           Delta_DOC_ugLd = (DOC_mgL_T - DOC_mgL_T0)*1000/Delta_T, #multiplying all by 1000 to put into ug/L
           Delta_DOC_ugLd_sd = (DOC_mgL_sd_T + DOC_mgL_sd_T0)*1000/Delta_T,
           Delta_DOC_13_ugLd = (DOC_13_mgL_T - DOC_13_mgL_T0)*1000/Delta_T,
           Delta_DOC_13_ugLd_sd = (DOC_13_mgL_sd_T + DOC_13_mgL_sd_T0)*1000/Delta_T,
           Delta_DOC_12_ugLd = (DOC_12_mgL_T - DOC_12_mgL_T0)*1000/Delta_T) %>%
    select(-ends_with(c('_T', '_T0')))

#### chunk of excess13CDOC code below = "excess13CDIC" in calculate_excess_CO2

excess13CDOC <- excess13CDOC %>%
    select( -inc_time) %>%
    pivot_wider(id_cols = c(site, trt, carbon, nuts), names_from = c('leachate'),
                values_from = c('Delta_DOC_ugLd', 'Delta_DOC_13_ugLd',
                                'Delta_DOC_ugLd_sd', 'Delta_DOC_13_ugLd_sd', 'Delta_DOC_12_ugLd')) %>%
    mutate(excess_13C_DOC = Delta_DOC_13_ugLd_1 - Delta_DOC_13_ugLd_0,
           excess_13C_sd = Delta_DOC_13_ugLd_sd_1 + Delta_DOC_13_ugLd_sd_0,
           excess_12C_DOC = Delta_DOC_12_ugLd_1 - Delta_DOC_12_ugLd_0) %>%
    select(-ends_with(c('_1','_0')))

write_csv(excess13CDOC, 'data/excessDOC_summer.csv')

for (i in 1:3) {
    # png(filename = 'figures/summerDOC_pres.png',
    #     width = 5, height = 5, res = 300, units = 'in')
    #
    plot_data <- subset(excess13CDOC, site == unique(excess13CDOC$site)[i])
    plot <- ggplot(plot_data,
           aes(carbon, excess_13C_DOC, col = nuts, fill = nuts, group = nuts)) +
    geom_line() + geom_point() +
    # geom_ribbon(aes(ymin = excess_13C_DIC - excess_13C_sd,
    #                   ymax = excess_13C_DIC + excess_13C_sd),col = 'transparent', alpha = 0.2,
    #             outline.type = 'full')+
    ggtitle(unique(excess13CDOC$site)[i])+
    xlab('Glucose Treatment') +theme_bw()+
    ylab('Excess 13C DIC')+
    guides(color = guide_legend(title = "Nutrient Treatment"))+
    guides(fill = FALSE) +
    facet_wrap(~site)

    print(plot)

    ggsave(filename = paste0("figures/", unique(excess13CDOC$site)[i], "_DOCplot.png"),
           plot = plot, width = 5, height = 5, units = "in", dpi = 300)

# dev.off()
}


# ggplot(excess13CDOC, aes(date, DOC_mgL, col = site)) +
# geom_point()



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
cal <- lm(DOC_mgL ~ CO2_12_Integral, doc_cal)
ggplot(doc_cal, aes(DOC_mgL, CO2_12_Integral, col = Delta_13C)) +
    geom_point() + geom_smooth(method = 'lm', se = FALSE)

coefs <- cal$coefficients

doc <- doc %>%
    mutate(DOC_12_mgL = CO2_12_Integral * coefs[2] + coefs[1],
           AF = deltoAF(Delta_13C),
           DOC_13_mgL = DOC_12_mgL * (AF/(1-AF))) %>%
    select( -Max_12CO2_ppm, -ends_with(c('Integral', 'Baseline'))) %>%
    filter(!is.na(Sample_ID))


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
    mutate(date = as.Date(Collection_date, format = '%m/%d/%Y'))
write_csv(dd, 'data/DOC_sample_data.csv')

filter(dd, date > as.Date('2023-07-11')) %>%
    group_by(date, trt, leachate, site) %>%
    summarize(across(everything(), mean)) %>%

ggplot( aes(date, DOC_13_mgL, col = trt, lty = factor(leachate))) +
    geom_point() + geom_line(size = 2) +
    facet_wrap(.~site, scales = 'free') + theme_classic()

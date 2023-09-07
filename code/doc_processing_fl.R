library(tidyverse)

doc <- read_csv('data/Kirby_flatheadlake/DOC_flatheadlake.csv')%>%
    select(Sample_ID, sample_type, Collection_date, trt, leachate, DOC_mgL,
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
    geom_point()

coefs <- cal$coefficients

doc <- doc %>%
    mutate(DOC_12_mgL = CO2_12_Integral * coefs[2] + coefs[1],
           AF = deltoAF(Delta_13C),
           DOC_13_mgL = DOC_12_mgL * (AF/(1-AF))) %>%
    select( -Max_12CO2_ppm, -ends_with(c('Integral', 'Baseline')))


# Determine the concentration of leachate:

gluc <- doc %>% filter(sample_type == 'conc') %>%
    select(-trt, -Collection_date, -sample_type) %>%
    filter(!(Sample_ID %in% c('glucose_A_REP', 'lake water'))) %>%
    mutate(Sample_ID = substr(Sample_ID, 1, nchar(Sample_ID)-2))

AF_gluc <- mean(gluc$AF[gluc$Sample_ID == 'glucose'])
del_LW <- doc$Delta_13C[doc$Sample_ID == 'lake water']
conc_LW <- doc$DOC_mgL[doc$Sample_ID == 'lake water']

leachate_dilution <- 0.3/100

leach <- gluc %>%
    mutate(glucose_13C = DOC_12_mgL * (AF_gluc/(1-AF_gluc)),
           extra_13C = DOC_13_mgL - glucose_13C,
           leachate_13C = extra_13C/leachate_dilution) %>%
    filter(leachate ==1)

leachate_conc_mgL = mean(leach$leachate_13C)

stocks <- data.frame(carbon_stock = c('leachate', 'glucose', 'flathead_lake'),
                     delta_13C = c(Inf, AFtodel(AF_gluc), del_LW),
                     AF = c(1, AF_gluc, deltoAF(del_LW)),
                     conc_stock_mgL = c(leachate_conc_mgL, 300, conc_LW),
                     conc_sample_mgL = c(leachate_conc_mgL/100, 1, conc_LW))

write_csv(stocks, 'data/DOC_stock_concentrations.csv')

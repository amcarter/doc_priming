library(tidyverse)

del_lw = -15
del_glucose = -12
AF_lw = deltoAF(del_lw)
AF_glucose = deltoAF(del_glucose)

C_lw = 1.2
C_glucose = 1

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

dat <- read_csv('processed_CO2/experiment1_processed.csv')

dd <- dat %>% #glimpse()
    select(sample_datetime, time_interval, treatment,
            DIC_molL, DIC_mgL, DIC_13C_umolL, AF, delCO2)


t0 <- filter(dd, time_interval == 'T0') %>%
    summarize(across(.cols = c(sample_datetime, DIC_molL, DIC_13C_umolL, delCO2),
                     .fns = mean))

end <- dd %>% filter( time_interval %in% c('T6', 'T5')) %>%
    select(-DIC_mgL, -AF, -delCO2, -time_interval)%>%
    mutate(delta_t = as.numeric(sample_datetime - t0$sample_datetime),
           delta_DIC_umolL = 10^6 * (DIC_molL - t0$DIC_molL)/delta_t,
           delta_DIC_13C_umolL = (DIC_13C_umolL - t0$DIC_13C_umolL)/delta_t,
           delta_DIC_12C_umolL = (delta_DIC_umolL - delta_DIC_13C_umolL)/delta_t)


ggplot(dd, aes(sample_datetime, DIC_molL*10^6, col = treatment))+
    geom_point()
ggplot(dd, aes(sample_datetime, DIC_13C_umolL, col = treatment))+
    geom_point()
ggplot(dd, aes(sample_datetime, delCO2, col = treatment))+
    geom_point()


lakewater_breakdown <- end %>%
    filter(treatment == 'LW') %>%
    mutate(lw_bd_umolCL = delta_DIC_12C_umolL/(1-AF_lw)) %>%
    select(lw_bd_umolCL)
lakewater_breakdown = mean(lakewater_breakdown$lw_bd_umolCL)


dd %>% filter(time_interval %in% c('T6', 'T5')) %>%
    group_by(treatment) %>%
    summarize(across(.cols = c(sample_datetime, DIC_molL, DIC_13C_umolL, delCO2),
                     .fns = mean)) %>%
    ungroup() %>%
    mutate(delta_DIC_umolL = (DIC_molL - t0$DIC_molL)*10^6,
           delta_DIC_13C = DIC_13C_umolL - t0$DIC_13C_umolL,
           delta_time = sample_datetime - t0$sample_datetime)

cw_12 <- 15
cs_12 <- 15
cl_13 <- 0.15

DICtrt_12 <- cw_12 + cs_12 + t0$DIC_12C_umolL
DICtrt_13 <- cw_12 * deltoAF(-27) + cs_12 * deltoAF(-12) + cl_13 + t0$DIC_13C_umolL

delDIC13 <- AFtodel(DICtrt_13/(DICtrt_12 + DICtrt_13))


print(paste('total DIC = ', DICtrt_12 + DICtrt_13))
print(paste('del13C = ', delDIC13))

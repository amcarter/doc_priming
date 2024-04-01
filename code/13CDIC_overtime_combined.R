library(tidyverse)
library(lme4)

fl <- read_csv('data/DIC_overtime_fl.csv')
fl2 <- read_csv('data/DIC_overtime_fl2.csv')
ybc <- read_csv('data/DIC_overtime_ybc.csv')
nsc <- read_csv('data/DIC_overtime_nsc.csv')
nmc <- read_csv('data/DIC_overtime_nmc.csv')

dat_sum <- bind_rows(fl2, ybc, fl, nsc, nmc)

p12 <- ggplot(dat_sum, aes(sample_datetime, DIC_mgL_mean, col = trt, lty = factor(leachate)))+
    geom_point() + geom_line(size = 1) +
    geom_errorbar(aes(ymin = DIC_mgL_mean - DIC_mgL_sd,
                      ymax = DIC_mgL_mean + DIC_mgL_sd), lty = 1)+
    facet_wrap(.~site, scales = "free") +
    xlab('Date') +theme_bw()

png(filename = 'figures/13CO2overtime_final.png',
width = 8, height = 5, res = 300, units = 'in')

ggplot(dat_sum, aes(sample_datetime, DIC_13C_ugL_mean, col = trt, lty = factor(leachate)))+
    geom_point() + geom_line(size = 1) +
    geom_errorbar(aes(ymin = DIC_13C_ugL_mean - DIC_13C_ugL_sd,
                      ymax = DIC_13C_ugL_mean + DIC_13C_ugL_sd), lty = 1)+
    facet_wrap(.~site, scales = "free") +
    xlab('Date') +theme_bw()

dev.off()

# png(filename = 'figures/ybc_expt_13CO2overtime.png',
    # width = 8, height = 10, res = 300, units = 'in')
ggpubr::ggarrange(p12, p13, common.legend = T)
# dev.off()

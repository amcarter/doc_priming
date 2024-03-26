library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)

fl <- read_csv("data/excessDIC_fl.csv")
nmc <- read_csv("data/excessDIC_nmc.csv")
nsc <- read_csv("data/excessDIC_nsc.csv")
fl2 <- read_csv("data/excessDIC_fl2.csv")
ybc <- read_csv("data/excessDIC_ybc.csv")
DOC1 <- read_csv("data/excessDOC_summer.csv")

# Final multipane plot with all DOC and DIC data (in progress) ####

combined <- rbind(fl, fl2, nmc, nsc, ybc)
combined$carbon <- factor(combined$carbon, levels = c("No Glucose", "Glucose"))

DOC1$carbon <- factor(DOC1$carbon, levels = c("No Glucose", "Glucose"))

DICplot <- ggplot(combined, aes(carbon, excess_13C_DIC, col = nuts, group = nuts)) +
    geom_line() +
    geom_point() +
    # geom_errorbar(aes(ymin = excess_13C_DIC - excess_13C_sd,
    #                   ymax = excess_13C_DIC + excess_13C_sd)) +
    facet_grid(cols = vars(site)) +  # Facet by site
    ggtitle('Excess 13C DIC')  # Title for the entire plot
DOCplot <- ggplot(DOC1, aes(carbon, excess_13C_DOC, col = nuts, group = nuts)) +
    geom_line() +
    geom_point()+
    facet_grid(cols = vars(site)) +
    ggtitle('Excess 13C DOC')

png(filename = 'figures/DOCDIC_excess13C.png',
    width = 10, height = 5, res = 300, units = 'in')
ggarrange(DICplot, DOCplot, ncol = 1)
dev.off()

# Excess 12C DOC plot ####

png(filename = 'figures/Excess_12C_DOC.png',
    width = 7, height = 5, res = 300, units = 'in')
ggplot(DOC1, aes(carbon, excess_12C_DOC, col = nuts, group = nuts)) +
    geom_line() +
    geom_point() +
    facet_grid(cols = vars(site)) +
    ggtitle('Excess 12C DOC (ug)')
dev.off()

# Excess DIC plot ####

png(filename = 'figures/Excess_DIC.png',
    width = , height = 5, res = 300, units = 'in')
ggplot(combined, aes(carbon, excess_DIC, col = nuts, group = nuts)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = excess_DIC - excess_DIC_sd,
                      ymax = excess_DIC + excess_DIC_sd)) +
    facet_grid(cols = vars(site)) +
    ggtitle('Excess DIC (ug)')
dev.off()

# Model bottle_dif combined table ####

fl_bottledif <- read_csv("data/DIC_change_by_bottle_fl.csv")
nmc_bottledif <- read_csv("data/DIC_change_by_bottle_nmc.csv")
nsc_bottledif <- read_csv("data/DIC_change_by_bottle_nmc.csv")

summer_bottledif <- rbind(fl_bottledif, nmc_bottledif, nsc_bottledif)
write_csv(summer_bottledif, 'data/summer_bottledif.csv')

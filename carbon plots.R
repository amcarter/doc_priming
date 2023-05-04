library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

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

setwd("/Users/niko/Dropbox/leachate_linder/processed_CO2")
experiment1_summary2 <- read_csv("~/Dropbox/leachate_linder/processed_CO2/experiment1_summary2.csv")
plot(experiment1_summary2$sample_datetime, experiment1_summary2$delCO2_norm, xlab="Time Interval",ylab="normalized delCO2")

sp <- ggplot(experiment1_summary2, aes(x = sample_datetime, y = CO2_conc, col = treatment))+
  geom_point()
sp

ggplot(experiment1_summary2, aes(x = sample_datetime, y = delCO2_norm, col = treatment))+
  geom_point()

filter(experiment1_summary2, treatment != 'LW_noleachate') %>% 
  ggplot( aes(x = sample_datetime, y = delCO2_norm, col = treatment))+
  geom_point()




#delCO2 plots
test <- experiment1_summary2 %>% 
  filter(treatment=="CNP" | time_interval=="T0")





x <- c(test$sample_datetime)
y <- c(test$delCO2_norm)
ggplot(test, aes(x,y), col = treatment)+
  geom_point()+
  geom_smooth(method = lm)

plot(test$rundate,test$delCO2)

plot(x,y)

lm(y ~ x)

x[12]-x[1]

#calculate rates and compare to literature values (deltaT in hours or days) micromol/L/day (multiply by volume of bottle)


#select() :- To select columns (variables)
#filter() :-To filter (subset) rows.
#mutate() :-To create new variables
#summarise() :- To summarize (or aggregate) data
#group_by() :- To group data
#arrange() :- To sort data



####Calculations#####
lake_DOC = 1.2
lake_del13C <- -25
lake_AF <- deltoAF(lake_del13C)

glu_DOC = 1
glu_del13C = -12
glu_AF = deltoAF(glu_del13C)

LC_DOC = lake_DOC + glu_DOC
LC_AF = (lake_AF * lake_DOC + glu_AF * glu_DOC)/LC_DOC
LC_del13C = AFtodel(LC_AF)
# calculate the fraction of CO2 from lake DOC:

exp1 <- experiment1_summary2 %>% 
  mutate(AF = deltoAF(delCO2_norm),
         CO2_13C_umol = CO2_conc *AF,
         CO2_12C_umol = CO2_conc * (1-AF),
         LC_CO2_13C_umol = case_when(treatment %in% c('LW', 'LW_noleachate', 'NP') ~
                                       (lake_AF * CO2_12C_umol)/(1-lake_AF),
                                     treatment %in% c('C', 'CNP') ~
                                       (LC_AF * CO2_12C_umol)/(1-LC_AF)),
         excess_CO2_13C_umol = CO2_13C_umol - LC_CO2_13C_umol)

t0 <- exp1 %>% filter(time_interval == 'T0')
t0_C <- t0 %>% mutate(treatment = 'C')
t0_CNP <- t0 %>% mutate(treatment = 'CNP')
t0_NP <- t0 %>% mutate(treatment = 'NP')

exp1 <- bind_rows(t0_C, t0_CNP, t0_NP, exp1)

exp1 %>%
  filter(treatment != 'LW_noleachate')%>%
ggplot(aes(sample_datetime, log(excess_CO2_13C_umol), col = treatment )) +
  geom_point()+
  geom_smooth(method = 'lm')
#lubridate::as.difftime()
exp1 <- exp1 %>% 
  mutate(time_days = as.numeric((sample_datetime - exp1$sample_datetime[1])/60/60/24))
# calculate decay coefficients ####
# LW 
LW <- exp1 %>% filter(treatment == 'LW')
mod_lin = lm(excess_CO2_13C_umol ~ time_days, LW)
mod_exp = lm(log(excess_CO2_13C_umol) ~ time_days, LW)
mod <- summary(mod_exp)
#summary(mod_lin)
k_LW <- mod$coefficients[2,1]
#LW_timedayscoef <- 5.359e-04

#C
C <- exp1 %>%  filter(treatment == 'C')
mod_lin = lm(excess_CO2_13C_umol ~ time_days, C)
mod_exp =lm(log(excess_CO2_13C_umol) ~ time_days, C)
mod <- summary(mod_exp)
#summary(mod_lin)
k_C <- mod$coefficients[2,1]
#C_timedayscoef <- 1.658e-03

#NP
NP <- exp1 %>%  filter(treatment == 'NP')
mod_lin = lm(excess_CO2_13C_umol ~ time_days, NP)
mod_exp =lm(log(excess_CO2_13C_umol) ~ time_days, NP)
mod <- summary(mod_exp)
k_NP <- mod$coefficients[2,1]
#summary(mod_lin)
#NP_timedayscoef <- 0.0016262

#CNP
CNP <- exp1 %>%  filter(treatment == 'CNP')
mod_lin = lm(excess_CO2_13C_umol ~ time_days, CNP)
mod_exp =lm(log(excess_CO2_13C_umol) ~ time_days, CNP)
mod <- summary(mod_exp)
k_CNP <- mod$coefficients[2,1]
#summary(mod_lin)
#CNP_timedayscoef <- 0.0036243

coef_data <- tibble(c(LW_timedayscoef,C_timedayscoef,NP_timedayscoef,CNP_timedayscoef))
View(coef_data)  

k_data <- tibble(c(k_LW,k_C,k_NP,k_CNP))
View(k_data)

names(k_data)[1] <- 'k_value'
 
##coef_data <- coef_data %>% add_column(new_col = c('LW', 'C', 'NP', 'CNP'))
#coef_data %>% 
  #rename(treatment = new_col)
#coef_data %>% 
  #ggplot(height = coef_data$linear_coefficient, width = 1, )

#####barplots####

# Create data
k_data$k_value[1]
k_dataframe <- data.frame(
  treatment = c('LW', 'C', 'NP', 'CNP') ,
  value = k_data$k_value ,
  normalized_value = k_data$k_value/k_data$k_value[1])
#normalized value is same as enhancement factor

  ggplot(k_dataframe, aes(x=treatment, y=value)) + 
  geom_bar(stat = "identity")

#normalized data plot for LW
ggplot(k_dataframe, aes(x=treatment, y=normalized_value)) + 
  geom_bar(stat = "identity")

#ggplot(datacoef, aes(x=treatment, y=normalized_value)) + 
  #geom_bar(stat = "identity")

##make it colored by treatment
coul <- brewer.pal(5, "Set2")
barplot(height=k_dataframe$normalized_value, names=k_dataframe$treatment, 
        col=coul,
        xlab = 'Treatment',
        ylab = 'Priming Enhancement Factor Compared to Control',
        main = 'Priming Effect by Treatment',
        ylim = c(0,4)
        )
View(k_dataframe)



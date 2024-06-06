# plots of hormone levels measured by GC-MS after DHT, DHT+ENZ, Fulvestrant in vitro cultures
# Petter Brodin, April 21th, 2024

library(tidyverse)
library(ggpubr)
library(stringr)
library(multcomp)
library(ggrepel)

#load data and preprocess
setwd("./data/ExtFigure4/")
d <- as.data.frame(read.csv("GAHT_in_vitro_240328eAN_HR.csv", header=T))

#Change gene formats to factors and numbers
d$Donor <- as.factor(d$Donor)
d$Pretreatment <- as.factor(d$Pretreatment) 
d$Stimulation <- as.factor(d$Stimulation) 
d$DHEA_pg_mL <- as.numeric(d$DHEA_pg_mL) 
d$Testo_pg_mL <- as.numeric(d$Testo_pg_mL) 
d$Estrone_pg_mL <- as.numeric(d$Estrone_pg_mL) 
d$DHT_pg_mL <- as.numeric(d$DHT_pg_mL) 
d$Androstendione_pg_mL <- as.numeric(d$Androstendione_pg_mL) 
d$Progesterone_pg_mL <- as.numeric(d$Progesterone_pg_mL) 
d$Estradiol_pg_mL <- as.numeric(d$Estradiol_pg_mL) 
d$X17alfa_Ohp_pg_mL <- as.numeric(d$X17alfa_Ohp_pg_mL) 

d.short <- d[, c("Donor", "Pretreatment", "Stimulation", "DHEA_pg_mL", "Testo_pg_mL", "Estrone_pg_mL", "DHT_pg_mL", "Androstendione_pg_mL", "X17alfa_Ohp_pg_mL")]

pivoted.dat <- pivot_longer(d.short, cols=!c('Donor', 'Pretreatment','Stimulation'), names_to = "hormone", values_to = "concentration_pg_ml")

--------------------------------------------
#For paper fig individual concentrations
plot1 <- ggplot(pivoted.dat, aes(x=hormone, y=log10(concentration_pg_ml), color=hormone)) + 
  geom_point(size=1, alpha=0.8) +
  theme_bw() + xlab("") + ylab("concentration (pg/ml)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ Pretreatment) +
  scale_colour_brewer(type="qual", palette = 2)
plot1


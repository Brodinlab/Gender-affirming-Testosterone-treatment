# Created by Camila Consiglio (camila.consiglio@med.lu.se)

## FIGURES 1B & 1C ##
library(ggplot2)
library(dplyr)


# meta
meta <- read.csv("Fig1_metadata.csv", row.names = 1)

# Bioavailable testosterone levels
kruskal.test(meta$BioavailTest~meta$VisitMonths) #p-value = 6.516e-10

G1 <- meta %>%
  ggplot(aes(x=VisitMonths, y=BioavailTest)) + 
  geom_rect(aes(ymin = 4.54, ymax = 23.64, xmin = -Inf, xmax = Inf), fill = "midnightblue", alpha = 0.006) +
  geom_point() +  
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  labs(x="Time (Months)", y="Bioavailable testosterone (nmol/L)", title="KWtest pvalue=6.516e-10") 
G1

# Estradiol levels
meta_E2 <- meta %>% 
  filter(!is.na(Estradiol)) %>% 
  filter(Estradiol < 900) #exclude samples that are being measured again

kruskal.test(meta_E2$Estradiol~meta_E2$VisitMonths) #p-value = 0.8101

G2 <- meta_E2 %>%
  ggplot(aes(x=VisitMonths, y=Estradiol)) + 
  geom_rect(aes(ymin = 110, ymax = 650, xmin = -Inf, xmax = Inf), fill = "deeppink4", alpha = 0.006) + 
  scale_y_continuous(limits = c(0,650)) +
  scale_x_continuous(breaks = c(0,3,12)) +
  theme_light() + #theme(legend.position = "none") +
  geom_path(aes(group=SubjectID), color="black", size=0.1) +
  geom_point() +  scale_shape_manual(values = c(16,21)) + 
  labs(x="Time (Months)", y="Estradiol (pmol/L)", title="KWtest pvalue=0.8101") 
G2
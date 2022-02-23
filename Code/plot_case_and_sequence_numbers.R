#plot cases over time with interventions

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory 
setwd("/Users/rhysinward/Documents/Master's_paper/EpiFilter")

#Load packages

library(tidyverse)
library(safejoin)
library(data.table)
library(lubridate)
library(dplyr)
library(zoo) 
library(countrycode) 
library(purrr) 
library(readr) 
library(stringr)   
library(tidyr) 
library(usdata)
library(seqinr)
library(ggplot2)
library(xts)
library(grid)

#Hong Kong

#cases

alldata = read.csv('Hong_Kong-COVID-19-data.csv')

#plot 

emergency <- grobTree(textGrob("State of\nEmergency Annouced", x=0.11,  y=0.88, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Border\nClosure", x=0.54,  y=0.88, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Ban on\nMass Gatherings", x=0.63,  y=0.88, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

alldata$Date_reported <- as.Date(alldata$Date_reported, '%d/%m/%Y')

HK_cases <- ggplot(alldata, aes(x = as.Date(Date_reported), y=New_cases)) + 
  geom_bar(stat = 'identity',fill="steelblue") +
  labs(x= "Time",y= 'Cases') +
  geom_vline(xintercept = as.Date("2020-01-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-25"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-29"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  scale_color_manual(values = c("deeppink", "royalblue")) + 
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15))

HK_cases

#sampling proportion

#Unsampled

load('unsampled.R')

freqs$names <- as.Date(freqs$names)

sampling_method <- grobTree(textGrob("Unsampled", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

unsampled <-ggplot(data=freqs, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,25) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

unsampled

#uniform

load('uniform.R')

uniform$names <- as.Date(uniform$names)

sampling_method <- grobTree(textGrob("Uniform", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

uniform <-ggplot(data=uniform, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,25) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

uniform

#proportional

load('proportional.R')

proportional$names <- as.Date(proportional$names)

sampling_method <- grobTree(textGrob("Proportional", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

proportional <-ggplot(data=proportional, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,25) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

proportional

#inverse

load('inverse.R')

inverse$names <- as.Date(inverse$names)

sampling_method <- grobTree(textGrob("Reciprocal-proportional", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

inverse <-ggplot(data=inverse, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,25) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

inverse

ggarrange(
  unsampled, proportional, uniform, inverse, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

#Brazil 

#Proportion of cases P.1

P1 <- read.csv ('P.1_proportion.csv')
P1$Group.1 <- as.Date (c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31',
                         '2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31'))

br_proportion <- ggplot(P1, aes(x= Group.1, y = x, fill = Identity
)) + geom_bar(stat = 'identity', position = 'fill') +
  labs(x="Date of Collection",y="Percentage of Sequences P.1/Gamma", fill = 'Varient') + theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

br_proportion

#cases

br_cases <- read.csv ('SIVEP_hospital_21-06-2021_SP_man.csv')

amazon <- filter(br_cases, Municip == 'MANAUS')

#Filter up to 7th feb

amazon$DateSymptoms <- as.Date(amazon$DateSymptoms)

amazon <- filter(amazon, DateSymptoms <=  '2021-02-07')

amazon_sum <- amazon %>% group_by (DateSymptoms) %>%
  summarize(New_cases = n())


#plot 



first <- grobTree(textGrob("First\nlockdown", x=0.16,  y=0.85, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
mask <- grobTree(textGrob("Mandatory\nuse of masks", x=0.31,  y=0.85, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
distancing <- grobTree(textGrob("Physical \ndistancing eased", x=0.49,  y=0.85, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))
elections <- grobTree(textGrob("Local\nelections", x=0.71,  y=0.85, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
second <- grobTree(textGrob("Second\nlockdown", x=0.80,  y=0.85, hjust=0,
                                   gp=gpar(col="Black", fontsize=12, fontface="italic")))

BR_cases <- ggplot(amazon_sum, aes(x = as.Date(DateSymptoms), y=New_cases)) + 
  geom_bar(stat = 'identity',fill="steelblue") +
  labs(x= "Time",y= 'Cases') +
  geom_vline(xintercept = as.Date("2020-03-17"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-06-02"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-07-06"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-11-15"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-25"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "2 month", date_labels = "%Y %b %d") +
  scale_color_manual(values = c("deeppink", "royalblue")) + 
  annotation_custom(first) +
  annotation_custom(mask) +
  annotation_custom(distancing) +
  annotation_custom(elections) +
  annotation_custom(second) +
  theme(text = element_text(size=15))

BR_cases

#sequences

#Unsampled

load('unsampled_br.R')

freqs$names <- as.Date(freqs$names)

sampling_method <- grobTree(textGrob("Unsampled", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

unsampled <-ggplot(data=freqs, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

unsampled

#uniform

load('uniform_br.R')

uniform$names <- as.Date(uniform$names)

sampling_method <- grobTree(textGrob("Uniform", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

uniform <-ggplot(data=uniform, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

uniform

#proportional

load('proportional_br.R')

by_case_sample$names <- as.Date(by_case_sample$names)

sampling_method <- grobTree(textGrob("Proportional", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

proportional <-ggplot(data=by_case_sample, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

proportional

#inverse

load('inverse_br.R')

inverse$names <- as.Date(inverse$names)

sampling_method <- grobTree(textGrob("Reciprocal-proportional", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))

inverse <-ggplot(data=inverse, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="steelblue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,25) +
  annotation_custom(sampling_method) +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") 

inverse

ggarrange(
  unsampled, proportional, uniform, inverse, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

  
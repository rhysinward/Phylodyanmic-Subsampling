######################################################################
## Apply EpiFilter to COVID data from Hong Kong
# From: Parag, KV, (2020) "Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves" BioRxiv.
######################################################################

# Notes and assumptions
# - load COVID incidence curves from WHO dashboard
# - estimate reproduction numbers with EpiFilter 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory 
setwd("/Users/rhysinward/Documents/Master's_paper/EpiFilter")

# Folder path for results
folres = paste0("./results/covid/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}
# Load data from WHO for New Zealand
alldata = read.csv('Hong_Kong-COVID-19-data.csv')
idcountry = which(alldata$Country == 'Hong Kong')

# Incidence and dates
Iday = alldata$New_cases[idcountry]
dates  = alldata$Date_reported[idcountry]
# Time series lengths
nday = length(dates); tday = 1:nday

# Approximate serial interval distribution from Ferguson et al
wdist = dgamma(tday, shape = 2.3669, scale = 2.7463)

# Total infectiousness
Lday = rep(0, nday) 
for(i in 2:nday){
  # Total infectiousness
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
}

######################################################################
## EpiFilter: provides formally smoothed and exact estimates
# Method based on Bayesian recursive filtering and smoothing
######################################################################

# Setup grid and noise parameters
Rmin = 0.01; Rmax = 10; eta = 0.1

# Uniform prior over grid of size m
m = 200; pR0 = (1/m)*rep(1, m)
# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)

# Filtered (causal) estimates as list [Rmed, Rhatci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredict(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025)
#have an extra argument for twice maximum of incidence 
# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025)

#smoothed growth rate predictions from filtered estimates 
#find the gamma parameters 
#mean serial interval = 6.6 with variation of 18.3
#alpha = 3.8, beta = 0.64
GR <- 0.26*(Rsmooth[[3]]**(1/1.87)-1)
GR
GRHPD <- 0.26*(Rsmooth[[2]]**(1/1.87)-1)
plot (GR)

#####EpiFilter done #####

#compare BDSKY with Epifilter 

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

# Set working directory 
setwd("/Users/rhysinward/Documents/Master's_paper/Data/Hong_Kong/BDSKY")


#Unsampled

load('unsampled/unsampled_new.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:106
colnames(time_series_BDSKY)[2] <- "BDSKY_mean"
time_series_BDSKY$mean <- rev(time_series_BDSKY$BDSKY_mean)
time_series_BDSKY$lb <- rev(time_series_BDSKY$`Lower Bound`)
time_series_BDSKY$ub <- rev(time_series_BDSKY$`Upper Bound`)
#HPD BDSKY
p <-ggplot(data=time_series_BDSKY, aes(x=time, y=mean)) + geom_point() + geom_line()
p <-p+geom_ribbon(aes(ymin=lb, ymax=ub), linetype=2, alpha=0.1)
p

#HPD epifilter
time_series_epifilter <- as.data.frame (Rsmooth[[2]][, 2:nday]) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (Rsmooth[[3]][2:nday]) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:106
q <-ggplot(data=time_series_epifilter, aes(x=time, y=Rsmooth[[3]][2:nday])) + geom_point() + geom_line()
q <-q+geom_ribbon(aes(ymin=V1, ymax=V2), linetype=2, alpha=0.1)
q
#combine datasets for epifilter and BDSKY 
#most likely will need to rbind them (i.e. make sure the rows align and then combine them)
#then will need to create another column with the group in it 

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter$`Rsmooth[[3]][2:nday]`
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'

#BDSKY 

bd <- as.data.frame(time_series_BDSKY$lb)
colnames(bd)[1] <- "lb"
bd$ub <- time_series_BDSKY$ub
bd$mean<- time_series_BDSKY$mean
bd$time <- time_series_BDSKY$time
bd$group <- 'BDSKY'

#check to see if confidence intervals overlap 

combined1 <- cbind(bd,ef)
combined1
rng = cbind(pmin(combined1[,1], combined1[,2]), pmax(combined1[,1], combined1[,2]),
            pmin(combined1[,6], combined1[,7]), pmax(combined1[,6], combined1[,7]))

olap = (rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3])

table (olap)


#add dates 

bd$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

ef$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

combined <- rbind(bd,ef)
combined$mean <- as.numeric(unlist(combined$mean))

# Create a text
sampling_method <- grobTree(textGrob("Unsampled", x=0.075,  y=0.97, hjust=0,
                          gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("State of\nEmergency Annouced", x=0.075,  y=0.88, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Border\nClosure", x=0.45,  y=0.88, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Ban on\nMass Gatherings", x=0.60,  y=0.88, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="italic")))

unsampled <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-01-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-25"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-29"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,4.0) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 

unsampled

#proportional

load('Proportional/proportional_new.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:106
colnames(time_series_BDSKY)[2] <- "BDSKY_mean"
time_series_BDSKY$mean <- rev(time_series_BDSKY$BDSKY_mean)
time_series_BDSKY$lb <- rev(time_series_BDSKY$`Lower Bound`)
time_series_BDSKY$ub <- rev(time_series_BDSKY$`Upper Bound`)
#HPD BDSKY
p <-ggplot(data=time_series_BDSKY, aes(x=time, y=mean)) + geom_point() + geom_line()
p <-p+geom_ribbon(aes(ymin=lb, ymax=ub), linetype=2, alpha=0.1)
p

#HPD epifilter
time_series_epifilter <- as.data.frame (Rsmooth[[2]][, 2:nday]) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (Rsmooth[[3]][2:nday]) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:106
library (ggplot2)
q <-ggplot(data=time_series_epifilter, aes(x=time, y=Rsmooth[[3]][2:nday])) + geom_point() + geom_line()
q <-q+geom_ribbon(aes(ymin=V1, ymax=V2), linetype=2, alpha=0.1)
q
#combine datasets for epifilter and BDSKY 
#most likely will need to rbind them (i.e. make sure the rows align and then combine them)
#then will need to create another column with the group in it 

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter$`Rsmooth[[3]][2:nday]`
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'

#BDSKY 

bd <- as.data.frame(time_series_BDSKY$lb)
colnames(bd)[1] <- "lb"
bd$ub <- time_series_BDSKY$ub
bd$mean<- time_series_BDSKY$mean
bd$time <- time_series_BDSKY$time
bd$group <- 'BDSKY'

#check to see if confidence intervals overlap 

combined1 <- cbind(bd,ef)
combined1
rng = cbind(pmin(combined1[,1], combined1[,2]), pmax(combined1[,1], combined1[,2]),
            pmin(combined1[,6], combined1[,7]), pmax(combined1[,6], combined1[,7]))

olap = (rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3])

table (olap)

#add - dates

bd$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

ef$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

combined <- rbind(bd,ef)
combined$mean <- as.numeric(unlist(combined$mean))

# Create a text
sampling_method <- grobTree(textGrob("Proportional", x=0.075,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("State of\nEmergency Annouced", x=0.075,  y=0.88, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Border\nClosure", x=0.45,  y=0.88, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Ban on\nMass Gatherings", x=0.60,  y=0.88, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

proportional <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-01-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-25"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-29"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,4.0) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 

proportional 

#uniform

load('Uniform/uniform_new.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:106
colnames(time_series_BDSKY)[2] <- "BDSKY_mean"
time_series_BDSKY$mean <- rev(time_series_BDSKY$BDSKY_mean)
time_series_BDSKY$lb <- rev(time_series_BDSKY$`Lower Bound`)
time_series_BDSKY$ub <- rev(time_series_BDSKY$`Upper Bound`)
#HPD BDSKY
p <-ggplot(data=time_series_BDSKY, aes(x=time, y=mean)) + geom_point() + geom_line()
p <-p+geom_ribbon(aes(ymin=lb, ymax=ub), linetype=2, alpha=0.1)
p

#HPD epifilter
time_series_epifilter <- as.data.frame (Rsmooth[[2]][, 2:nday]) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (Rsmooth[[3]][2:nday]) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:106
q <-ggplot(data=time_series_epifilter, aes(x=time, y=Rsmooth[[3]][2:nday])) + geom_point() + geom_line()
q <-q+geom_ribbon(aes(ymin=V1, ymax=V2), linetype=2, alpha=0.1)
q
#combine datasets for epifilter and BDSKY 
#most likely will need to rbind them (i.e. make sure the rows align and then combine them)
#then will need to create another column with the group in it 

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter$`Rsmooth[[3]][2:nday]`
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'

#BDSKY 

bd <- as.data.frame(time_series_BDSKY$lb)
colnames(bd)[1] <- "lb"
bd$ub <- time_series_BDSKY$ub
bd$mean<- time_series_BDSKY$mean
bd$time <- time_series_BDSKY$time
bd$group <- 'BDSKY'

#check to see if confidence intervals overlap 

combined1 <- cbind(bd,ef)
combined1
rng = cbind(pmin(combined1[,1], combined1[,2]), pmax(combined1[,1], combined1[,2]),
            pmin(combined1[,6], combined1[,7]), pmax(combined1[,6], combined1[,7]))

olap = (rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3])

table (olap)

#add dates 

bd$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

ef$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

combined <- rbind(bd,ef)
combined$mean <- as.numeric(unlist(combined$mean))

# Create a text
sampling_method <- grobTree(textGrob("Uniform", x=0.075,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("State of\nEmergency Annouced", x=0.075,  y=0.88, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Border\nClosure", x=0.45,  y=0.88, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Ban on\nMass Gatherings", x=0.60,  y=0.88, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

uniform <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-01-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-25"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-29"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,4.0) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 

uniform

#inverse

load('Inverse/inverse_new.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:106
colnames(time_series_BDSKY)[2] <- "BDSKY_mean"
time_series_BDSKY$mean <- rev(time_series_BDSKY$BDSKY_mean)
time_series_BDSKY$lb <- rev(time_series_BDSKY$`Lower Bound`)
time_series_BDSKY$ub <- rev(time_series_BDSKY$`Upper Bound`)
#HPD BDSKY
p <-ggplot(data=time_series_BDSKY, aes(x=time, y=mean)) + geom_point() + geom_line()
p <-p+geom_ribbon(aes(ymin=lb, ymax=ub), linetype=2, alpha=0.1)
p

#HPD epifilter
time_series_epifilter <- as.data.frame (Rsmooth[[2]][, 2:nday]) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (Rsmooth[[3]][2:nday]) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:106
q <-ggplot(data=time_series_epifilter, aes(x=time, y=Rsmooth[[3]][2:nday])) + geom_point() + geom_line()
q <-q+geom_ribbon(aes(ymin=V1, ymax=V2), linetype=2, alpha=0.1)
q
#combine datasets for epifilter and BDSKY 
#most likely will need to rbind them (i.e. make sure the rows align and then combine them)
#then will need to create another column with the group in it 

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter$`Rsmooth[[3]][2:nday]`
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'

#BDSKY 

bd <- as.data.frame(time_series_BDSKY$lb)
colnames(bd)[1] <- "lb"
bd$ub <- time_series_BDSKY$ub
bd$mean<- time_series_BDSKY$mean
bd$time <- time_series_BDSKY$time
bd$group <- 'BDSKY'

#check to see if confidence intervals overlap 

combined1 <- cbind(bd,ef)
combined1
rng = cbind(pmin(combined1[,1], combined1[,2]), pmax(combined1[,1], combined1[,2]),
            pmin(combined1[,6], combined1[,7]), pmax(combined1[,6], combined1[,7]))

olap = (rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3])

table (olap)

#add_dates

#add dates 

bd$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

ef$dates <- seq(as.Date('2020-01-27'), as.Date('2020-05-11'), by = "1 days")

combined <- rbind(bd,ef)
combined$mean <- as.numeric(unlist(combined$mean))

# Create a text
sampling_method <- grobTree(textGrob("Reciprocal-proportional", x=0.075,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("State of\nEmergency Annouced", x=0.075,  y=0.88, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Border\nClosure", x=0.45,  y=0.88, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Ban on\nMass Gatherings", x=0.60,  y=0.88, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

inverse <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-01-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-25"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-03-29"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,4.0) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 

inverse

#add JS distance 

sampling_scheme <- c("Unsampled", "Proportional", "Uniform", "Reciprocal-proportional")
js <- c(0.24, 0.13, 0.27, 0.28)

df <- data.frame(sampling_scheme, js)

df$sampling_scheme <- factor(df$sampling_scheme,levels = c("Proportional","Unsampled", "Uniform", "Reciprocal-proportional"))


hk_plot <- ggplot(df, aes(x = sampling_scheme,y = js )) +  
  geom_bar(stat = "identity" ) +
  theme(legend.position="none") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none") +
  labs(x="Sampling Scheme", y= "Jensen-Shannon Distance") 



library(ggpubr)
ggarrange(
  unsampled, proportional,uniform, inverse,ncol =2, nrow =2, labels = c("A", "B", 'C','D','E'),
  common.legend = TRUE
)

ggarrange (hk_plot, ncol =1, nrow =2,labels = c('E'),
           common.legend = FALSE
)



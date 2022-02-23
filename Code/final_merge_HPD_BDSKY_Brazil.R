######################################################################
## Apply EpiFilter to COVID data from New Zealand
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
alldata = read.csv('Brazil-COVID-19-data_updated.csv')
idcountry = which(alldata$Country == 'Brazil')

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

recursPredict1 <- function(Rgrid, pR, Lday, Rmean, a){
  
  # Grid size and length of time series
  nday = nrow(pR); m = ncol(pR)
  # Test lengths of inputs
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }
  
  # Mean prediction: Lday[i] => Iday[i+1]
  pred = Lday*Rmean; pred = pred[1:length(pred)-1]
  
  # Discrete space of possible predictions
  Igrid = 0:50000; lenI = length(Igrid);
  
  # Check if close to upper bound
  if (any(pred > 0.9*max(Igrid))){
    stop("Epidemic size too large")  
  }
  
  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predInt = matrix(0, 4, nday-1)
  
  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    # Compute rate from Poisson renewal
    rate = Lday[i]*Rgrid
    # Prob of any I marginalised over Rgrid
    pI = rep(0, lenI)
    
    # Probabilities of observations 1 day ahead
    for(j in 1:lenI){
      # Raw probabilities of Igrid
      pIset = dpois(Igrid[j], rate)
      # Normalised by probs of R
      pI[j] = sum(pIset*pR[i, ])
    }
    
    # Quantile predictions and CDF at i+1
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)
    
    # Assign prediction results
    predInt[1, i] = Igrid[id1[1]]; predInt[2, i] = Igrid[id2[1]]
    predInt[3, i] = Igrid[id3[1]]; predInt[4, i] = Igrid[id4[1]]
  }
  # Main outputs: mean and 95% predictions
  recursPredict = list(pred, predInt)
}
Ifilt = recursPredict1(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025)

#have an extra argument for twice maximum of incidence 
# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict1(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025)

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
setwd("/Users/rhysinward/Documents/Master's_paper/Data/Brazil/BDSKY")


#Unsampled

load('Unsampled/BDSKY_Trail_2/unsampled.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:69
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
time_series_epifilter$time <- 1:69
library (ggplot2)
q <-ggplot(data=time_series_epifilter, aes(x=time, y=Rsmooth[[3]][2:nday])) + geom_point() + geom_line()
q <-q+geom_ribbon(aes(ymin=V1, ymax=V2), linetype=2, alpha=0.1)
q
#combine datasets for epifilter and BDSKY 
#most likely will need to rbind them (i.e. make sure the rows align and then combine them)
#then will need to create another column with the group in it 

#making sure the rows align

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter$`Rsmooth[[3]][2:nday]`
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

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

bd$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))


#combine

combined <- rbind(bd,ef)

df.new = combined[seq(1, nrow(combined), 3), ]



#plot graph 

sampling_method <- grobTree(textGrob("Unsampled", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("Suspension of\nCommercial \nActivities", x=0.18,  y=0.85, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Resumption of\nCommercial \nActivities", x=0.44,  y=0.85, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Restrictions \nre-introducted", x=0.63,  y=0.85, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

unsampled <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,5.5) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 

unsampled

#proportional

load('Proportional/trial_2/proportional.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:69
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
time_series_epifilter$time <- 1:69
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

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

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

bd$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))

#combine

combined <- rbind(bd,ef)

df.new = combined[seq(1, nrow(combined), 3), ]

#plot graph 

sampling_method <- grobTree(textGrob("Proportional", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("Suspension of\nCommercial \nActivities", x=0.18,  y=0.85, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Resumption of\nCommercial \nActivities", x=0.44,  y=0.85, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Restrictions \nre-introducted", x=0.63,  y=0.85, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

proportional <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,5.5) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 


proportional

#uniform

load('Uniform/trail_2/uniform.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:69
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
time_series_epifilter$time <- 1:69
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

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

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

bd$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))

#combine

combined <- rbind(bd,ef)

df.new = combined[seq(1, nrow(combined), 3), ]


#plot graph 

sampling_method <- grobTree(textGrob("Uniform", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("Suspension of\nCommercial \nActivities", x=0.18,  y=0.85, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Resumption of\nCommercial \nActivities", x=0.44,  y=0.85, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Restrictions \nre-introducted", x=0.63,  y=0.85, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

uniform <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,5.5) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15))

uniform

#inverse

load('Inverse/trail_2/inverse.R')

time_series_BDSKY <- as.data.frame (Re_gridded_hpd)
time_series_BDSKY <- as.data.frame(t(time_series_BDSKY))
time_series_BDSKY$time <- 1:69
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
time_series_epifilter$time <- 1:69
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

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

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

bd$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30,','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))

#combine

combined <- rbind(bd,ef)

df.new = combined[seq(1, nrow(combined), 3), ]


#plot graph 

sampling_method <- grobTree(textGrob("Reciprocal-proportional", x=0.05,  y=0.97, hjust=0,
                                     gp=gpar(col="Black", fontsize=12, fontface="bold")))
emergency <- grobTree(textGrob("Suspension of\nCommercial \nActivities", x=0.18,  y=0.85, hjust=0,
                               gp=gpar(col="Black", fontsize=12, fontface="italic")))
border <- grobTree(textGrob("Resumption of\nCommercial \nActivities", x=0.44,  y=0.85, hjust=0,
                            gp=gpar(col="Black", fontsize=12, fontface="italic")))
gatherings <- grobTree(textGrob("Restrictions \nre-introducted", x=0.63,  y=0.85, hjust=0,
                                gp=gpar(col="Black", fontsize=12, fontface="italic")))

inverse <- ggplot(combined, aes(x = dates, y=mean, color = group, fill = group)) + 
  geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.3)) +
  scale_fill_manual(values=c("royalblue", "deeppink"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(0,5.5) + 
  scale_color_manual(values = c("royalblue", "deeppink")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15))

inverse

#add JS distance 

sampling_scheme <- c("Unsampled", "Proportional", "Uniform", "Reciprocal-proportional")
js <- c(0.53, 0.42, 0.40, 0.43)

df <- data.frame(sampling_scheme, js)

df$sampling_scheme <- factor(df$sampling_scheme,levels = c("Uniform", "Proportional",  "Reciprocal-proportional","Unsampled"))


hk_plot <- ggplot(df, aes(x = sampling_scheme,y = js)) +  
  geom_bar(stat = "identity" ) +
  theme(legend.position="none") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none") +
  labs(x="Sampling Scheme", y= "Jensen-Shannon Distance") 

library(ggpubr)
ggarrange(
  unsampled, proportional, uniform, inverse, labels = c("A", "B", 'C','D'),
  common.legend = TRUE
)

ggarrange (hk_plot, ncol =1, nrow =2,labels = c('E'),
           common.legend = FALSE
)


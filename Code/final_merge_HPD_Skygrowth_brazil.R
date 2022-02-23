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

#Extract the R0 values and HPD 

library (ggplot2)

time_series_epifilter <- as.data.frame (Rsmooth[[2]][2:nday]) # Median Values

class(time_series_epifilter$`Rsmooth[[2]][2:nday]`)

#####EpiFilter done #####

#compare skygrowth with Epifilter 

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
setwd("/Users/rhysinward/Documents/Master's_paper/Data/Brazil/skygrowth")

#skygrowth

#unsampled

load('Unsampled/unsamapled_skygrowth.R')
q <- as.data.frame(growth[["data"]][["med"]])
q$lb <- as.data.frame(growth[["data"]][["lb"]])
q$ub <- as.data.frame(growth[["data"]][["ub"]])
q <- 0.26*(q**(1/1.87)-1)
q <- q[1:6,]


time_series_epifilter <- as.data.frame (GRHPD) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (GR) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:70

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter1$GR
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'
ef <- ef[1:69,]

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

#skygrowth

sy <- as.data.frame(q$`growth[["data"]][["lb"]]`)
colnames(sy)[1] <- "lb"
sy$ub <- q$`growth[["data"]][["ub"]]`
colnames(sy)[2] <- "ub"
sy$mean <- q$`growth[["data"]][["med"]]`
colnames(sy)[3] <- "mean"
sy <- sy[1:6,]
sy$time <- c(1,14,28,42,56,69)
sy$group <- 'Skygrowth'


#unsampled alone

unsampledsg <- ggplot(sy, aes(x = time, y=mean)) + 
  geom_line()+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.6)) +
  labs(x= "Time",y="rT")+
  theme_bw()

unsampledsg

#check to see if confidence intervals overlap 
df <- sy
df <- df[1:7,]
df <- df[rep(seq_len(nrow(df)), each = 13), ]

#add_dates

sy$dates <- as.Date( c('2020-11-30','2020-12-13',"'2020-12-27'",'2021-01-10','2021-01-24','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'))

df <- ef %>% slice(which(row_number() %% 2 == 1))


combined <- rbind(sy,df)

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
  scale_fill_manual(values=c("deeppink", "royalblue"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(-0.2,1)  + 
  scale_color_manual(values = c("deeppink", "royalblue")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 

unsampled

#uniform

load('Uniform/uniform_skygrowth.R')
q <- as.data.frame(growth[["data"]][["med"]])
q$lb <- as.data.frame(growth[["data"]][["lb"]])
q$ub <- as.data.frame(growth[["data"]][["ub"]])
q <- 0.26*(q**(1/1.87)-1)
q <- q[1:6,]
q$time <- 1:6

time_series_epifilter <- as.data.frame (GRHPD) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (GR) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:70

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter1$GR
ef$time <- time_series_epifilter$time
ef$group <- 'epifilter'
ef <- ef[1:69,]

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))


#skygrowth

sy <- as.data.frame(q$`growth[["data"]][["lb"]]`)
colnames(sy)[1] <- "lb"
sy$ub <- q$`growth[["data"]][["ub"]]`
colnames(sy)[2] <- "ub"
sy$mean <- q$`growth[["data"]][["med"]]`
colnames(sy)[3] <- "mean"
sy <- sy[1:6,]
sy$time <- c(1,14,28,42,56,69)
sy$group <- 'skygrowth'


#unsampled alone

unsampledsg <- ggplot(sy, aes(x = time, y=mean)) + 
  geom_line()+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.6)) +
  labs(x= "Time",y="rT")+
  theme_bw()

unsampledsg

#check to see if confidence intervals overlap 
df <- sy
df <- df[1:7,]
df <- df[rep(seq_len(nrow(df)), each = 13), ]

#add_dates

sy$dates <- as.Date( c('2020-11-30','2020-12-13',"'2020-12-27'",'2021-01-10','2021-01-24','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'
                        ))

df <- ef %>% slice(which(row_number() %% 2 == 1))


combined <- rbind(sy,df)

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
  scale_fill_manual(values=c("deeppink", "royalblue"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(-0.2,1)  + 
  scale_color_manual(values = c("deeppink", "royalblue")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15))


uniform

#proportionl

load('Proportional/proportional_skygrowth.R')
q <- as.data.frame(growth[["data"]][["med"]])
q$lb <- as.data.frame(growth[["data"]][["lb"]])
q$ub <- as.data.frame(growth[["data"]][["ub"]])
q <- 0.26*(q**(1/1.87)-1)
q <- q[1:6,]
q$time <- 1:6

time_series_epifilter <- as.data.frame (GRHPD) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (GR) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:70

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter1$GR
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'
ef <- ef[1:69,]

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

#skygrowth

sy <- as.data.frame(q$`growth[["data"]][["lb"]]`)
colnames(sy)[1] <- "lb"
sy$ub <- q$`growth[["data"]][["ub"]]`
colnames(sy)[2] <- "ub"
sy$mean <- q$`growth[["data"]][["med"]]`
colnames(sy)[3] <- "mean"
sy <- sy[1:6,]
sy$time <- c(1,14,28,42,56,69)
sy$group <- 'Skygrowth'


#unsampled alone

unsampledsg <- ggplot(sy, aes(x = time, y=mean)) + 
  geom_line()+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.6)) +
  labs(x= "Time",y="rT")+
  theme_bw()

unsampledsg

#check to see if confidence intervals overlap 
df <- sy
df <- df[1:7,]
df <- df[rep(seq_len(nrow(df)), each = 13), ]

#add_dates

sy$dates <- as.Date( c('2020-11-30','2020-12-13',"'2020-12-27'",'2021-01-10','2021-01-24','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'
                       ))

df <- ef %>% slice(which(row_number() %% 2 == 1))


combined <- rbind(sy,df)

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
  scale_fill_manual(values=c("deeppink", "royalblue"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(-0.2,1)  + 
  scale_color_manual(values = c("deeppink", "royalblue")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15)) 


proportional

#inverse

load('Inverse_proportional/inverse_skygrowth.R')
q <- as.data.frame(growth[["data"]][["med"]])
q$lb <- as.data.frame(growth[["data"]][["lb"]])
q$ub <- as.data.frame(growth[["data"]][["ub"]])
q <- 0.26*(q**(1/1.87)-1)
q <- q[1:6,]
q$time <- 1:6

time_series_epifilter <- as.data.frame (GRHPD) # Median Values
time_series_epifilter <- as.data.frame(t(time_series_epifilter))
time_series_epifilter1 <- as.data.frame (GR) # Median Values
time_series_epifilter <- cbind(time_series_epifilter,time_series_epifilter1)
time_series_epifilter$time <- 1:70

#making sure the rows allign

#epifilter
ef <- as.data.frame (time_series_epifilter$V1)
colnames(ef)[1] <- "lb"
ef$ub <- time_series_epifilter$V2
ef$mean <- time_series_epifilter1$GR
ef$time <- time_series_epifilter$time
ef$group <- 'EpiFilter'
ef <- ef[1:69,]

#rolling average 

ef$mean <- as.numeric(rollapply(ef[,3],5,fill=NA,partial = TRUE,FUN=mean))
ef$lb <- as.numeric(rollapply(ef[,1],5,fill=NA,partial = TRUE,FUN=mean))
ef$ub <-  as.numeric(rollapply(ef[,2],5,fill=NA,partial = TRUE,FUN=mean))

#skygrowth

sy <- as.data.frame(q$`growth[["data"]][["lb"]]`)
colnames(sy)[1] <- "lb"
sy$ub <- q$`growth[["data"]][["ub"]]`
colnames(sy)[2] <- "ub"
sy$mean <- q$`growth[["data"]][["med"]]`
colnames(sy)[3] <- "mean"
sy <- sy[1:6,]
sy$time <- c(1,14,28,42,56,69)
sy$group <- 'Skygrowth'


#unsampled alone

unsampledsg <- ggplot(sy, aes(x = time, y=mean)) + 
  geom_line()+
  geom_ribbon(aes(ymin=lb, ymax=ub, alpha=0.6)) +
  labs(x= "Time",y="rT")+
  theme_bw()

unsampledsg

#check to see if confidence intervals overlap 
df <- sy
df <- df[1:7,]
df <- df[rep(seq_len(nrow(df)), each = 13), ]

#add_dates

sy$dates <- as.Date( c('2020-11-30','2020-12-13',"'2020-12-27'",'2021-01-10','2021-01-24','2021-02-06'))
ef$dates <- as.Date( c('2020-11-30','2020-12-01','2020-12-02','2020-12-03','2020-12-04','2020-12-05',
                       '2020-12-06','2020-12-07','2020-12-08','2020-12-09','2020-12-10','2020-12-11','2020-12-12',
                       '2020-12-13','2020-12-14','2020-12-15','2020-12-16','2020-12-17','2020-12-18','2020-12-19',
                       '2020-12-20','2020-12-21','2020-12-22','2020-12-23','2020-12-24','2020-12-25','2020-12-26',
                       '2020-12-27','2020-12-28','2020-12-29','2020-12-30','2020-12-31','2021-01-01','2021-01-02',
                       '2021-01-03','2021-01-04','2021-01-05','2021-01-06','2021-01-07','2021-01-08','2021-01-09',
                       '2021-01-10','2021-01-11','2021-01-12','2021-01-13','2021-01-14','2021-01-15','2021-01-16',
                       '2021-01-17','2021-01-18','2021-01-19','2021-01-20','2021-01-21','2021-01-22','2021-01-23',
                       '2021-01-24','2021-01-25','2021-01-26','2021-01-27','2021-01-28','2021-01-29','2021-01-30',
                       '2021-01-31','2021-02-01','2021-02-02','2021-02-03','2021-02-04','2021-02-05','2021-02-06'
                       ))

df <- ef %>% slice(which(row_number() %% 2 == 1))

combined <- rbind(df,sy)

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
  scale_fill_manual(values=c("deeppink", "royalblue"), name="Model") +
  guides(alpha = "none") +
  guides(color = "none") +
  labs(x= "Time",y= expression(R[t]), fill = "Model") +
  geom_vline(xintercept = as.Date("2020-12-23"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-12-28"), color = "black", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-01-12"), color = "black", linetype = "dashed") +
  theme_bw() +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  ylim(-0.2,1)  + 
  scale_color_manual(values = c("deeppink", "royalblue")) + 
  annotation_custom(sampling_method) +
  annotation_custom(emergency) +
  annotation_custom(border) +
  annotation_custom(gatherings) +
  theme(text = element_text(size=15))

inverse

#add JS distance 

sampling_scheme <- c("Unsampled", "Proportional", "Uniform", "Reciprocal-proportional")
js <- c(0.57, 0.47, 0.45, 0.46)

df <- data.frame(sampling_scheme, js)

df$sampling_scheme <- factor(df$sampling_scheme,levels = c("Uniform", "Reciprocal-proportional", "Proportional", "Unsampled"))


hk_plot <- ggplot(df, aes(x = sampling_scheme,y = js)) +  
  geom_bar(stat = "identity" ) +
  theme(legend.position="none") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none") +
  labs(x="Sampling Scheme", y= "Jensen-Shannon Distance") 


library(ggpubr)
ggarrange(unsampled, proportional, uniform, inverse, labels = c("A", "B", 'C','D'),
  common.legend = TRUE)

ggarrange (hk_plot, ncol =1, nrow =2,labels = c('E'),
           common.legend = FALSE
)

#convert meta data to epi week
#install.packages('lubridate')
library(lubridate)

#epiweeks as according to this
#https://www.cmmcp.org/sites/g/files/vyhlif2966/f/uploads/epiweekcalendar2020.pdf



get_epiweek_starts <- function() {
  
  epiweek_jan <- c("2019-12-29","2020-01-05","2020-01-12","2020-01-19","2020-01-26")
  epiweek_feb <- paste("2020-02-",c("02","09","16","23"),sep="")
  epiweek_mar <- paste("2020-03-",c("01","08","15","22","29"),sep="")
  epiweek_apr <- paste("2020-04-",c("05","12","19","26"),sep="")
  epiweek_may <- paste("2020-05-",c("03","10","17","24","31"),sep="")
  epiweek_jun <- paste("2020-06-",c("07","14","21","28"),sep="")
  epiweek_jul <- paste("2020-07-",c("05","12","19","26"),sep="")
  epiweek_aug <- paste("2020-08-",c("02","09","16","23","30"),sep="")
  epiweek_sep <- paste("2020-09-",c("06","13","20","27"),sep="")
  epiweek_oct <- paste("2020-10-",c("04","11","18","25"),sep="")
  epiweek_nov <- paste("2020-11-",c("01","08","15","22","29"),sep="")
  epiweek_dec <- paste("2020-12-",c("06","13","20","27"),sep="")
  epiweek_jan21 <- paste("2021-01-",c("03","10","17","24","31"),sep="")
  epiweek_feb21 <- paste("2021-02-",c("07","14","21","28"),sep="")
  epiweek_mar21 <- paste("2021-03-",c("07","14","21","28"),sep="")
  epiweek_apr21 <- paste("2021-04-",c("04","11","18","25"),sep="")
  epiweek_may21 <- paste("2021-05-",c("02","09","16","23","30"),sep="")
  epiweek_jun21 <- paste("2021-06-",c("06","13","20","27"),sep="")
  
  epiweek_starts <- c(epiweek_jan,epiweek_feb,epiweek_mar,epiweek_apr,epiweek_may,epiweek_jun,
                      epiweek_jul,epiweek_aug,epiweek_sep,epiweek_oct,epiweek_nov,epiweek_dec,
                      epiweek_jan21,epiweek_feb21,epiweek_mar21,epiweek_apr21,epiweek_may21,
                      epiweek_jun21)
  return( epiweek_starts)
}

get_epiweek_ends <- function() {
  
  epiweek_jan <- c(paste("2020-01-",c("04","11","18","25"),sep=""),"2020-02-01")
  epiweek_feb <- paste("2020-02-",c("08","15","22","29"),sep="")
  epiweek_mar <- c(paste("2020-03-",c("07","14","21","28"),sep=""),"2020-04-04")
  epiweek_apr <- c(paste("2020-04-",c("11","18","25"),sep=""),"2020-05-02")
  epiweek_may <- c(paste("2020-05-",c("09","16","23","30"),sep=""),"2020-06-06")
  epiweek_jun <- c(paste("2020-06-",c("13","20","27"),sep=""),"2020-07-04")
  epiweek_jul <- c(paste("2020-07-",c("11","18","25"),sep=""),"2020-08-01")
  epiweek_aug <- c(paste("2020-08-",c("08","15","22","29"),sep=""),"2020-09-05")
  epiweek_sep <- c(paste("2020-09-",c("12","19","26"),sep=""),"2020-10-03")
  epiweek_oct <- paste("2020-10-",c("10","17","24","31"),sep="")
  epiweek_nov <- c(paste("2020-11-",c("07","14","21","28"),sep=""),"2020-12-05")
  epiweek_dec <- c(paste("2020-12-",c("12","19","26"),sep=""),"2021-01-02")
  epiweek_jan21 <- paste("2021-01-",c("09","16","23","30"),sep="")
  epiweek_feb21 <- paste("2021-02-",c("06","13","20","27"),sep="")
  epiweek_mar21 <- paste("2021-03-",c("06","13","20","26"),sep="")
  epiweek_apr21 <- paste("2021-04-",c("03","10","17","25"),sep="")
  epiweek_may21 <- paste("2021-05-",c("01","08","15","22","29"),sep="")
  epiweek_jun21 <- c(paste("2021-06-",c("05","12","19","26"),sep=""),"2021-07-03")
  
  epiweek_ends <- c(epiweek_jan,epiweek_feb,epiweek_mar,epiweek_apr,epiweek_may,epiweek_jun,
                    epiweek_jul,epiweek_aug,epiweek_sep,epiweek_oct,epiweek_nov,epiweek_dec,
                    epiweek_jan21,epiweek_feb21,epiweek_mar21,epiweek_apr21,epiweek_may21,
                    epiweek_jun21)
  return( epiweek_ends)
}

# uses package lubridate
get_week <- function(dateTxt, useISOWeek=TRUE) {
  if (useISOWeek) {
    return( isoweek(ymd(dateTxt)) )
  } else {
    return( week(ymd(dateTxt)) )
  }
}


# uses package lubridate
get_epiweek <- function(dateTxt) {
  return( epiweek(ymd(dateTxt)) )
}


#this creates an object with the epi week
library(data.table)
epiISL<-as.data.frame(fread("ref_metadata.tsv"))
dateTxt <- as.data.frame(epiISL$'Collection date')
epiweek <- as.numeric(apply(as.data.frame(dateTxt),1,get_epiweek))
epiweek <- as.data.frame(epiweek)
epiweek_under_20 <- filter(epiweek, epiweek < 20)
epiweek_under_20 <- epiweek_under_20 + 52
epiweek_over_20 <- filter(epiweek, epiweek > 20)
epiweek <- rbind(epiweek_over_20,epiweek_under_20)
summary(epiweek)
#this creates a histogram showing the number of sequences over time
library(ggplot2)

ggplot()+geom_histogram(data=epiweek,aes(x=epiweek),bins=17,fill="light blue",col="black")+labs(x="Date of Collection (Epi-Weeks) ",y="Number of Sequences") +
  scale_x_continuous(breaks=seq(0,17,2))

#want to extract the location of each sequence 
library(stringr)
location <- as.data.frame(epiISL$'Location')
location <- data.frame(do.call('rbind',strsplit(as.character(location$`epiISL$Location`),'/',fixed = TRUE)))
g <- ggplot(location) + geom_bar(aes(y = X3), fill="light blue",col="black")
g
#what is the mean number of sequences in each epi-week
epiweek <- as.numeric(apply(as.matrix(dateTxt),1,get_epiweek))
summary(epiweek)
##mean number of sequences = 9.17
#need to calculate number of sequences per-week 
epiweek <- as.data.frame(epiweek)
df <- data.frame(epiweek,epiweek$epiweek)
ew <- df [,'epiweek.epiweek']
ew <- as.factor(df[,'epiweek.epiweek'])
summary(ew)
#need to calculate the number of cases per epiweek 
case <- read.csv('Hong_Kong-COVID-19-data.csv')
#need to change the ordering of the dates
library(dplyr)
library (tidyr)
library(data.table)
library(tidyverse)
b <- data.frame(do.call('rbind',strsplit(as.character(case$Date_reported),'/',fixed = TRUE)))
b$b1 <- paste(b$X3,'-',b$X2)
case$b
df1 <- data.frame(b) %>% unite(b,c('X3','X2','X1'),sep="-",remove =FALSE)
#want to assign epi-week to each case 
dateTxt1 <- as.matrix(df1$b)
epiweek1 <- as.numeric(apply(as.matrix(dateTxt1),1,get_epiweek))
epiweek1 <- as.data.frame(epiweek1)
#want to add epiweek to the cases
epicase <- cbind (epiweek1,case)
#figure out the number of cases per epiweek 
cpw <- aggregate(epicase[,6], list(epicase$epiweek1),sum)
sum (cpw$x)
#total = 1012
#calaulate the % of cases per week 
cpw$percentage <- (cpw[[2]]/1012)
#need to apply the sampling stratgies to the genetic data 
####Method 1  - use the mean number of sequences per week####
#Mean number of sequences per week = 9.17 
#filter all those with sequences above 9
epi <- cbind(epiISL,epiweek)
four <- filter(epi,epiweek == 4)
five <- filter(epi,epiweek == 5)
seven <- filter (epi,epiweek == 7)
twelve <- filter (epi,epiweek ==12)
thirteen <- filter (epi,epiweek == 13)
nine <- rbind (four,five,seven,twelve,thirteen)
#sample all those above 9
setseed = 1000
s1 <- nine %>% group_by(epiweek) %>% sample_n(9, replace = FALSE)
#want to add this to the rest of the data set 
six <- filter(epi,epiweek == 6)
eight <- filter(epi,epiweek == 8)
nin <- filter(epi,epiweek == 9)
ten <- filter(epi,epiweek == 10)
eleven <- filter(epi,epiweek == 11)
fourteen <- filter(epi,epiweek == 14)
fifteen <- filter(epi,epiweek == 15)
sixteen <- filter(epi,epiweek == 16)
uniform <- rbind(s1, six, eight, nin, ten, eleven,fourteen,fifteen,sixteen)
uniform
ggplot()+geom_histogram(data=uniform,aes(x=epiweek),bins=13,fill="light blue",col="black")+labs(x="Date of Collection (Epi-Weeks) ",y="Number of Sequences") +
  scale_x_continuous(breaks=seq(0,16,2))
write.table(uniform, file = "01mar2021_EPI_ISL_uniform_sample_summary.txt",row.names=TRUE,sep=",", quote = TRUE)



####Method 2 - use the number of sequences in proportion to the cases####
cpw$num_seq <-117*cpw$percentage
cpw
#This is too match the number of sequences to each week 
four <- four %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
five <- five %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
six <- six %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
seven <- seven %>% group_by(epiweek) %>% sample_n(3, replace = FALSE)
eight <- eight %>% group_by(epiweek) %>% sample_n(2, replace = FALSE)
nin <- nin %>% group_by(epiweek) %>% sample_n(3, replace = FALSE)
eleven <- eleven %>% group_by(epiweek) %>% sample_n(4, replace = FALSE)
twelve <- twelve %>% group_by(epiweek) %>% sample_n(15, replace = FALSE)
thirteen <- thirteen %>% group_by(epiweek) %>% sample_n(35, replace = FALSE)
fourteen <- fourteen %>% group_by(epiweek) %>% sample_n(32, replace = FALSE)
fifteen <- fifteen %>% group_by(epiweek) %>% sample_n(15, replace = FALSE)
sixteen <- sixteen %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
#10,13,14,15 (the case values for these were larger than the number of sequences available)
#merge these data frames together
by_case_sample <- rbind(four,five,six,seven,eight,nin,ten,eleven,twelve,thirteen,fourteen,fifteen,sixteen)
ggplot()+geom_histogram(data=by_case_sample,aes(x=epiweek),bins=13,fill="light blue",col="black")+labs(x="Date of Collection (Epi-Weeks) ",y="Number of Sequences") +
  scale_x_continuous(breaks=seq(0,16,2))
by_case_sample
write.table(by_case_sample, file = "01mar2021_EPI_ISL_case_sample_summary.txt",row.names=TRUE,sep=",", quote = TRUE)
#Plot graph of num_seq from the case data (in the same proportion as the genetic data)
cpw$round <- round(cpw$num_seq)
barplot(cpw$round)
p<-ggplot(cpw, aes(x=Group.1, y=round)) +
  geom_bar(stat="identity",fill="light blue",col="black") +
 labs(x="Date of Collection (Epi-Weeks) ",y="Number of Sequences")
p
###Method 3: Inverse Proportion of the sequence data
cpw$inverse_sequences <- (1/cpw[[3]])
cpw
#the sequences represent the maximum number of sequences in the weeks
four <- four %>% group_by(epiweek) %>% sample_n(11, replace = FALSE)
five <- five %>% group_by(epiweek) %>% sample_n(13, replace = FALSE)
six <- six %>% group_by(epiweek) %>% sample_n(6, replace = FALSE)
seven <- seven %>% group_by(epiweek) %>% sample_n(20, replace = FALSE)
eight <- eight %>% group_by(epiweek) %>% sample_n(6,replace = FALSE)
nin <- nin %>% group_by(epiweek) %>% sample_n(6, replace = FALSE)
ten <- ten %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
eleven <- eleven %>% group_by(epiweek) %>% sample_n(5, replace = FALSE)
twelve <- twelve %>% group_by(epiweek) %>% sample_n(8, replace = FALSE)
thirteen <- thirteen %>% group_by(epiweek) %>% sample_n(3, replace = FALSE)
fourteen <- fourteen %>% group_by(epiweek) %>% sample_n(4, replace = FALSE)
fifteen <- fifteen %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
sixteen <- sixteen %>% group_by(epiweek) %>% sample_n(1, replace = FALSE)
#merge these data frames together
by_case_sample <- rbind(four,five,six,seven,eight,nin,ten,eleven,twelve,thirteen,fourteen,fifteen,sixteen)
ggplot()+geom_histogram(data=by_case_sample,aes(x=epiweek),bins=13,fill="light blue",col="black")+labs(x="Date of Collection (Epi-Weeks) ",y="Number of Sequences") +
  scale_x_continuous(breaks=seq(0,16,2))
by_case_sample
write.table(by_case_sample, file = "01mar2021_EPI_ISL_case_sample_inverse_summary.txt",row.names=TRUE,sep=",", quote = TRUE)

###Method 4: Flaxman 



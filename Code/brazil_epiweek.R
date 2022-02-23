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
  
  epiweek_starts <- c(epiweek_jan,epiweek_feb,epiweek_mar,epiweek_apr,epiweek_may,epiweek_jun,
                      epiweek_jul,epiweek_aug,epiweek_sep,epiweek_oct,epiweek_nov,epiweek_dec)
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
  
  epiweek_ends <- c(epiweek_jan,epiweek_feb,epiweek_mar,epiweek_apr,epiweek_may,epiweek_jun,
                    epiweek_jul,epiweek_aug,epiweek_sep,epiweek_oct,epiweek_nov,epiweek_dec)
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

library (seqinr)
library (dplyr)
#this creates an object with the epi week
data <- read.fasta("unsampled_Brazil.fasta")
taxa <- as.data.frame(as.matrix(attributes(data)$names))
taxa_split <- data.frame(do.call('rbind',strsplit(as.character(taxa$V1),'|',fixed = TRUE)))
head(taxa_split)
dateTxt <- as.matrix(taxa_split$X6)
epiweek <- as.numeric(apply(as.matrix(dateTxt),1,get_epiweek))
epiweek <- as.data.frame(epiweek)
library(dplyr)
epiweek_big <- filter(epiweek, epiweek > 20)
epiweek_small <- filter(epiweek, epiweek < 20)
epiweek_small <- epiweek_small + 53
epiweek <- rbind (epiweek_big,epiweek_small)
summary (epiweek)
table(epiweek)

#plot proportion of P.1 sequences 

P1 <- read.csv ('P.1_proportion.csv')
P1$Group.1 <- as.Date (c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31',
                 '2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31'))
ggplot(P1, aes(x= Group.1, y = x, fill = Identity
)) + geom_bar(stat = 'identity', position = 'fill') +
  labs(x="Date of Collection",y="Percentage of P.1 (%)") + theme_bw() +
  scale_x_date(date_breaks = "2 week", date_labels = "%Y %b %d") 

#this creates a histogram showing the number of sequences over time

library(ggplot2)
freqs <- aggregate(epiweek$epiweek, by=list(epiweek$epiweek), FUN=length)
freqs$names <- c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31')
unsampled <-ggplot(data=freqs, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="light blue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) +scale_x_discrete(breaks = c('2020-11-29', '2021-01-03','2021-01-31'))
unsampled
#need to calculate number of sequences per-week 
epiweek <- as.data.frame(epiweek)
df <- data.frame(epiweek,epiweek$epiweek)
ew <- df [,'epiweek.epiweek']
ew <- as.factor(df[,'epiweek.epiweek'])
summary(ew)
#need to calculate the number of cases per epiweek 
case <- read.csv('Brazil-COVID-19-data_updated.csv')
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
epiweek_big <- filter(epiweek1, epiweek1 > 20)
epiweek_small <- filter(epiweek1, epiweek1 < 20)
epiweek_small <- epiweek_small + 53
epiweek1 <- rbind (epiweek_big,epiweek_small)
#want to add epiweek to the cases
epicase <- cbind (epiweek1,case)
epicase <- epicase[1:69,]
#figure out the number of cases per epiweek 
cpw <- aggregate(epicase[,5], list(epicase$epiweek1),sum)
deaths <- aggregate(epicase[,7], list(epicase$epiweek1),sum) 
filter <- as.data.frame (cbind(cpw,deaths))
write.csv(filter,'data_for_epifilter.csv')

sum (cpw$x)
#total = 1012
#calaulate the % of cases per week 
cpw$percentage <- (cpw[[2]]/ 8246)
#need to apply the sampling stratgies to the genetic data 
####Method 1  - use the mean number of sequences per week####
#Mean number of sequences per week = 20
#filter all those with sequences above 20
epiISL <- taxa_split$X1
epi <- cbind(epiISL,epiweek)
table (epiweek)
fourty_nine <- filter(epi,epiweek == 49)
fifty <- filter(epi,epiweek == 50)
fifyy_one <- filter(epi,epiweek == 51)
fifyy_two <- filter(epi,epiweek == 52)
fifyy_three <- filter(epi,epiweek == 53)
fifyy_four <- filter(epi,epiweek == 54)
fifyy_five <- filter (epi,epiweek == 55)
fifyy_six <- filter (epi,epiweek ==56)
fifyy_seven <- filter (epi,epiweek == 57)
fifyy_eight <- filter (epi,epiweek == 58)
#sample all those above 9
setseed = 1000
sample <- rbind (fifyy_two,fifyy_three,fifyy_four,fifyy_five,fifyy_six)
s1 <- sample %>% group_by(epiweek) %>% slice_sample(n = 20, replace = FALSE)
#want to add this to the rest of the data set 
uniform <- rbind(s1, fourty_nine, fifty, fifyy_one, fifyy_seven, fifyy_eight)
uniform <- aggregate(uniform$epiweek, by=list(uniform$epiweek), FUN=length)
uniform$names <- c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31')
uniform <-ggplot(data=uniform, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="light blue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) + scale_x_discrete(breaks = c('2020-11-29', '2021-01-03','2021-01-31'))
uniform
#write.table(uniform, file = "01mar2021_EPI_ISL_uniform_sample_summary.txt",row.names=TRUE,sep=",", quote = TRUE)
####Method 2 - use the number of sequences in proportion to the cases####
cpw$num_seq <-196*cpw$percentage
cpw
#This is too match the number of sequences to each week 
fourty_nine <- fourty_nine %>% group_by(epiweek) %>% slice_sample(n = 1, replace = FALSE)
fifty <- fifty %>% group_by(epiweek) %>% slice_sample(n = 2, replace = FALSE)
fifyy_one <- fifyy_one %>% group_by(epiweek) %>% slice_sample(n = 6, replace = FALSE)
fifyy_two <- fifyy_two %>% group_by(epiweek) %>% slice_sample(n = 15, replace = FALSE)
fifyy_three <- fifyy_three %>% group_by(epiweek) %>% slice_sample(n = 24, replace = FALSE)
fifyy_four <- fifyy_four %>% group_by(epiweek) %>% slice_sample(n = 32, replace = FALSE)
fifyy_five <- fifyy_five %>% group_by(epiweek) %>% slice_sample(n = 38, replace = FALSE)
fifyy_six <- fifyy_six %>% group_by(epiweek) %>% slice_sample(n = 30, replace = FALSE)
fifyy_seven <- fifyy_seven %>% group_by(epiweek) %>% slice_sample(n = 22, replace = FALSE)
fifyy_eight <- fifyy_eight %>% group_by(epiweek) %>% slice_sample(n = 25, replace = FALSE)
#merge these data frames together
by_case_sample <- rbind(fourty_nine,fifty,fifyy_one,fifyy_two,fifyy_three,fifyy_four,fifyy_five,fifyy_six,fifyy_seven,fifyy_eight)
by_case_sample <- aggregate(by_case_sample$epiweek, by=list(by_case_sample$epiweek), FUN=length)
by_case_sample$names <- c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31')
proportional <-ggplot(data=by_case_sample, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="light blue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) +scale_x_discrete(breaks = c('2020-11-29', '2021-01-03','2021-01-31'))
proportional
#write.table(by_case_sample, file = "01mar2021_EPI_ISL_case_sample_summary.txt",row.names=TRUE,sep=",", quote = TRUE)
###Method 3: Inverse Proportion of the sequence data
cpw$inverse_sequences <- (1/cpw[[3]])
cpw
#the sequences represent the maximum number of sequences in the weeks
fourty_nine <- fourty_nine %>% group_by(epiweek) %>% slice_sample(n = 145, replace = FALSE)
fifty <- fifty %>% group_by(epiweek) %>% slice_sample(n = 109, replace = FALSE)
fifyy_one <- fifyy_one %>% group_by(epiweek) %>% slice_sample(n = 31, replace = FALSE)
fifyy_two <- fifyy_two %>% group_by(epiweek) %>% slice_sample(n = 12, replace = FALSE)
fifyy_three <- fifyy_three %>% group_by(epiweek) %>% slice_sample(n = 8, replace = FALSE)
fifyy_four <- fifyy_four %>% group_by(epiweek) %>% slice_sample(n = 6, replace = FALSE)
fifyy_five <- fifyy_five %>% group_by(epiweek) %>% slice_sample(n = 5, replace = FALSE)
fifyy_six <- fifyy_six %>% group_by(epiweek) %>% slice_sample(n = 7, replace = FALSE)
fifyy_seven <- fifyy_seven %>% group_by(epiweek) %>% slice_sample(n = 9, replace = FALSE)
fifyy_eight <- fifyy_eight %>% group_by(epiweek) %>% slice_sample(n = 8, replace = FALSE)
#merge these data frames together
inverse <- rbind(fourty_nine,fifty,fifyy_one,fifyy_two,fifyy_three,fifyy_four,fifyy_five,fifyy_six,fifyy_seven,fifyy_eight)
inverse <- aggregate(inverse$epiweek, by=list(inverse$epiweek), FUN=length)
inverse$names <- c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31')
inverse <-ggplot(data=inverse, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="light blue",col="black")+labs(x="Date of Collection",y="Number of Sequences") +
  theme_bw() + ylim (0,40) +scale_x_discrete(breaks = c('2020-11-29', '2021-01-03','2021-01-31'))
inverse
#write.table(inverse, file = "01mar2021_EPI_ISL_case_sample_inverse_summary.txt",row.names=TRUE,sep=",", quote = TRUE)
#merge all graphs into a figure 
cpw$names <- c('2020-11-29','2020-12-06','2020-12-13','2020-12-20','2020-12-27','2021-01-03','2021-01-10','2021-01-17','2021-01-24','2021-01-31')
case <- ggplot(data=cpw, aes(x=names, y=x)) +
  geom_bar(stat="identity", fill="light blue",col="black")+labs(x="Date of Collection",y="Number of Cases") +
  theme_bw()+scale_x_discrete(breaks = c('2020-11-29', '2021-01-03','2021-01-31'))
library(ggpubr)
ggarrange(case,                                                 # First row with scatter plot
          ggarrange(unsampled, proportional, inverse, uniform, labels = c("B", "C",'D','E')),
          labels = "A"                                        # Labels of the scatter plot
) 

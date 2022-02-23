#read the tsv file
library(data.table)
data<-as.data.frame(fread("gisaid_hcov-19_2020_12_22_18.tsv"))
#this is to generate a txt file you want to have seporators 
library(dplyr)
library (tidyr)
df1 <- data%>% unite (hk,'Accession ID': 'Virus name',sep=",")
df2 <- data %>% unite (hk2,'Accession ID': 'Collection date',sep=",")
df3 <- df2 %>% unite (hk3,'hk2': 'Location',sep=",")
df5 <- bind_cols (df1$hk,df3$hk3)
df6 <- df5 %>% unite (hk,'...1':'...2' ,sep="|")
colnames(df6) <- "EpiID,SeqName,Date,Location"
df7 <- df6[,1, drop=FALSE]
#to export as text
write.table(df7, file = "22dec2020_EPI_ISL_summary.txt",row.names=FALSE,sep="\t", quote = FALSE)

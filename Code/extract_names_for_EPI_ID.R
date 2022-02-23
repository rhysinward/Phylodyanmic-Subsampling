# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#set WD
setwd("/Users/rhysinward/Documents/Master's_paper/Data")
# Folder path for results
folres = (path = "./Results/")
# Main functions to run script
files.sources = list.files(path = "./Main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}



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
library(ape)
library(plyr)

#load trees 

data_unsampled <- as.data.frame(fread("Brazil/BDSKY/Unsampled/BDSKY_Trail_2/unsampled_Brazil.trees",sep = "",header = FALSE))

#extract ID for tree

tree_ID_unsampled <-  data.frame(do.call('rbind',strsplit(as.character(data_unsampled$V1),'|',fixed = TRUE)))

tree_ID_unsampled <- tree_ID_unsampled %>% slice(6:201)

tree_ID_unsampled$X8 <- 'EPI_ISL_' 

tree_ID_unsampled$X9 <- paste(tree_ID_unsampled$X8, tree_ID_unsampled$X1, sep = "" )

tree_ID_unsampled_AI <- as.data.frame(tree_ID_unsampled$X9)

write.csv(tree_ID_unsampled_AI, file = "519_treemmer_sequences.csv",row.names=TRUE, quote = TRUE)

#Proportional

#load trees 

data_proportional <- as.data.frame(fread("Brazil/BDSKY/Proportional/trail_2/proportional_Brazil.trees",sep = "",header = FALSE))

#extract ID for tree

data_proportional <-  data.frame(do.call('rbind',strsplit(as.character(data_proportional$V1),'|',fixed = TRUE)))

data_proportional <- data_proportional %>% slice(5:174)

data_proportional$X8 <- 'EPI_ISL'

data_proportional$X9 <- paste(data_proportional$X8, data_proportional$X1, sep = "_" )

tree_ID_proportional_AI <- as.data.frame(data_proportional$X9)

write.csv(tree_ID_proportional_AI, file = "519_treemmer_sequences.csv",row.names=TRUE, quote = TRUE)

#uniform

#load trees 

data_uniform <- as.data.frame(fread("Brazil/BDSKY/Uniform/trail_2/uniform_Brazil.trees",sep = "",header = FALSE))

#extract ID for tree

data_uniform <-  data.frame(do.call('rbind',strsplit(as.character(data_uniform$V1),'|',fixed = TRUE)))

data_uniform <- data_uniform %>% slice(5:155)

data_uniform$X8 <- 'EPI_ISL'

data_uniform$X9 <- paste(data_uniform$X8, data_uniform$X1, sep = "_" )

tree_ID_uniform_AI <- as.data.frame(data_uniform$X9)

write.csv(tree_ID_uniform_AI, file = "519_treemmer_sequences.csv",row.names=TRUE, quote = TRUE)


#load trees 

data_inverse <- as.data.frame(fread("Brazil/BDSKY/Inverse/trail_2/inverse_Brazil.trees",sep = "",header = FALSE))

#extract ID for tree

tree_ID_inverse <-  data.frame(do.call('rbind',strsplit(as.character(data_inverse$V1),'|',fixed = TRUE)))

tree_ID_inverse <- tree_ID_inverse %>% slice(5:72)

tree_ID_inverse$X8 <- 'EPI_ISL'

tree_ID_inverse$X9 <- paste(tree_ID_inverse$X8, tree_ID_inverse$X1, sep = "_" )

tree_ID_inverse_AI <- as.data.frame(tree_ID_inverse$X9)

tree_ID_inverse_AI$`tree_ID_inverse$X9`

write.csv(tree_ID_inverse_AI, file = "519_treemmer_sequences.csv",row.names=TRUE, quote = TRUE)







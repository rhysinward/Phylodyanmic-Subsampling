library (ape)
library(dplyr)
getEl	<- function( line, sep=",", ind=-1, final=FALSE, reconstruct=FALSE, ex=-1, fromEnd=FALSE ) {
  els	<- strsplit(line, sep)[[1]]
  
  if (ind[1] != -1) {
    if (fromEnd) {
      ind <- length(els)-(ind-1)
    }
  }
  
  if (final) {
    return( els[length(els)] )
  } else {
    
    if (reconstruct) {
      if (ex[1] > 0) {
        if (fromEnd) {
          ex <- length(els)-(ex-1)
        }
        ind <- setdiff((1:length(els)),ex)
      }
      
      newLine <- els[ind[1]]
      if (length(ind) > 1) {
        for (i in 2:length(ind)) {
          newLine <- paste(newLine, els[ind[i]], sep=sep)
        }
      }
      return ( newLine )
    } else {
      if ( ind[1] == -1 ) {
        return( els )
      } else {
        return( els[ind] )
      }
    }
  }
}
path = ''
info <- read.csv(paste(path,"sequences_to_remove_locations_part_2.csv",sep=""))
print(info)
combinedFname <- paste("Brazil_locations_post_tempest_450_seq.fasta",sep="")
seqs <- read.dna(paste(path,combinedFname,sep=""),format="fasta", as.matrix=FALSE)
taxa <- as.data.frame(as.matrix(attributes(seqs)$names))
isName <- apply(taxa, 1, getEl, ind=1, sep="\\|")
epiISL <- apply(taxa, 1, getEl, ind=2, sep="\\|")
taxa$minds  <- as.numeric(match(isName, info$EpiID))
dim(taxa)
taxa$minds <- as.numeric(taxa$minds)
taxa <- filter(taxa, minds > 0)
#load package
library(seqinr)
#load file containing sequences
data<-read.fasta("Brazil_locations_post_tempest_450_seq.fasta")
#create a vector containing species names: these are the last part of the string
vec.names<-unlist(lapply(strsplit(names(data), ";"), function(x)x[length(x)]))
#find names to keep: indices which are not in the species to remove
species.to.remove<-taxa$V1
vec.tokeep<-which(! vec.names %in%  species.to.remove)
length(vec.tokeep)
#write the final output
write.fasta(sequences=data[vec.tokeep], names=names(data)[vec.tokeep], file.out="Brazil_locations_post_tempest_419_seq.fasta")

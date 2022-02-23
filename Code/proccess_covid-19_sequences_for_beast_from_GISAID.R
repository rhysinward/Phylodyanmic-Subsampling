# script to rename sars-cov-2
#load ape library 
library (ape)
#Split sequence name 
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
#calulate decimal date 
calcDecimalDate	<- function(day, month, year, defaultMonth=6, defaultDay=15) {
  cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334)
  
  if (month==0) {
    if (defaultMonth >= 1) {
      month <- defaultMonth
    } else {
      month	<- ceiling(runif(1)*12)
    }
  }
  
  if (day==0) {
    if (defaultDay >= 1) {
      day	<- defaultDay
    } else {
      day	<- ceiling(runif(1)*30)
    }
  }
  
  dd	<- cd[month] + day - 1
  
  decDate <- year + (dd/365)
  
  return ( decDate )
}


invertDecimalDate <- function( decDate, formatAsTxt=FALSE, ddmmyy=FALSE ) {
  cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
  fractD<- cd/365
  
  year		<- floor(decDate)
  fractYear 	<- decDate-year
  month		<- which(fractD >= fractYear)[1]-1
  
  if (month > 0) {
    fractMonth  <- fractYear-fractD[month]
    day		<- round((fractMonth*365)+1)
  } else {
    month <- 1
    day   <- 1
  }
  
  res <- c(day,month,year)
  
  if (formatAsTxt) {
    if (month < 10) {
      mm  <- paste("0",month,sep="")
    } else {
      mm <- month
    }
    res <- paste(year,mm,day,sep="-")
  }
  
  if (ddmmyy) {
    months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    if (day < 10) {
      dd <- paste("0",day,sep="")
    } else {
      dd <- day
    }
    res <- paste(dd,months[month],year,sep="-")
  }
  return( res )
  
}



calcDecimalDate_fromTxt	<- function( dateTxt, sep="/", namedMonths=FALSE, dayFirst=FALSE) {
  els 	<- strsplit(dateTxt, sep)[[1]]
  if (dayFirst) {
    if (length(els) > 1) {
      els <- els[length(els):1]
    }
  }
  
  year 	<- as.integer(els[1])
  
  if (length(els)==1) {
    month <- 6  #7
    day	<- 15 #2
    decDate <- year + 0.5
  } else {
    
    if (length(els)==2) {
      if (nchar(els[2]) > 0) {
        if (namedMonths) {
          month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
        } else {
          month <- as.integer(els[2])
        }
        day	<- 15
        decDate <- calcDecimalDate(day,month,year)
      } else {
        month <- 6 #7
        day   <- 15 #2
        decDate <- year + 0.5
      }
    } else {
      if (namedMonths) {
        month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
      } else {
        month <- as.integer(els[2])
      }
      
      if (nchar(els[3]) > 0) {
        day 	<- as.integer(els[3])
      } else {
        day <- 15
      }
      decDate <- calcDecimalDate(day,month,year)
    }
  }
  
  
  return ( decDate )
}


calcDecimalDate_from_yymmdd	<- function( dateTxt, sep="/", ycutoff=15, defaultMonth=6, defaultDay=15 ) {
  els	<- strsplit(dateTxt, sep)[[1]]
  yy	<- as.integer(els[1])
  mm	<- as.integer(els[2])
  dd	<- as.integer(els[3])
  
  if (!is.finite(yy)) {
    return( -1 )
  } else {
    if (yy <= ycutoff) {
      yy <- yy+2000
    }
    if ((yy > ycutoff) & (yy < 99)) {
      yy <- yy+1900
    }
    
    if (!is.finite(mm)) {
      mm <- 0
    }
    if (!is.finite(dd)) {
      dd <- 0
    }
    return ( calcDecimalDate( dd, mm, yy, defaultMonth=defaultMonth, defaultDay=defaultDay ) )
  }
  
}
#make sure to point to the correct directory 
path <- ""
fnames<- dir()
fnames<- fnames[grep(".fasta",fnames)]
#remember to rename the file with all the sequences BetaCoV_Wuhan_177seqs_22dec2020 
# create one file with all of the sequences
#need to rename file
combinedFname <- paste("BetaCoV_Hong_Kong_",177,"seqs_22dec2020.fasta",sep="")
pos <- grep(combinedFname,dir(),fixed=TRUE)
if (length(pos)==0) {
  for (i in 1:length(fnames)) {
    #seqs <- read.dna( paste(path,fnames[i],sep=""), format="fasta", as.matrix=FALSE)
    #write.dna(seqs,file=paste(path,combinedFname,sep=""),format="fasta",nbcol=-1,colsep="",append=(i>0))
    temp <- readLines(paste(path,fnames[i],sep=""))
    write(temp,file=paste(path,combinedFname,sep=""),append=(i>0))
  }
  print(paste("Written",combinedFname))	
} else {
  print(paste("Already done",combinedFname))
}

#read the sequence information file
info <- read.csv(paste(path,"11juneinverseproportion2020_EPI_ISL_summary.txt",sep=""))
print(info)
#read file and alter the sequence names 
rootname <-	paste(gsub("\\.fas","",combinedFname),"_beastNames",sep="")
pos <- grep(paste(rootname,".fas",sep=""),dir(),fixed=TRUE)
if (length(pos)==0) {
  
  seqs <- read.dna(paste(path,combinedFname,sep=""),format="fasta", as.matrix=FALSE)
  taxa <- as.matrix(attributes(seqs)$names)
  isName <- apply(taxa, 1, getEl, ind=1, sep="\\|")
  epiISL <- apply(taxa, 1, getEl, ind=2, sep="\\|")
  minds  <- match(epiISL, info$EpiID)
  all( epiISL==info$EpiID[minds] )
  dateTxt <- as.matrix(info$Date[minds])
  decDate <- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
  location<- as.matrix(info$Location[minds])
  country <- apply(as.matrix(location), 1, getEl, ind=1, sep=" / ")
  state   <- apply(as.matrix(location), 1, getEl, ind=2, sep=" / ")
  place   <- apply(as.matrix(location), 1, getEl, ind=3, sep=" / ")
  
  pos     <- which((state=="Kanagawa Prefecture") & is.na(place))
  place[pos] <- "Yokohama"
  print(paste("Changed place of",info$SeqName[pos],"from NA to Yokohama (capital city)"))
  
  newTaxa <- paste(gsub("EPI_ISL_","",epiISL),country,state,place,isName,dateTxt,format(decDate,digits=7),sep="|")
  newTaxa <- gsub(" ","_",newTaxa)
  attributes(seqs)$names <- newTaxa
  write.dna(seqs, file=paste(path,rootname,".fas",sep=""), format="fasta", nbcol=-1, colsep="")
  
  newInfo <- cbind(newTaxa,epiISL,isName,dateTxt,country,state,place,decDate)
  colnames(newInfo) <- c("SeqName","EPI_ISL","IsolateName","CollectionDate","Country","State","Place","decDate")
  write.table(newInfo,file=paste(path,rootname,"_infoTbl.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
} else {
  print("Already done renaming")
}		


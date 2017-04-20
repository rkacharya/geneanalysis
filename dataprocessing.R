#Data Tranformation Functions
#Written by Rajesh Acharya
#Description: sets of functions that can be used to simplify expression discovery
#Uncomment below to install if necessary
##source('http://bioconductor.org/biocLite.R')
##biocLite('preprocessCore')
#load package
library(preprocessCore)
library(plotly)

####DATA TRANSFORMATION FUNCTIONS####
#Quantile normalize data (from RAW)
#var:data = rawdata. var:omitcols = columns to omit in calculations (nonpatientcol)
selfquantnorm <- function(data, omitcols){
  quantnormdata <- data.frame(normalize.quantiles(as.matrix(data[,!omitcols])))
  colnames(quantnormdata) <- colnames(data[!omitcols]) #add column identifiers back in
  lognormeddata <- cbind(data[,omitcols],log10(quantnormdata)) #take log10 and merge gene identifiers
  return(lognormeddata)
}
#log10 of Quantile Normalized data. Shifted +1 to account for 0.
log10datashift <- function(data2, omitcols){
  logdata <- cbind(data2[,omitcols],log10(data2[,!omitcols]+1)) #shifted +1 to deal with 0 counts
  return(logdata)
}
#log10 of Quantile Normalized data. 


#Z-score data (RAW,log or normalized)
z_scored_data <- function(data1, omitcols){
  meanexpression <- rowMeans(data1[ ,!omitcols])
  sdexpression <- c()
  for (x in rownames(data1)){
    sdexpression[x] <- sd(data1[x,!omitcols]) #Change: method should take sample data also
  }
  
  zdata <- cbind(data1[,omitcols]
                 ,(data1[rownames(data1),!omitcols] 
                   - meanexpression[names(meanexpression)])
                 /sdexpression[names(sdexpression)])
  return(zdata)
}



####Readability Functions####
#Function to split multiple values in a cell by its separator
#e.g. cellsplit("10/24/2016","/") > c(10,24,2016) as string
cellsplit <- function(cellvalue,sep){
  unlist(strsplit(as.character(cellvalue), sep))
}

#Substring n characters from the end of string "x"
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Import patient response data
patientfileread <- function(file){
  read.table(file,
             sep="\t",
             header=TRUE)
}

#difference in date

diff_in_days <- function (strtdt,progdt,dod,lastdt){
  surv <- c()
  alv0dead1prog2 <- c()
  for (x in 1:length(strtdt)){
    if (is.na(as.Date(progdt[x],'%m/%d/%Y'))){
      if (is.na(as.Date(dod[x],'%m/%d/%Y'))){
        surv[x] <- as.numeric(difftime(as.Date(lastdt[x],'%m/%d/%Y'),as.Date(strtdt[x],'%m/%d/%Y'),units="days"))
        alv0dead1prog2[x] <- as.numeric(0)
      }else{
        surv[x] <- as.numeric(difftime(as.Date(dod[x],'%m/%d/%Y'),as.Date(strtdt[x],'%m/%d/%Y'),units="days"))
        alv0dead1prog2[x] <- as.numeric(1)
      }
    }else{
      surv[x] <- as.numeric(difftime(as.Date(progdt[x],'%m/%d/%Y'),as.Date(strtdt[x],'%m/%d/%Y'),units="days"))
      alv0dead1prog2[x] <- as.numeric(2)
    }
  }
  cbind(surv,alv0dead1prog2)
}

####PLOT FORMATING FUNCTIONS####

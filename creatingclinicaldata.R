#Response Data identifier merge
library(dplyr)
library(tidyr)
#Map filenames
responsefile <- "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\Immunotherapy Patient List w- Project Tabs and RESPONSE v5 Edited.txt"
samplelistfile <- "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\Sample list.txt"
tumorpurityfile <- "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\Tumor purity 20170127.txt"
savefile <- "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\Immunotherapy Patient List DeIdentified RA 20170127.txt"



#Read in tab delimited files to dataframes
responsedata <- read.table(responsefile,
                           sep="\t",
                           header=TRUE
)

samplelistdata <- read.table(samplelistfile, 
                             sep="\t", 
                             header=TRUE
)

tumorpuritydata <- read.table(tumorpurityfile,
                              sep="\t",
                              header=TRUE
                              )

responsedata2 <- responsedata %>%
  mutate(Sample = as.character(Sample)) %>%
  mutate(Sample = replace(Sample,which(Sample==""),"N/A")) %>%
  mutate(Sample = replace(Sample,which(is.na(Sample)),"NA")) %>%
  mutate(Sample = strsplit(as.character(Sample),",")) %>%
  unnest(.drop=FALSE) %>%
  select(Name,Patient.ID,Sample,everything()) %>%
  mutate(Sample = replace(Sample,which(Sample=="N/A"),NA)) %>%
  mutate(Sample = as.numeric(Sample)) %>%
  mutate(Tumor.Purity = tumorpuritydata$RNA_Per[match(Sample, tumorpuritydata$SampleNum)])

# responsedata2 <- responsedata2 %>%
#   mutate(Tumor.Purity = tumorpuritydata$RNA_Per[match(Sample, tumorpuritydata$SampleNum)])


deidentified <- responsedata2[,!(colnames(responsedata2) %in% c("Name", "MRN"))]
  

# samplesplit <- strsplit(as.character(responsedata$Sample),",")
#   
# expandedvec <- samplesplit %>%
#   unlist %>%
#   as.numeric
# 
# sapply(samplesplit,length) %>%
  


write.table(deidentified, savefile, sep="\t", na = "",row.names = FALSE)



#match the contents of Name in responsedata with all similar data in samplelist$Name. Prob with grep
#Save row matches as MatchVector variable
#for the column samplelistdata$sample#forjounce line up with MatchVector and extract all values into Samples_vec
#Use string merging to turn all values in the Samples_vec into a string_of_numbers separated by commas
#Place string_of_numbers into responsedata$Sample at row of currently being matched

#Create list of samples to use associated with names. match and place in Responsedata$Sample.to.use

#Create patient identifiers




##nanostring data filter
##Written By Rajesh Acharya
##
##
##To Implement: 1. Selection of target variables/range. Bind publishable names to genes
##2.Read CSV and format based on headings and selections. Use $ methods
##3.Call gplots and draw heatmap
##
##
##Future: Implement user friendly data selector, search algorithm for genes of interest,
##Directory selection, 


this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#setwd("C:/Users/z33ky/Documents/UC stuff") #Set Working Directory
filename = "CHI02MesoNanostringData.csv"
patientrow = 2 #row where the names of patients are located
datacol = 4 #column where patient data first appears

#Set up data
colsize <- ncol(read.csv(filename, #Gather # of columns for future calc
                       skip = patientrow-1, 
                       header = TRUE, 
                       nrows = 1
                       )
                   )
rawdata <- read.csv(filename,  #Import bulk of raw data & format
                    skip = patientrow-1,
                    header = TRUE, 
                    sep = ",",
                    row.names = NULL,
                    colClasses = c("factor", #SortableData
                                   rep("character",2), #Character Data
                                   rep("numeric",colsize-3) #rest of data is numeric
                                   )
                    
                    )

Response <- as.character(read.csv(filename, #Get Response attribute seperately per patient
                     sep = ",",
                     nrows = 1,
                     header = FALSE,
                     colClasses = "character",
                     fill = TRUE
                        )
                     )

#Set Genes of interest
#TODO: set up function to search Gene column and add to GEP
GEP <- c("TIGIT","CD27","CD8A","PDCD1LG2","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
name_vector <- c("TIGIT","CD27","CD8A","PDCD1LG2(PD-L2)","LAG3","CD274(PD-L1)","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA.DQA1","CD276","STAT1","HLA.DRB1","HLA.E")
match_list <- match(GEP, rawdata$Name)
sortdata <- rawdata[match_list, ]
#sortdata.response <- Response[4:colsize]

num_data <- sortdata[,4:colsize]
rownames(num_data) <- name_vector
final_data <- data.matrix(num_data)
library("gplots")
heatmap.2(final_data, col=bluered(100), key=TRUE, scale="row", trace="none")
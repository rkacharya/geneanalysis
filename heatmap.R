library("gplots")
#setwd// set working directory to same as where script is
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#Read in CSV. Row 1 & 2 are headers. Column 1-3 are headers.
rawdata <- read.csv(file="CHI02MesoNanostringData.csv", header=FALSE, sep="," )
#Reformat Raw Data

#labeling of categories
columnnamesraw <- rawdata[2,4:as.numeric(length(rawdata))]
columnnames<- apply(as.vector(columnnamesraw),2,as.character) #had to use conversion for format

#Specify/Extract data for genes of interest
#TIGIT, CD27, CD8A, PDCD1LG2(PD-L2), LAG3, CD274 (PD-L1), CXCR6, CMKLR1, NKG7, CCL5, PSMB10, IDO1, CXCL9, HLA.DQA1, CD276, STAT1, HLA.DRB1, HLA.E
GEP <- c("TIGIT","CD27","CD8A","PDCD1LG2","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
name_vector <- c("TIGIT","CD27","CD8A","PDCD1LG2(PD-L2)","LAG3","CD274(PD-L1)","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA.DQA1","CD276","STAT1","HLA.DRB1","HLA.E")

#Prepopulate variables and get rows of interest
rowofinterest <- 0
for (i in 1:as.numeric(length(GEP))) {
  rowofinterest[i] <- as.numeric(grep(GEP[i], rawdata[ ,2]))
}

filtered_data <- rawdata[rowofinterest,4:as.numeric(length(rawdata)) ] #created GEP filtered data
num_data <- apply(as.matrix(filtered_data),2,as.numeric) #conversion to numeric for heatmap

rownames(num_data) <- name_vector #labeling rows
colnames(num_data) <- columnnames #labeling columns

#Sort first by response then pearson correlation co-efficient

#



#heatmap.2
heatmap.2(num_data, col=bluered(100), key=TRUE, scale="row", trace="none")

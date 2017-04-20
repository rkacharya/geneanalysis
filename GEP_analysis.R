#GEP Analysis and Plot Script

#Data preparation
#Module: Use GEP algorithm to analyze all genes and export P value from student t-test for each gene in GEP
#Module: Use data formatter (nanostringdataimport.R) to construct working matrix on gene signature
#Module: Suplement data formatter to append accessory data to subjects(response, survival)
#Output: GEP matrix with normalized data, supplemental data per subject, P value from GEP calc.
#Output: Classification for control groups


#Ideas: User input for perform calculations or select genes manually

##########################
##Nanostring Data Import##
##########################
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#setwd("C:/Users/z33ky/Documents/UC stuff") #Set Working Directory
filename = file.choose(new=FALSE)
patientrow = 2 #row where the names of patients are located
datacol = 4 #column where patient data first appears

rawdata <- read.table(filename,
                      sep="\t",
                      header=TRUE
                      )
nonpatientcol <- grepl('Code|Name|Accession', names(rawdata))
#Genes of interest
geneset <- c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")
#subset raw data using geneset
genesetdata <- rawdata[match(geneset, rawdata$Name), ]

#response <- #get data from new spreadsheet
meanexpression <- rowMeans(genesetdata[ ,!nonpatientcol])

sdexpression <- c()
for (x in rownames(genesetdata)){
  sdexpression[x] <- sd(genesetdata[x,!nonpatientcol]) #Change: method should take sample data also
}

#calculate z scores
zscoreddata <- cbind(genesetdata[,nonpatientcol]
                     ,(genesetdata[rownames(genesetdata),!nonpatientcol] 
                      - meanexpression[names(meanexpression)])
                     /sdexpression[names(sdexpression)])

sigscore <- colMeans(zscoreddata[!nonpatientcol])
additivescore <- colSums((zscoreddata[!nonpatientcol]))






#Test to see if matrix can be truncated by excluding column names
truncate_test <- rawdata[ ,!grepl('Code|Name|Accession', names(rawdata))]

#Module: Scoring
#Module performs scoring on data. Input(GEP,supp data,P-value)
#
# Function: z-score
#   for x = (1:#subjects)
#     mean_exp <- mean(GEP')
#     sd_exp <- sd(GEP')
#     ztotal <- 0
#     for y = (1:#genes)
#       zscore <- ((selected gene value) - mean_exp)./sd_exp
#       ztotal =+ zscore
#     return (ztotal)
#     score(x) = ztotal./(#Genes)
#   return (score) ##Output vector size of #subjects
#
# Function: pearson-scoring
#   incomplete: will obtain control values and use pearson correlation to control to score.
# 
# Function: signature-score (S-score)
#   Relevant paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3725832/
#   input: control_P-value_vector ## P values of genes in control experiment used to obtain GEP
#   input: control_expr_vector ## expr values of genes from control experiment
#   
#   
#   for x = (1:#subjects)
#     mean_exp <- mean(GEP') #mean and sd calculated across all samples
#     if test sample profiled on same platform as case-control
#       sd_exp <- (((#testsamples-1)*sd(gene,testsamples)^2 + (#controlsamples-1)*sd(gene,testsamples)^2)/(#testsamples+#controlsamples-2))^(.5)
#     else
#       sd_exp <- sd(GEP')
#     #**if fold vector must be reset before each iteration add here**
#     fold <- GEP(:,x)./(control_exp_vector)
#     fold_sign <- sign(fold)
#     ztotal <- 0
#     pstar <- 0
#     sum_pstar <- 0
#     sscore_total <- 0
#     for y = (1:#genes)
#       fold(y) <- GEP(y,x)/control_exp_vector(y) ##no need to store fold long term after troubleshooting
#       fold_sign(y) <- sign(fold)
#       zscore <- (GEP(y,x) - mean_exp)./sd_exp
#       if -log_10(control_P-value_vector(y)) < 1
#         Pstar <- 0
#       elif -log_10(control_P-value_vector(y)) > 1 && -log_10(control_P-value_vector(y)) <= 4
#         pstar <- -log_10(control_P-value_vector(y)) - 1
#       elif -log_10(control_P-value_vector(y)) > 4
#         pstar <- 3
#       return pstar
#       sum_pstar =+ abs(pstar)
#       sscore <- (fold_sign(y))*(pstar)*zscore
#       sscore_total =+ sscore
#     score(x) <- sscore_total./sum_pstar
#   return (score)
#
#
#Module: Plotting
#Module plots the scored data with labels
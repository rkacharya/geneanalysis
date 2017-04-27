#HNC Acquired Resistance script
#Created by Rajesh Acharya on Apr 20th 2017
#Desciption: Use expression data to find differences between acquired resistance patientsk

#Load Packages
#source("https://bioconductor.org/biocLite.R") #biocondoctor source. use biocLite("pckgname") to install
library(dplyr)
library(limma)
library(edgeR)
library(DESeq)
library(edge)
library(ggplot2)


####SEMI-STATIC VARIABLES####

#Expression data path. Use "file.choose(new=FALSE)" for file browser
##filename = file.choose(new=FALSE)
exp_file = "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\HNSCC\\CHI-02.HNSCC.noOutliers.NanostringData.QuantileNormalized.txt"
#File save path. Leave as "" for file browser or declare path.
filesave = ""
#Clinical data file
clin_file = "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\Immunotherapy Patient List DeIdentified RA 20170127.txt"
ARsamp = list(AA=c(1), HE=c(4,81), MH=c(80),PG=c(51),TS=c(48,79)) #Hard list of acquired resistance samples to compare
BLsamp = list(AA=c(9), HE=c(78), MH=c(7),PG=c(77),TS=c(47)) #Hard list of corresponding baseline samples
Nsamp = list(DF=c(84))
#--------------------------------------------------------------------------------------------------#

####DATA FORMAT & CONDITIONING####

expdata <- read.table(exp_file,
                      sep="\t",
                      header = TRUE)
clindata <- read.table(clin_file,
                       sep="\t",
                       header = TRUE)


#exclude columns that have non sample data
npc <- grepl('Code|Name|Accession', names(expdata))




####FUNCTIONS####

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

`%!in%` <- Negate(`%in%`)


####Limma####


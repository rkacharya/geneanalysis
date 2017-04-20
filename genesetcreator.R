#Gene Set Creator
#Rajesh Acharya
#Created 2017 February 02
#Description: Compile and process Meso data. Looking for genes of interest.
library(dplyr)
library(ggplot2)
library(limma)
library(edgeR)
library(qvalue)

####NOTES & TODO####
#Import new clinical data and relate old functions to it
#Find new categorizing metric
#Learn significance of plots corelating genes to genes
#Learn known pathways in Meso. See if signatures needed for expression rather than gene
#Learn known biology in meso
#find proper way to apply DGE to datasets.







####USER VARIABLES####
dfg <- "C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\Meso - Jan2017\\CHI-02.Meso.NanostringData.Normalized.txt"
#dfg <- file.choose(new=FALSE) #Comment if hardcoding data in
dfsd <- "C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\Survival Data.txt"
dfpsd <- "C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\patient sample data.txt"
dfrdr <- "C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\response data raw.txt"


####GENE SIGNATURES####

IFNG.6.Gene <-  list(genes = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"), Name = "IFN-g 6-gene Signature")




####DATA IMPORT & SELECTION####

rd <- read.table(dfg,
                  sep="\t",
                  header=TRUE)
sd <- read.table(dfsd,
                 sep="\t",
                 header=TRUE)
psd <- read.table(dfpsd,
                  sep="\t",
                  header=TRUE)
rdr <- read.table(dfrdr,
                  sep="\t",
                  header=TRUE)
#Clinical info that Tanguy sent. Imported to view but not use.
mcd <- read.table("C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\Mainclinicaldata.txt",
                  sep="\t",
                  header=TRUE)
colnames(mcd)[1] <- "Patient.ID"

compiled <- psd %>%
  cbind(sd[match(psd$Study,sd$Sample),colnames(sd)]) %>%
  cbind(rdr[match(psd$Study,rdr$Patient.ID),1:length(colnames(rdr))])

compiled[,18] <- NULL #Removes second name vector. Uncomment to check
compiled$Best.Patient.Response[grep("SD ",compiled[,"Best.Patient.Response"])] <- "SD"
compiled$Best.Patient.Response <- factor(compiled$Best.Patient.Response, levels=c("SD","PR","PD","Early death"))
#compiled <- select(compiled,Name,Screen,MRN,Study,respond.status.1.20.17,responses.expanded,Best.Patient.Response,PFS,Study.Time,PFS.details,MPS,mps)
  
####FUNCTIONS####
npc <- function(dta){
  grepl('Code|Name|Accession', names(rd))
}

#Edit colMeans section to include
gene_sig <- function(dfr,setid,gv,np){
  if (!all(gv %in% dfr[,"Name"])){
    print(gv[!(gv %in% dfr[,"Name"])], " does not exist in Names for ", setid)
  }else if (all(gv %in% dfr[,"Name"])){
    colMeans(dfr[match(gv,dfr[,"Name"]),!np])
  }else{
    print("Something went wrong!")
  }
}


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

deidentify <- function(fdf){
  fdf[,!(colnames(fdf) %in% c("Name", "MRN"))]
}

`%!in%` <- Negate(`%in%`)

mesoID2num <- function(df){
  as.numeric(gsub("\\..*","",gsub("M|ML","",colnames(df))))
}


####CODES THAT DEPEND ON FUNCTIONS####

as.numeric(as.Date(mcd$LAST.RESPONS.DATE,'%m/%d/%Y')-as.Date(mcd$ENROLLED.DATE,'%m/%d/%Y'))
as.numeric(as.Date(mcd$Death_Date,'%m/%d/%Y')-as.Date(mcd$ENROLLED.DATE,'%m/%d/%Y'))
as.numeric(as.Date(mcd$OFF_STUDY.DATE,'%m/%d/%Y')-as.Date(mcd$ENROLLED.DATE,'%m/%d/%Y'))



nnc <- npc(rd) #make list of non-numeric columns in nanostring data
patientIDs <- gsub("M|ML","",colnames(rd[,!nnc]))
samplenumbers <- as.numeric(gsub("\\..*","",gsub("M|ML","",colnames(rd[,!nnc]))))

IFNGsig <- gene_sig(rd,IFNG.6.Gene[[2]],IFNG.6.Gene[[1]],nnc)
IFNGsn <- gsub("M|ML","",names(IFNGsig))
IFNGPID <- as.numeric(gsub("\\..*","",IFNGsn))
IFNGdf <- compiled[match(IFNGPID, compiled$Study),]
IFNGdf$Nano.String.Sample <- IFNGsn
IFNGdf$IFNG.Gene.Sig <- IFNGsig

IFNGdf$CD163 <- unlist(rd[match("CD163",rd[,"Name"]),!nnc])
IFNGdf$CD206 <- unlist(rd[match("MRC1",rd[,"Name"]),!nnc])
IFNGdf$CD274 <- unlist(rd[match("CD274",rd[,"Name"]),!nnc])
IFNGdf$CD8A <- unlist(rd[match("CD8A",rd[,"Name"]),!nnc])
IFNGdf$CSF2 <- unlist(rd[match("CSF2",rd[,"Name"]),!nnc])
IFNGdf$FOXP3 <- unlist(rd[match("FOXP3",rd[,"Name"]),!nnc])
IFNGdf$IL10 <- unlist(rd[match("IL10",rd[,"Name"]),!nnc])
IFNGdf$IL12A <- unlist(rd[match("IL12A",rd[,"Name"]),!nnc])

#Filtering down to necessary data
Plotthis <- IFNGdf %>%
  filter(Screen=="") %>%
  select(-Screen) %>%
  arrange(Nano.String.Sample) %>%
  filter(Study %!in% c(3,21,31,36,38)) %>%
  deidentify() #%>%


#Grouping
#If study time is greater than pfs split (PFS mean/median?) then add to high survival group
pfscutoff <- median(Plotthis$PFS[!is.na(Plotthis$PFS)])
PFS_Upper <- filter(Plotthis, Study.Time >= pfscutoff | PFS > pfscutoff)
PFS_Lower <- filter(Plotthis, PFS <= pfscutoff)

unique(PFS_Upper[,"Study"])

rdnum <- rd[,!nnc]
rownames(rdnum) <- rd[,"Name"]
#Filter out duplicate names before this point

rd_Upper <- rdnum[,match(unique(PFS_Upper[,"Study"]),mesoID2num(rdnum))]
rd_Lower <- rdnum[,match(unique(PFS_Lower[,"Study"]),mesoID2num(rdnum))]
rd_FULL <- cbind(rd_Upper,rd_Lower)
grp <- c(rep("HIGH",length(colnames(rd_Upper))),rep("LOW",length(colnames(rd_Lower))))

cds <-DGEList(rd_FULL,group=grp)
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
cds <- calcNormFactors(cds)
plotMDS.DGEList(cds,labels=colnames(cds$counts))
cds <- estimateCommonDisp(cds)
cds.com <- estimateCommonDisp(cds)
cds <- estimateTagwiseDisp( cds , prior.n = 10 )
cds.tgw <- estimateTagwiseDisp(cds)

meanVarPlot <- plotMeanVar( cds , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

de.cmn <- exactTest(cds.com, pair = c("HIGH","LOW"), dispersion = "auto")
de.tgw <- exactTest(cds.tgw, pair = c("HIGH","LOW"), dispersion = "auto")

topTags( de.cmn , n = 20 , sort.by = "p.value" )
rownames(topTags(de.cmn, n=20)$table)

resultsByFC <- topTags(de.cmn , n = nrow(de.cmn$table), sort.by = "logFC")$table
rownames(resultsByFC)[1:20]

qs <- qvalue(de.cmn$table$PValue, fdr.level = 0.05)


log2FC1 <- log2(rowMeans(rd_Upper))-log2(rowMeans(rd_Lower))
log2FC2 <- rowMeans(log2(rd_Upper))-rowMeans(log2(rd_Lower))
FC <- 2^log2FC1
log2FC1[order(abs(log2FC1), decreasing = TRUE)][1:20]

sort(exp(rowMeans(log(rd_Upper))) / exp(rowMeans(log(rd_Lower)))) #georowmean


design <- model.matrix(~factor(grp))
fit <- lmFit(log2(rd_FULL), design)
ebayesfit <- eBayes(fit)
tab <- topTable(ebayesfit, coef=2, adjust="fdr", n=150)
qsl <- qvalue(tab$adj.P.Val,fdr.level = 0.05)
tab$qvalues <- qsl$qvalues
hm <- heatmap(as.matrix(rd_FULL[rownames(tab),]), labCol = grp)
#heatmap(as.matrix(rd_FULL[rownames(tab),]), labCol = colnames(rd_FULL))
#unsupordertbl <- Plotthis[match(mesoID2num(rd_FULL)[tony$colInd],Plotthis$Study),]
#tab[order(tab$qvalues),][1:20,]

print(hm)


####PLOT CODES####
responsecolors <- c("blue","red","green","grey")
#Test Plot
#pfsdataonly <- Plotthis[is.na()]
#P <- ggplot(Plotthis,aes(x=IFNG.Gene.Sig,y=PFS)) 
#P + geom_point(aes(colour=factor(Best.Patient.Response)))
Plotthis$ResColors <- responsecolors[Plotthis$Best.Patient.Response]
#plot(Plotthis$CD274, Plotthis$CD8A)

for (x in 12:19){
  jpgfile <- paste0("C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\plots\\",colnames(Plotthis)[x],".jpeg")
  jpeg(file=jpgfile)
  plot(Plotthis[,colnames(Plotthis)[x]],Plotthis$PFS,
       col=Plotthis$ResColors,
       xlab=colnames(Plotthis)[x],
       ylab="PFS in months")
  
  
  dev.off()
}



####WORKING WITH CODES####
#duplicated(colnames(compiled))
#colnames(compiled[which(duplicated(colnames(compiled)))]) <- paste0(colnames(compiled[which(duplicated(colnames(compiled)))]),".dup")
#write.table(arrange(IFNGdf,Nano.String.Sample), "C:\\Users\\z33ky\\Documents\\UC stuff\\Meso proj\\Meso data compiled.txt", sep="\t", na = "",row.names = FALSE)
#avg_pfs <- mean(Plotthis$PFS[!is.na(Plotthis$PFS)])
#sum(!duplicated(PFS_lower$Study))
#setdiff(IFNGdf$Study,Plotthis$Study)
#PFS_Upper[which(!duplicated(PFS_Upper$Study)),"Study"]
#lotthis[is.na(Plotthis$PFS) & is.na(Plotthis$Study.Time),] #find which colums no PFS or study time
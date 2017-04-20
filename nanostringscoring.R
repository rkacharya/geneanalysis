#Scoring test script
#Written by Rajesh Acharya
#Uncomment below to install if necessary
##source('http://bioconductor.org/biocLite.R')
##biocLite('preprocessCore')
#load package
library(preprocessCore)
source("dataprocessing.R")
library(plotly)
library(dplyr)


####USER VARIABLES####

#==============================================================================================
#Genes of interest
geneset <- c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")#,"CCR5","PRF1","CXCL11","GZMA")
#File path. Use "file.choose(new=FALSE)" for file browser
##filename = file.choose(new=FALSE)
filename = "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\HNSCC\\CHI-02.HNSCC.noOutliers.NanostringData.QuantileNormalized.txt"
#File save path. Leave as "" for file browser or declare path.
filesave = ""
#Patient data file for response and PFS.
patientfile = "C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\Immunotherapy Patient List DeIdentified RA 20170127.txt"

#Parameter to be exported as tab delimited file
parameter = "" # subset(patientdata, select=c(1:14))

#OS days cutoff. Patients who are not tracked for a certain amount of time truncated
checkup_cutoff = 180 #days
#Turn Plots and file exports on/off
expo_data=FALSE #mean quantile normalized log10 data
expo_sigPFS_plot=FALSE #PFS/signature plot << Trigger to be implemented


#==============================================================================================



####DATA IMPORT AND CONDITIONING####
#input raw data count sheet
rawdata <- read.table(filename,
                      sep="\t",
                      header=TRUE
)
#exclude columns that have non sample data
nonpatientcol <- grepl('Code|Name|Accession', names(rawdata))
#create filtered data with genes of interest only
genesetdata <- rawdata[match(geneset, rawdata$Name), ]
#Transform nanostring data as needed. For Quantile Normalized use log
transformeddata <- genesetdata #log10datashift(genesetdata,nonpatientcol)
#Mean of signature scores
sigscore <- colMeans(transformeddata[,!nonpatientcol])

#Alternative sig using jounce numbers
jouncesigscores <- patientfileread("C:\\Users\\z33ky\\Documents\\UC stuff\\IFNg plot\\sigs_for_tanguy.txt")
rownames(jouncesigscores) <- as.character(jouncesigscores$Sample)
jouncesigscores <- jouncesigscores[,-1]
jouncesigscores <- t(jouncesigscores)
sigscores <- unlist(jouncesigscores["IFNG.6gene.signature",])
#sigscores <- sigscores/min(sigscores)
sigscore <- sigscores

#Import patient response data if needed "function:patientfileread"
patientimportdata <- patientfileread(patientfile)
#Match patientdata to names. Functions:substrRight
samplematch <- c()
for (x in 1:length(patientimportdata$Sample)){
  samplematch[x] <- match(as.numeric(patientimportdata$Sample[x]), 
                          as.numeric(substrRight(names(sigscore),2)))
  
}
#Append patient data with sigscores and IDs & subset
patientimportdata$Sig.Score <- sigscore[samplematch]
patientimportdata$Chi.ID <- names(sigscore[samplematch])

#Add Gene Values of Seq data to patient data
justnanostringdata <- genesetdata[,!nonpatientcol]
rownames(justnanostringdata) <- genesetdata$Name
for (x in rownames(justnanostringdata)){
  patientimportdata[x]<- unlist(justnanostringdata[x,])[samplematch]
}

#Truncate to only usable samples

#Add check to make sure each Patient identifier is unique after truncation.
patientdata <- patientimportdata %>%
  filter(Sample==Sample.to.use) %>% #Remove Samples we did not select to use
  filter(!is.na(Sample)) %>% # Remove samples without actual Samples
  filter(!is.na(Sig.Score)) # Remove Samples which did not have nanostring data

#Bind/filter survival days and follow-up category
patientdata <- patientdata %>%  
  cbind(as.data.frame(diff_in_days(patientdata$Start.Date,patientdata$Date.of.Progression,patientdata$Date.of.Death,patientdata$Date.of.Last.Follow.Up))) %>%
  filter(!((alv0dead1prog2 %in% c(0,1)) & (as.numeric(difftime(as.Date(Date.of.Last.Follow.Up,'%m/%d/%Y'),as.Date(Start.Date,'%m/%d/%Y'),units="days")) < checkup_cutoff)))
  #select(Patient.ID,Sample,Chi.ID,Sig.Score,Start.Date,Date.of.Progression,Date.of.Death,Date.of.Last.Follow.Up,`PFS..days.,surv`,alv0dead1prog2,Best.Overall.Response,HPV,Acquired.Resistance,)


# patientimportdata <- patientimportdata[as.numeric(patientimportdata$Sample)==as.numeric(patientimportdata$Sample.to.use),]
# patientimportdata <- patientimportdata[!is.na(patientimportdata$Sample)] #remove N/A columns

#Add flag for Using PFS or OS and put into set $Survival
#If OS too low and there is no PFS AND patient didn't die then remove point






#ADD LATER: list of non-matches from patient data side and nanostring data side
#patientdata <- subset(patientimportdata,!is.na(patientimportdata$Sig.Score),)
#Remove rows without PFS
#patientdata <- subset(patientdata, !is.na(patientdata$PFS..days.),)






####PLOT AND DATA EXPORT METHODS####
#Export Data
if (expo_data){
  if (filesave==""){
    print("Type a name for save file (*.txt) in dialogue box")
    savedir <- file.choose(new=TRUE)
  }
  else {
    savedir <- filesave
  }
  write.table(parameter, savedir,sep="\t")
}



####PLOTTING AND GRAPHICS####
#Plotly scatter of signature score vs PFS including patient response to treatment

#Turn into function.
#For loop, if new column contains characters for level then replace and relevel
Responselevels <- c("PD","SD","PR","CR","Not evaluable")
patientdata$Response.Flat <- patientdata$Best.Overall.Response
for (x in Responselevels){
  patientdata$Response.Flat[grep(x,patientdata$Response.Flat)] <- x
}
#as.factor(as.character(patientdata$Response.Flat))
plot_title <- "IFNG Signature and PFS" #expression(paste("IFN-",gamma, " Signature (GEP) and PFS"))
xlabel <- "Signature Score (Quantile Normalized and log10 Tranformed)"
ylabel <- "Progression Free Survival (days)"

patientdata$Response.Flat <- factor(patientdata$Response.Flat, levels=Responselevels)

shapes <- c("triangle-up","triangle-up","square","square","circle")
shapes <- setNames(shapes, levels(patientdata$Response.Flat))
pal <- c("#8067fc","#4a26ff","#ff2e26","#ff6851","#00ffd4")#Might have to <-setNames(pal,Responselevels)
pal <- setNames(pal, levels(patientdata$Response.Flat))
sizing <- (!is.na(patientdata$Acquired.Resistance))*5+10
sizing <- setNames(sizing, names(patientdata$Response.Flat))
patientdata$shapes <- shapes[patientdata$Response.Flat]
patientdata$pal <- pal[patientdata$Response.Flat]
patientdata$sizing <- sizing
patientdata$AcqRes[1:length(patientdata$Acquired.Resistance)] <- ""
patientdata$AcqRes[patientdata$Acquired.Resistance==TRUE] <- "</br>Acquired Resistance"
conditions <- c("</br>Followed: ", "</br>OS: ","</br>PFS: ")
conditions_ind <- c(0,1,2)
patientdata$survtxt <- conditions[match(patientdata$alv0dead1prog2,conditions_ind)]


hovertext= paste('Sample: ',patientdata$Chi.ID,
                 '</br>Patient ID: ', patientdata$Patient.ID,
                  '</br>Score: ',round(patientdata$Sig.Score,2),
                  patientdata$survtxt ,patientdata$surv, ' days',
                 '</br>Response: ',patientdata$Best.Overall.Response
                 ,
                  patientdata$AcqRes)

#Graphic Color Pallet
color1 = "#1d2120" #usually plot BG
color2 = "#5a5c51" #usually color surrounding plot
color3 = "#ba9077" #Text color
color4 = "#bcd5d1" #grid colors
color5 = "" #Rectangle 1
color6 = "" #Rectangle 2
color7 = "" #Rectangle 3
color8 = "" #Do androids dream of electric sheep?
color9 = "" #Danger Zone

sigPFSplot <- plot_ly(patientdata,x = ~Sig.Score, y = ~surv,type='scatter',mode='markers',
                      hoverinfo='text',
                      text=hovertext,
                      symbol=~Response.Flat, symbols=shapes,
                      color=~Response.Flat, colors=pal,
                      #size=~sizing,
                      marker=list(size=12),
                      showlegend=T
                      ) %>%
  # add_trace(patientdata, x=~Sig.Score, y=~PFS..days., type='scatter',mode='markers',
  #           visible="legendonly", showlegend=T, legendgroup="Response"
  # 
  #           ) %>%

  layout(title=plot_title,
         titlefont=list(
           color=color4, 
           size="20", 
           family="Arial"),
         paper_bgcolor=color2,
         plot_bgcolor=color1,
         
         margin=list(
           l=80, r= 120, b=80, t=100, pad=1),
         
         xaxis = list(
           gridcolor = color3,
           #linecolor = "color" #sets axis line color
           tickfont = list(
             family="Arial",
             size="14",
             color=color4),
           #nticks = integer #max number ticks
           #linewidth = #width of axis line
           #autorange = TRUE #TRUE|FALSE|"reversed"
           title=xlabel,
           titlefont=list(
             family="Arial",
             size="14",
             color=color4)),
         
         yaxis = list(
           gridcolor = color3,
           tickfont = list(
             family="Arial",
             size = 8,
             color=color4),
           title=ylabel,
           titlefont=list(
             family="Arial",
             size="14",
             color=color4)),
         
         # shapes = list(
         #   list(type="rect",
         #        x0="", x1="", xref="x")
         #   list(type="rect",
         #        x0="", x1="," xref="x")),
         
         legend = list(
           font = list(
             family = "Arial",
             size = 12,
             color = color4),
           bgcolor = color1,
           bordercolor = "#000000",
           borderwidth = 1)
             
           )

print(sigPFSplot)           
print("Script Completed")



#workin on codes
##fullfilepath <- unlist(strsplit(filename,'\\\\'))
#match(unlist(strsplit(as.character(patientimportdata$Sample[3]), ",")), substr(colnames(rawdata),nchar(colnames(rawdata))-2+1,nchar(colnames(rawdata))))

#finds elements of first vector and see if they exist in second
#setdiff(as.numeric(substrRight(names(sigscore),2)), patientimportdata$Sample)

#as.numeric(difftime(as.Date(patientdata$Start.Date[1],'%m/%d/%Y'),as.Date(patientdata$Date.of.Progression[1],'%m/%d/%Y'), units="days"))

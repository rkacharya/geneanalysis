#Convert G group list to readable list format
#Created by Rajesh Acharya
#Made: 2017,Mar,17. Update: 



####READ FILES####
Gdf <- read.table("C:\\Users\\z33ky\\Documents\\UC stuff\\HLA Stuff\\HLA_nom_g_grps.txt",
                  header = FALSE,
                  sep = "")

GL <- strsplit(as.character(Gdf[,1]), ";")

idf <- read.table("C:\\Users\\z33ky\\Documents\\UC stuff\\HLA Stuff\\HLAVBseq_out_normal.txt",
                  header = FALSE,
                  sep = "\t")

####FUNCTIONS####

matchGgrp <- function(gene,fst,snd){
  matchstr <- paste(gene, fst, ":", snd, sep = "")
  for (ii in 1:length(GL)){
    tmp1 <- strsplit(strsplit(GL[[ii]][2],"/")[[1]],":")
    tmp2 <- vector(mode = "character", length = length(tmp1))
    for (jj in 1:length(tmp1)){
      tmp2[jj] <- paste(GL[[ii]][1],tmp1[[jj]][1],":",tmp1[[jj]][2], sep="")
    }
    didstrmtch <- grepl(matchstr, tmp2)
    if (any(didstrmtch)){
      outpt <- paste(gene, GL[[ii]][3], sep="")#Output matching G name
    }else{
      outpt <- paste(gene, fst,":", snd, sep="")
    }
  }
  return(outpt)
}

#convert from HLA nomelcature 3 back to 2. Input as.character
conv_nom_3to2 <- function(gene, fst, snd, G){
  
  if (missing(G)){
    newset <- c(gene, fst, snd)
  }else if (G %in% c("G","g")){
    newset <- c(gene, fst, snd, "g")
  }else{
    print("Error in G group format")
    newset <- c(gene, fst, snd, "ERROR")
  }
  
  if ((gene=="A*") && (fst=="02")){
    newfst <- "92"
    newsnd <- as.character(as.numeric(snd)-100)
  }else if((gene=="B*") && (fst=="15")){
    newfst <- "95"
    newsnd <- as.character(as.numeric(snd)-100)
  }else{
    newfst <- fst
    newsnd <- snd
  }
  newset[2] <- newfst
  newset[3] <- newsnd
  return(newset)
}




####MAIN CODES####

#Re-write table with G groups
wdf <- idf
for (ii in 3:8){
  for(jj in length(idf[1])){
    tmp1 <- strsplit(wdf[ii,jj])
    wdf[jj,ii]
  }
}

#Convert HLA_nom 3 to HLA_nom 2



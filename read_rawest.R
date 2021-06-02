library(stringr)

############# read cond table
cond1 <- read.csv("data/cond1.csv", stringsAsFactors = FALSE)
cond2 <- read.csv("data/cond2.csv", stringsAsFactors = FALSE)
cond3 <- read.csv("data/cond3_qpcr.csv", stringsAsFactors = FALSE)
cond4 <- read.csv("data/cond4.csv", stringsAsFactors = FALSE)
ob_cond <- merge(cond1,cond2)
ob_cond <- merge(ob_cond, cond3, all.x=TRUE)
ob_cond <- merge(ob_cond, cond4)

ob_cond <- ob_cond[order(ob_cond$libname),]
rownames(ob_cond) <- ob_cond$libname

## remove 2 failed libs
keep <- !(rownames(ob_cond) %in% c("N704_S506", "N704_S507"))
ob_cond <- ob_cond[keep,]

ob_cond$infected <- factor(as.character(ob_cond$infected))
ob_cond$donorid <- factor(ob_cond$donorid)
ob_cond$Treatment <- factor(ob_cond$Treatment)
ob_cond$covidqpcr <- log10(1+ob_cond$covidqpcr)
ob_cond$sex <- factor(ob_cond$sex)
ob_cond$ismock <- ob_cond$Treatment=="none" & ob_cond$infected=="FALSE"

ob_cond$inflevel2 <- factor(ob_cond$inflevel2)

ob_cond <- ob_cond[mixedorder(ob_cond$Umuwell),]


############ read counts; sum up transcripts to genes

read_featurecount <- function(fname){
  dat <- read.table(fname, stringsAsFactors = FALSE,sep="\t", header = TRUE, row.names = 1)
  dat <- dat[,-(1:5)]
  dat
}


if(FALSE){
  ################## kallisto style
  ob_counts <- read.csv("data/counts.kallisto.csv", stringsAsFactors = FALSE, row.names = 1)
  sp <- str_split_fixed(colnames(ob_counts), "_", n = 3)
  colnames(ob_counts) <- sprintf("%s_%s", sp[,1], sp[,2])
  
  #Make cond and counts to match up
  ob_cond <- ob_cond[rownames(ob_cond) %in% colnames(ob_counts),]
  ob_counts <- ob_counts[,rownames(ob_cond)]
  
  #sort out kallisto names
  rownames(ob_counts) <- str_split_fixed(rownames(ob_counts),"\\.",2)[,1]
  
  #merge kallisto transcripts to genes
  ob_counts$transcript <- rownames(ob_counts)
  dd <- melt(ob_counts, id.vars = "transcript")
  dd <- merge(dd,map_ens_trans)
  dd <- sqldf("select variable, sum(value) as cnt, geneid from dd group by variable, geneid")
  
  dd <- dcast(dd, variable~geneid, value.var = "cnt")
  rownames(dd) <- dd$variable
  dd <- dd[,-1]
  ob_counts <- round(t(dd))
  remove(dd)
  
} else {
  ##################### STAR style
  
  ob_counts <- read_featurecount("data/counts.star.csv")
  
  #sum up for all 4 lanes in each sample
  summat <- matrix(nrow=nrow(ob_counts), ncol=ncol(ob_counts)/4)
  for(i in (1:(ncol(ob_counts)/4))-1){
    istart <- i*4+1
    iend <- i*4+4
    summat[,i+1] <- rowSums(ob_counts[istart:iend])
  }
  colnames(summat) <- colnames(ob_counts)[((1:ncol(summat))-1)*4+1]  
  rownames(summat) <- rownames(ob_counts)
  
  #Name properly
  sp <- str_split_fixed(colnames(summat), "_", n = 3)
  colnames(summat) <- sprintf("%s_%s", sp[,1], sp[,2])
  ob_counts <- summat
  
  #Make cond and counts to match up
  ob_cond <- ob_cond[rownames(ob_cond) %in% colnames(ob_counts),]
  ob_counts <- ob_counts[,rownames(ob_cond)]
}




write.csv(ob_cond, "data/final_cond.csv")
write.csv(ob_counts, "data/final_counts.csv")


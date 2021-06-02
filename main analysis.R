library(biomaRt)
library(GO.db)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(reshape2)
library(sqldf)
library(Rtsne)
library(viridis)
library(gtools)
library(DESeq2)

###########################################################################
##################### Helper functions ####################################
###########################################################################


unrep <- function(r){
  r[!duplicated(r$symbol),]
}




###########################################################################
##################### Load the data #######################################
###########################################################################


list_ylva <- read.csv("ylva.csv",stringsAsFactors = FALSE)
colnames(list_ylva) <- c("ylva_uniprot","symbol")

list_hya <- read.csv("list_hya.csv",stringsAsFactors = FALSE)

list_markergenes <- read.csv("list_markergenes.csv",stringsAsFactors = FALSE)


map_ens_trans <- read.csv("data/ensid_transcript.csv", stringsAsFactors = FALSE, sep="\t")
colnames(map_ens_trans) <- c("geneid","transcript")

ensid <- read.csv("genesym_human.csv", stringsAsFactors = FALSE,sep="\t")
colnames(ensid) <- c("geneid","symbol")

ensid_covid <- read.csv("genesym_covid.csv", stringsAsFactors = FALSE,sep="\t")
colnames(ensid_covid) <- c("geneid","symbol")

ensid <- rbind(ensid,ensid_covid)


genelocation <- read.csv("data/gene_location.csv", stringsAsFactors = FALSE,sep="\t")
colnames(genelocation) <- c("geneid","chrom","transid")

sex_geneid <- unique(genelocation$geneid[genelocation$chrom %in% c("X","Y")])


## Read the main data
ob_cond   <- read.csv("data/final_cond.csv")
rownames(ob_cond) <- ob_cond[,1]
ob_cond <- ob_cond[,-1]
ob_counts <- read.csv("data/final_counts.csv")
rownames(ob_counts) <- ob_counts[,1]
ob_counts <- ob_counts[,-1]

## keep total counts
ob_cond$totc <- colSums(ob_counts)


## reduced ensid list for our count table
ensid_red <- ensid[ensid$geneid %in% rownames(ob_counts),]
ensid_red <- ensid_red[!duplicated(ensid_red$geneid),]
rownames(ensid_red) <- ensid_red$geneid


###########################################################################
##################### Deal with virus qPCR ################################
###########################################################################

### the qPCR data makes most sense on log-scale
#hist(ob_cond$covidqpcr[!is.na(ob_cond$covidqpcr)])
#hist(log10(1+ob_cond$covidqpcr[!is.na(ob_cond$covidqpcr)]))

######################### Quick check of levels among donors
#x <- ob_cond[ob_cond$Treatment=="none" & ob_cond$donorid!="o" & ob_cond$infected=="TRUE",]
x <- ob_cond[ob_cond$donorid!="o" & ob_cond$infected=="TRUE",]

x$donor_inf <- sprintf("%s_%s",x$donorid, x$Treatment)
p<-ggplot(data=x, aes(x=donor_inf, y=covidqpcr, fill=sex)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 45))
p
ggsave("out/out.covidlevel/covidlev_qpcr.pdf",p)


#### Measure of covid infection
dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
ob_cond$covidrnaseq <- log10(1+colSums(dat[ensid_covid$geneid,]))
ob_cond$covidrnaseq[ob_cond$infected=="FALSE"] <- 0
ob_cond$log_newpcr <- log10(1+ob_cond$newpcr)

pdf("out/out.covidlevel/qpcr vs rnaseq.pdf")
plot(
  ob_cond$covidrnaseq, log10(1+ob_cond$newpcr),
  xlab="log10(1+sum normalized RNAseq covid genes)",ylab="log10(1+rtQPCR readout)")
dev.off()


#### Define a severity score for each donor
x <- ob_cond[ob_cond$donorid!="o" & ob_cond$infected=="TRUE" & ob_cond$Treatment=="none",]
#ob_cond <- merge(ob_cond,sqldf("select donorid, avg(covidqpcr) as severity from x group by donorid"), all.x=TRUE)
ob_cond <- merge(ob_cond,sqldf("select donorid, avg(newpcr) as severity_pcr, avg(covidrnaseq) as severity_rnaseq from x group by donorid"), all.x=TRUE)
rownames(ob_cond) <- ob_cond$libname
ob_cond <- ob_cond[colnames(ob_counts),]



###########################################################################
##################### Write data for other users ##########################
###########################################################################


dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
rownames(ensid_red) <- ensid_red$geneid
rownames(dat) <- ensid_red[rownames(dat),]$symbol
colnames(dat) <- ob_cond$Umuwell

write.csv(dat, file = "out/out.counts/all_counts_sf.csv")
write.csv(ob_cond, file = "out/out.counts/all_cond.csv")
write.csv(counts(dds, normalized=FALSE), file = "out/out.counts/all_counts_raw.csv")



###########################################################################
################## GO stuff ###############################################
###########################################################################

list_genes_ecm <- unique(c(
  read.csv("QuickGO-0005615.tsv",sep="\t",stringsAsFactors = FALSE)$SYMBOL,
  read.csv("QuickGO-0031012.tsv",sep="\t",stringsAsFactors = FALSE)$SYMBOL))

list_genes_hyaluronan_metabolic_process <- read.csv("GO-0030212.tsv",sep="\t",stringsAsFactors = FALSE)$SYMBOL








###########################################################################
##################### Individual gene exp plots ###########################
###########################################################################


##############################################################
## Function to produce a heatmap
plot_heatmap_forgenes <- function(dat, glist, keep){
  
  ## Name the samples
  colnames(dat) <- sprintf("%s_%s_%s_%s_%s",ob_cond$Umuwell, ob_cond$infected,ob_cond$Treatment,ob_cond$donorid, ob_cond$sex)
  
  ## Subset the samples
  subgenes <- ensid_red[ensid_red$symbol %in% glist,]
  subdat <- dat[subgenes$geneid,keep]
  
  ## Name by gene symbol  
  mm <- as.matrix(subdat)
  rownames(mm) <- subgenes$symbol
  mm <- mm[,ncol(mm):1]
  
  ## Rescale
  for(i in 1:nrow(mm)){
    # mm[i,] <- mm[i,] - min(mm[i,])
    # mm[i,] <- mm[i,]/max(mm[i,])
    mm[i,] <- mm[i,]/max(abs(mm[i,]))
  }
  mm[is.na(mm)]<-0
  
  mm <- mm[,mixedorder(colnames(mm),decreasing = TRUE)]
  
  mmm <- melt(mm)
  colnames(mmm) <- c("X1","X2","value")
  mmm$X1 <- factor(as.character(mmm$X1), levels = glist)
  mmm$X2 <- factor(as.character(mmm$X2), levels = colnames(mm))
  
  p <- ggplot(mmm, aes(X1, X2, fill= value)) + 
    geom_tile() + theme(axis.text.x = element_text(angle=90)) + scale_fill_gradient2()
  # if(!is.null(fname)){
  #   ggsave(fname, p)
  # } 
  
  p
  # heatmap(
  #   Colv = NA, Rowv = NA, 
  #   t(mm),
  #   #scale = "col",
  #   col= colorRampPalette(brewer.pal(8, "Blues"))(25))
}



## Lists of marker genes
unique(list_markergenes$genelist)

## SF-normalize
dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)
dat_plain <- log10(1+counts(dds, normalized=TRUE))



######################################################################
## Heatmaps for all marker genes
for(glistname in unique(list_markergenes$genelist)){
  #pdf(sprintf("out/out.heatmap.plain/heatmap_tubes_%s.pdf",glistname), width = 15)
  plot_heatmap_forgenes(
    dat_plain,
    list_markergenes$symbol[list_markergenes$genelist==glistname], 
    ob_cond$donorid=="o")
  #dev.off()
  ggsave(sprintf("out/out.heatmap.plain/heatmap_tubes_%s.pdf",glistname))
  
  #pdf(sprintf("out/out.heatmap.plain/heatmap_planes_%s.pdf",glistname), width = 15)
  plot_heatmap_forgenes(
    dat_plain,
    list_markergenes$symbol[list_markergenes$genelist==glistname], 
    ob_cond$donorid!="o" & ob_cond$Treatment=="none")
  ggsave(sprintf("out/out.heatmap.plain/heatmap_planes_%s.pdf",glistname))
  #dev.off()
}



###########################################################################
######### Individual gene exp plots -- with mock subtraction ##############
###########################################################################


## SF-normalize
dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)

## Remove mock
dat <- log10(1+counts(dds, normalized=TRUE))
dat_nomock <- log10(1+counts(dds, normalized=TRUE))
for(i in 1:ncol(dat_nomock)){
  mock_i <- ob_cond$donorid==ob_cond$donorid[i] & ob_cond$Treatment=="none" & ob_cond$infected=="FALSE"
  dat_nomock[,i] <- dat_nomock[,i] - rowMeans(dat[,mock_i,drop=FALSE])
}
## Only keep donors for which there is a mock left... turns out to be all of them
#donor_with_mock <- unique(ob_cond$donorid[ob_cond$Treatment=="none" & ob_cond$infected=="FALSE"])


#################################################
## Heatmaps for all marker genes, mock subtracted
for(glistname in unique(list_markergenes$genelist)){
  
  #pdf(sprintf("out/out.heatmap.mocksubtracted/heatmap_tubes_%s.pdf",glistname), width = 15)
  plot_heatmap_forgenes(
    dat_nomock,
    list_markergenes$symbol[list_markergenes$genelist==glistname], 
    ob_cond$donorid=="o" & !ob_cond$ismock)
  #dev.off()
  ggsave(sprintf("out/out.heatmap.mocksubtracted/heatmap_tubes_%s.pdf",glistname))
  
  
  #pdf(sprintf("out/out.heatmap.mocksubtracted/heatmap_planes_%s.pdf",glistname), width = 15)
  plot_heatmap_forgenes(
    dat_nomock,
    list_markergenes$symbol[list_markergenes$genelist==glistname], 
    ob_cond$donorid!="o" & !ob_cond$ismock & ob_cond$Treatment=="none")
  #dev.off()
  ggsave(sprintf("out/out.heatmap.mocksubtracted/heatmap_planes_%s.pdf",glistname))
  
}



###########################################################################
##################### Covid levels for each donor #########################
###########################################################################

dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)
dat <- log10(1+counts(dds, normalized=TRUE))

##Only keep covid genes
dat <- dat[ensid_covid$geneid,]

for(donorid in unique(ob_cond$donorid)){
  
  #keep <- ob_cond$donorid==donorid & ob_cond$Treatment %in% c("none","beta03","beta3")
  keep <- ob_cond$donorid==donorid & ob_cond$Treatment=="none"
  #keep <- ob_cond$donorid==donorid & ob_cond$infected=="TRUE" & ob_cond$Treatment=="none"
  
  dd <- melt(dat[,keep])#, id.vars = "transcript")
  colnames(dd) <- c("geneid","libname","transcript")
  dd <- merge(dd, ensid)
  dd <- merge(dd, ob_cond)
  
  p<-ggplot(data=dd, aes(x=symbol, y=transcript, fill=Umuwell)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 45))+
    ylim(0,7)
  p
  ggsave(sprintf("out/out.covidlevel/untreated_%s.pdf",donorid),plot = p)
}


for(donorid in c("o")){
  
  keep <- ob_cond$donorid==donorid & ob_cond$Treatment %in% c("beta03","beta3")
  #keep <- ob_cond$donorid==donorid & ob_cond$Treatment=="none"
  #keep <- ob_cond$donorid==donorid & ob_cond$infected=="TRUE" & ob_cond$Treatment=="none"
  
  dd <- melt(dat[,keep])#, id.vars = "transcript")
  colnames(dd) <- c("geneid","libname","transcript")
  dd <- merge(dd, ensid)
  dd <- merge(dd,ob_cond)
  
  p<-ggplot(data=dd, aes(x=symbol, y=transcript, fill=Umuwell)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 45))+
    ylim(0,7)
  p
  ggsave(sprintf("out/out.covidlevel/treated_%s.pdf",donorid),plot = p)
}








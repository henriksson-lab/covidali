# Note - this code is not used for the analysis presented in the paper


library(NMF)
library(fastICA)
#BiocManager::install("fgsea")
#http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
library(fgsea)
library(org.Hs.eg.db)


###############################################################################
########## get normalized counts and dispersion ###############################
###############################################################################

####
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,ob_cond$donorid!="o"], colData = ob_cond[ob_cond$donorid!="o",], design = ~1)
dds <- DESeq(dds)
dat <- getVarianceStabilizedData(dds)
colnames(dat) <- ob_cond[ob_cond$donorid!="o",]$Umuwell
rownames(ensid_red) <- ensid_red$geneid
rownames(dat) <- ensid_red[rownames(dat),]$symbol


############## get normalized counts, rescaled
nc <- counts(dds, normalized=TRUE)
nc <- log10(1+nc)
mean(rowMeans(nc)>0.01)
nc <- nc[rowMeans(nc)>0.01,]
nc <- t(scale(t(nc),center=FALSE))

############## Find highly dispersed genes


#plotDispEsts(dds)

dispGene <- mcols(dds)$dispGeneEst
dispFit <- mcols(dds)$dispFit

plot(log10(1+rowMeans(counts(dds,normalized=TRUE))),log10(dispGene)-log10(dispFit))
plot(log10(1+rowMeans(counts(dds,normalized=TRUE))),log10(dispGene))
plot(log10(1+rowMeans(counts(dds,normalized=TRUE))),log10(dispFit))

dispDiff <- log10(dispGene)-log10(dispFit)
dispDiff[is.na(dispDiff)]<-0
plot(log10(1+rowMeans(counts(dds,normalized=TRUE))),dispDiff)

#rownames(dat)[dispDiff> 1.5]
#rownames(dat)[order(dispDiff, decreasing = TRUE)[1:100]]

##############

#Pull out genes with particularly high dispersion?
#subdat <- dat[rowMeans(nc)>0.02,]


###############################################################################
############################ nmf #############################################
###############################################################################

subdat <- dat[order(dispDiff, decreasing = TRUE)[1:2000],]

#b[order(b[,5],decreasing = TRUE)[1:100],1]
#rownames(b)[order(b[,1],decreasing = TRUE)[1:100]]

nmfTopGenesMatrix <- function(res, topn=200){
  b <- basis(res)
  mo <- matrix(nrow=topn, ncol=ncol(b))
  for(i in 1:ncol(b)){
    mo[,i] <- rownames(b)[order(b[,i],decreasing = TRUE)[1:topn]]
  }
  colnames(mo) <- sprintf("base_%s",1:ncol(b))
  mo
}

res5 <- nmf(subdat, 5)
write.csv(nmfTopGenesMatrix(res5),"our_bulk/out.nmf/base5.csv")
write.csv(nmfTopGenesMatrix(nmf(subdat, 10)),"our_bulk/out.nmf/base10.csv")


rowSums(coefficients(res5))


#Can pull out the top 100 genes in each list; only do a heatmap from these for all these
#basismap(res, subsetRow=TRUE)


#res.multi.method <- nmf(subdat, 3, list('brunet', 'lee', 'ns'), seed=123456, .options='t')
#compare(res.multi.method)



# estimate quality measures from the shuffled data (use default NMF algorithm)
#V.random <- randomize(subdat)
#estim.r.random <- nmf(V.random, 2:6, nrun=10, seed=123456)
# plot measures on same graph
#plot(estim.r, estim.r.random)



###############################################################################
############################ ica #############################################
###############################################################################


icaTopGenesMatrix <- function(res, topn=nrow(res$S)){
  b <- abs(res$S)
  mo <- matrix(nrow=topn, ncol=ncol(b))
  for(i in 1:ncol(b)){
    mo[,i] <- rownames(b)[order(b[,i],decreasing = TRUE)[1:topn]]
  }
  colnames(mo) <- sprintf("base_%s",1:ncol(b))
  mo
}


doGseaOnRanks <- function(exampleRanks){
  hs <- org.Hs.eg.db
  ov <- select(hs, 
               keys = names(exampleRanks),
               columns = c("ENTREZID", "SYMBOL"),
               keytype = "SYMBOL")
  ov <- merge(ov,data.frame(SYMBOL=names(exampleRanks),value=as.double(exampleRanks)))
  ov <- ov[!is.na(ov$ENTREZID),]
  ovm <- ov$value
  names(ovm) <- ov$ENTREZID
  
  pathways <- reactomePathways(names(ovm))
  fgseaRes <- fgsea(pathways, ovm, maxSize=500, nperm=1000)
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
  # topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  # topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ovm, fgseaRes, gseaParam=0.5)
  fgseaRes
}


subdat.ica <- dat[order(dispDiff, decreasing = TRUE)[1:10000],]
#subdat.ica <- dat[rowMeans(nc)>sort(rowMeans(nc),decreasing = TRUE)[10000],]
dim(subdat.ica)

for(i in 1:6){
  print(i)
  ti <- fastICA(subdat.ica, i, alg.typ = "parallel", fun = "logcosh", alpha = 1,
                method = "R", row.norm = FALSE, maxit = 200,
                tol = 0.0001, verbose = TRUE)
  write.csv(icaTopGenesMatrix(ti),sprintf("our_bulk/out.ica/ica_%s.csv",i))
  
  for(j in 1:i){
    pdf(sprintf("our_bulk/out.ica/ica_%s_reactome_%s.pdf",i,j), width = 12)
    doGseaOnRanks(ti$S[,j])  #hya comes out if i=5
    dev.off()
  }
  
}



i<-5
ti <- fastICA(subdat.ica, i, alg.typ = "parallel", fun = "logcosh", alpha = 1,
              method = "R", row.norm = FALSE, maxit = 200,
              tol = 0.0001, verbose = TRUE)
tg <- icaTopGenesMatrix(ti)
View(tg)
dev.off()

dev.off()
doGseaOnRanks(ti$S[,2])

dev.off()
doGseaOnRanks(ti$S[,3])

dev.off()
doGseaOnRanks(ti$S[,4])

dev.off()
doGseaOnRanks(ti$S[,5])






###############################################################################
############################ mara #############################################
###############################################################################

dat_sitecount <- read.csv("mara/mara_sitecount.csv",sep="\t")
rownames(ensid_red) <- ensid_red$geneid
dat_sitecount <- dat_sitecount[rownames(dat_sitecount) %in% ensid_red$geneid,]
rownames(dat_sitecount) <- ensid_red[rownames(dat_sitecount),]$symbol

plot(sort(sqrt(rowVars(dat))))

subdat.mara <- log10(1+counts(dds, normalized=TRUE))
subdat.mara <- subdat.mara[order(dispDiff,decreasing = TRUE),][1:3000,]
#  dat#[order(dispDiff, decreasing = TRUE)[1:10000],]
#subdat.mara <- subdat.mara[sqrt(rowVars(subdat.mara))>0.5,]


rownames(ensid_red) <- ensid_red$geneid
rownames(subdat.mara) <- ensid_red[rownames(subdat.mara),]$symbol
subdat.mara <- subdat.mara[rownames(subdat.mara) %in% rownames(dat_sitecount),]


dat_sitecount.mara <- dat_sitecount[rownames(subdat.mara),]

write.table(dat_sitecount.mara,"mara/send/mara_sitecount.csv",sep="\t",quote = FALSE)
write.table(subdat.mara,"mara/send/mara_exp.csv",sep="\t",quote = FALSE)

dim(dat_sitecount.mara)
dim(subdat.mara)


###############################################################################
############################ mara analysis #############################################
###############################################################################

dat_act <- read.csv("mara/ismara_report/activity_table",sep="\t")
dat_act

dat_mat <- read.csv("mara/ismara_report/active_matrices",sep="\t")
dat_mat

dat_delta <- read.csv("mara/ismara_report/delta_table",sep="\t")
dat_delta

delta_cond <- ob_cond[ob_cond$Umuwell %in% rownames(dat_delta),]
rownames(delta_cond) <- delta_cond$Umuwell
rownames(delta_cond) == rownames(dat_delta)


rownames(delta_cond)


dispDiff


dd <- colMeans(dat_delta[delta_cond$infected=="TRUE" & delta_cond$Treatment=="none",]) - 
  colMeans(dat_delta[delta_cond$infected=="FALSE" & delta_cond$Treatment=="none",])
dd <- data.frame(tf=names(dd), diff=dd, stringsAsFactors = FALSE)
dd <- dd[order(dd$diff),]
head(dd$tf,n=20)
tail(dd$tf,n=20)

### Attempt at some basic network analysis. Not used for the paper


###########################################################################
##################### Clustering of samples ###############################
###########################################################################


dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)
dat <- log10(1+counts(dds, normalized=TRUE))
colnames(dat) <- sprintf("%s_%s_%s_%s",ob_cond$infected,ob_cond$Treatment,ob_cond$donorid, ob_cond$sex)
hist(as.double(rowMeans(dat)),breaks = 100)
dat <- dat[rowMeans(dat)>0.2,]


### ordinary hclust
pdf("out/out.clust/samples_hclust.pdf", width = 15, height = 15)
plot(hclust(dist(t(dat))))
dev.off()



### tSNE
set.seed(1)
tsne <- Rtsne(t(dat), dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)

pdf("out/out.clust/samples_tsne.pdf")
#use_for_col <- ob_cond$donorid
use_for_col <- ob_cond$infected
#use_for_col <- ob_cond$Treatment
colors = rainbow(length(unique(use_for_col)))
plot(tsne$Y, t='n', main="tSNE", xlab="tSNE 1", ylab="tSNE 2", "cex.main"=2, "cex.lab"=1.5)
text(tsne$Y, labels=colnames(dat),cex=0.5, col=colors[use_for_col])
dev.off()




###########################################################################
##################### Clustering of genes #################################
###########################################################################



dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~1)
dds <- estimateSizeFactors(dds)
dat <- log10(1+counts(dds, normalized=TRUE))
#hist(as.double(rowMeans(dat)),breaks = 100)
dat <- dat[rowMeans(dat)>0.2,]
rownames(ensid_red) <- ensid_red$geneid
rownames(dat) <- ensid_red[rownames(dat),]$symbol


dat <- as.data.frame(dat) %>% distinct()  ### not sure which genes are kicked out, but some annoying ones?  300 genes
dat <- t(scale(t(dat)))

set.seed(1)
tsne <- Rtsne(dat, dims = 2, perplexity=50, verbose=TRUE, max_iter = 500)

pdf("out/out.clust/tsne gene network.pdf", width = 20, height = 20)
thecol <- rep("black", nrow(dat))
#thecol[rownames(dat) %in% ylva$symbol] <- "blue"
thecol[rownames(dat) %in% list_hya$symbol] <- "blue"
thecol[rownames(dat) %in% c("ISG15","MX1","MX2","IFITM1","IFITM3","RARRES1","PLAU","TNFAIP2","STAT1","SAA1","SAA2","CXCL2","CCL20","NFKBIZ","CXCL1")] <- "red"
plot(tsne$Y, pch=19, cex=0, col=thecol)
#points(tsne$Y[thecol!="black",], pch=19, cex=0.7, col=thecol[thecol!="black"])
#text(tsne$Y, labels=rownames(dat),cex=0.8)
text(tsne$Y[thecol!="black",], labels=rownames(dat)[thecol!="black"],cex=0.8, col=thecol[thecol!="black"])
dev.off()






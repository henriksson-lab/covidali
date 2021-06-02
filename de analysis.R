#### this code assumes that "main analysis.R" has been run

###################################################################################################################
##################### DE testing ##################################################################################
###################################################################################################################


plotVolcano <- function(r, title="", removesex=TRUE, minexp=50,removecovid=TRUE){

  listcolor <- data.frame(
    symbol=list_ylva$symbol,
    stringsAsFactors = FALSE
  )  
  listcolor <- data.frame(
    symbol=list_markergenes$symbol[list_markergenes$genelist=="hya_go"],
    stringsAsFactors = FALSE
  )  
  listcolor <- data.frame(
    symbol=list_genes_ecm,
    stringsAsFactors = FALSE
  )  
  listcolor$col <- "red"
  
  
  r <- r[!is.na(r$padj),]
  if(removesex){
    r <- r[!(r$geneid %in% sex_geneid),] 
  }
  if(removecovid){
    r <- r[!(r$geneid %in% ensid_covid$geneid),] 
  }
  r <- r[r$baseMean>=minexp,]
  r <- unrep(r)  
  
  r<-merge(r,listcolor, all.x=TRUE)
  r$col[is.na(r$col)]<-"black"
  r$col[r$symbol %in% list_hya$symbol] <- "green"
  
  sp <- r$col!="black"
  plot(r$log2FoldChange, -log10(r$pvalue),cex=0, main=title)
  text(r$log2FoldChange[!sp], -log10(r$pvalue)[!sp], labels = r$symbol[!sp], cex=0.8, col=r$col[!sp])
  text(r$log2FoldChange[sp], -log10(r$pvalue)[sp], labels = r$symbol[sp], cex=0.8, col=r$col[sp])
}


plotVolcano.save <- function(r, title, removesex=FALSE, minexp=50,removecovid=FALSE){
  plotVolcano(r, title, removesex, minexp, removecovid)
  pdf(sprintf("out/out.de/all_%s.pdf",title), width = 15, height = 15)
  plotVolcano(r, title, removesex)
  dev.off()
  write.csv(
    unrep(r)[,c("geneid","baseMean","log2FoldChange","pvalue","padj","symbol")],
    sprintf("out/out.de/all_%s.csv",title),row.names = FALSE)
  
  ### Secretome
  r_sub <- r[r$geneid %in% dat_secretome$Ensembl.gene.id,]
  
  pdf(sprintf("out/out.de/secretome_%s.pdf",title), width = 15, height = 15)
  plotVolcano(r_sub, title, removesex)
  dev.off()
  write.csv(
    unrep(r_sub)[,c("geneid","baseMean","log2FoldChange","pvalue","padj","symbol")],
    sprintf("out/out.de/secretome_%s.csv",title),row.names = FALSE)
  
  ### Hya
  r_sub <- r[r$symbol %in% list_genes_hyaluronan_metabolic_process,]
  
  pdf(sprintf("out/out.de/hya_%s.pdf",title), width = 15, height = 15)
  plotVolcano(r_sub, title, removesex)
  dev.off()
  write.csv(
    unrep(r_sub)[,c("geneid","baseMean","log2FoldChange","pvalue","padj","symbol")],
    sprintf("out/out.de/hya_%s.csv",title),row.names = FALSE)
}

dat_secretome <- read.csv("secretome.csv",stringsAsFactors = FALSE)



####################################################
####################################################   2021-05-26: compare infected vs control, group high, no treatment, men and female sep
####################################################

keep <- ob_cond$inflevel2=="high" & ob_cond$Treatment %in% c("none") & ob_cond$sex=="male"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="2021-05-26 high inf vs mock men")



keep <- ob_cond$inflevel2=="high" & ob_cond$Treatment %in% c("none") & ob_cond$sex=="female"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="2021-05-26 high inf vs mock female")

####################################################
####################################################   2021-05-21: compare infected vs control, group high, no treatment treat
####################################################
#män vs kvinnor i gruppen High på de infekterade proverna - oinfekterade prover
#Dvs ((a5 + a6 minus a2), (k5 + k6 minus k2)) vs ((d5+d6 minus d2), (e5+e6 minus e2), (f5+f6 minus f2)).

# ob_cond$infected_woman <- "FALSE"
# ob_cond$infected_woman[ob_cond$infected=="TRUE" & ob_cond$sex=="female"] <- "TRUE"
# ob_cond$infected_man <- "FALSE"
# ob_cond$infected_man[ob_cond$infected=="TRUE" & ob_cond$sex=="male"] <- "TRUE"

# keep <- ob_cond$inflevel2=="high" & ob_cond$Treatment %in% c("none")
# dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected_woman+infected_man)
# dds <- DESeq(dds)
# resultsNames(dds)
# r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c(0,1,-1)))
# r$geneid <- rownames(r)
# r <- merge(r, ensid)
# r <- r[order(r$pvalue),]
# plotVolcano.save(r,title="2021-05-21 diff in infection men and woman compared to uninf in high")

ob_cond$inf_sex <- paste(ob_cond$infected, ob_cond$sex,sep="_")

keep <- ob_cond$inflevel2=="high" & ob_cond$Treatment %in% c("none")
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~inf_sex+0)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c(-1,1,1,-1)))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="2021-05-21 diff in infection men and woman compared to uninf in high")



####################################################
####################################################   2021-05-21: compare infected vs control, group high, no treatment treat
####################################################
#män vs kvinnor i gruppen High på de infekterade proverna - oinfekterade prover
#Dvs ((a5 + a6 minus a2), (k5 + k6 minus k2)) vs ((d5+d6 minus d2), (e5+e6 minus e2), (f5+f6 minus f2)).
# 
# ob_cond$name <- rownames(ob_cond)
# donors_group_high <- as.character(unique(ob_cond$donor[ob_cond$inflevel2=="high"]))
# thevec <- rep(0, nrow(ob_cond))
# thevec[ob_cond$inflevel2=="high" &  ob_cond$infected=="TRUE" & ob_cond$sex=="male"] <-  1/sum(ob_cond$inflevel2=="high" &  ob_cond$infected=="TRUE" & ob_cond$sex=="male")
# thevec[ob_cond$inflevel2=="high" &  ob_cond$infected=="FALSE" & ob_cond$sex=="male"] <- -1/sum(ob_cond$inflevel2=="high" & !ob_cond$infected=="FALSE" & ob_cond$sex=="male")
# thevec[ob_cond$inflevel2=="high" &  ob_cond$infected=="TRUE" & ob_cond$sex=="female"] <- -1/sum(ob_cond$inflevel2=="high" &  ob_cond$infected=="TRUE" & ob_cond$sex=="female")
# thevec[ob_cond$inflevel2=="high" &  ob_cond$infected=="FALSE" & ob_cond$sex=="female"] <-  1/sum(ob_cond$inflevel2=="high" & !ob_cond$infected=="FALSE" & ob_cond$sex=="female")
# 
# #keep <- ob_cond$donor %in% c(donors_group_high) & ob_cond$infected=="TRUE" & ob_cond$inflevel2=="high" & ob_cond$Treatment %in% c("none")
# dds <- DESeqDataSetFromMatrix(countData = ob_counts, colData = ob_cond, design = ~name+0)
# dds <- DESeq(dds)
# r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = thevec))
# r$geneid <- rownames(r)
# r <- merge(r, ensid)
# r <- r[order(r$pvalue),]
# plotVolcano.save(r,title="2021-05-21 diff in infection men and woman compared to uninf in high")


####################################################
####################################################   2021-05-13: compare men vs women, group high, infected samples
####################################################
#Dvs (a5 +a6, k5 + k6) vs (d5+d6, e5+e6, f5+f6).
#a,k,d,e,f

keep <- ob_cond$donor!="o" & ob_cond$infected=="TRUE" & ob_cond$inflevel2=="high" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~sex)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("sex","female","male")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="2021-05-13 infected men vs women in high")



####################################################
####################################################   beta, on infection
####################################################

############### Treatment effect: beta
keep <- ob_cond$donor=="o" & ob_cond$infected=="TRUE"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~Treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("Treatment","beta03","none")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
r_beta03 <- data.frame(
  symbol=r$symbol,
  fc_03=r$log2FoldChange,
  p_03=r$pvalue
)

plotVolcano.save(r,title="SUBSET tubes infected TEST beta 0.3nM vs none")


############### Treatment effect: beta -- somewhat dealing with infection
keep <- ob_cond$donor=="o" & ob_cond$infected=="TRUE"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~Treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("Treatment","beta3","none")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET tubes infected TEST beta 3nM vs none")

r_beta3 <- data.frame(
  symbol=r$symbol,
  fc_3=r$log2FoldChange,
  p_3=r$pvalue
)


############### Dosage dependence on beta concentration on infection
r_allbeta <- merge(unrep(r_beta03), unrep(r_beta3))
r_allbeta <- r_allbeta[log10(r_allbeta$p_03)< -2 | log10(r_allbeta$p_3)< -2,]
r_allbeta <- na.omit(r_allbeta)
pdf("out/out.de/comparison beta3 and beta0.3.pdf")
plot(r_allbeta$fc_03, r_allbeta$fc_3,cex=0)
text(r_allbeta$fc_03, r_allbeta$fc_3, labels = r_allbeta$symbol)
dev.off()



####################################################
####################################################   Impact of infection
####################################################

############### Infection effect in tubes
keep <- ob_cond$donor=="o" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET tubes untreated TEST infected vs uninfected")




############### Infection effect overall
keep <- ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~donorid+infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET tubesAND9ind untreated TEST infected vs uninfected")


############### Infection effect in 9ind
keep <- ob_cond$donor!="o" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected+donorid)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET 9ind untreated TEST infected vs uninfected")



############### Sex differences non-inf
keep <- ob_cond$donor!="o" & ob_cond$Treatment=="none" & ob_cond$infected=="FALSE"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~sex)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("sex","male","female")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET 9ind untreated uninfected TEST male vs female")


############### Sex differences - inf
keep <- ob_cond$donor!="o" & ob_cond$Treatment=="none" & ob_cond$infected=="TRUE"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~sex)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("sex","male","female")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET 9ind untreated infected TEST male vs female")



#####################
#####################  with 2 infection levels
#####################

keep <- ob_cond$donor!="o" & ob_cond$infected=="TRUE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~inflevel2)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("inflevel2","high","low")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r[r$baseMean>50,],title="inflevel_2 SUBSET 9ind infected untreated MODEL inflevel2 TEST infectionlevel(2) high vs low")


keep <- ob_cond$donor!="o" & ob_cond$infected=="FALSE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~inflevel2)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("inflevel2","high","low")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r[r$baseMean>50,],title="inflevel_2 SUBSET 9ind uninfected untreated MODEL inflevel2 TEST infectionlevel(2) high vs low")


keep <- ob_cond$donor!="o" & ob_cond$infected=="FALSE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~inflevel2)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("inflevel2","high","low")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r[r$baseMean>50,],title="inflevel_2 SUBSET 9ind uninfected untreated MODEL inflevel2 TEST infectionlevel(2) high vs low")


##############################

keep <- ob_cond$donor!="o" & ob_cond$donor!="c" & ob_cond$Treatment=="none"

ob_cond_red <- ob_cond[keep,]
ob_cond_red$donorid <- factor(as.character(ob_cond_red$donorid))
ob_cond_red$inflevel2 <- factor(as.character(ob_cond_red$inflevel2))
ob_cond_red$inflevel2_diff <- 0 + (ob_cond_red$infected=="TRUE" & ob_cond_red$inflevel2=="high")

dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond_red, design = ~donorid+infected+inflevel2_diff)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r[r$baseMean>50,],title="inflevel_2 SUBSET 9ind untreated MODEL donor infected inflevel2 TEST infectionlevel(2)_diff high vs low")

################################


#Uninfected high vs infected high
keep <- ob_cond$donor!="o" & ob_cond$inflevel2=="high" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r[r$baseMean>50,],title="inflevel_2 SUBSET 9ind high untreated MODEL inflevel2 TEST infected vs noninfected")


#Uninfected low vs infected low
keep <- ob_cond$donor!="o" & ob_cond$inflevel2=="low" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~infected)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("infected","TRUE","FALSE")))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r[r$baseMean>50,],title="inflevel_2 SUBSET 9ind low untreated MODEL inflevel2 TEST infected vs noninfected")


#####################
#####################  with 3 infection levels
#####################

combo_inf <- data.frame(
  lev1=c("high","high","mid"),
  lev2=c("mid","low","low"),
  stringsAsFactors = FALSE
)


keep <- ob_cond$donor!="o" & ob_cond$infected=="TRUE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~inflevel3)
dds <- DESeq(dds)
for(i in 1:3){
  r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("inflevel3",combo_inf$lev1[i],combo_inf$lev2[i])))
  r$geneid <- rownames(r)
  r <- merge(r, ensid)
  r <- r[order(r$pvalue),]
  plotVolcano.save(r[r$baseMean>50,],
                   title=sprintf("inflevel_3 SUBSET 9ind infected untreated TEST infectionlevel(3) %s vs %s",combo_inf$lev1[i],combo_inf$lev2[i]))
}


keep <- ob_cond$donor!="o" & ob_cond$infected=="FALSE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~inflevel3)
dds <- DESeq(dds)
for(i in 1:3){
  r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds, contrast = c("inflevel3",combo_inf$lev1[i],combo_inf$lev2[i])))
  r$geneid <- rownames(r)
  r <- merge(r, ensid)
  r <- r[order(r$pvalue),]
  plotVolcano.save(r[r$baseMean>50,],
                   title=sprintf("inflevel_3 SUBSET 9ind uninfected untreated TEST infectionlevel(3) %s vs %s",combo_inf$lev1[i],combo_inf$lev2[i]))
  
}




####################################################
####################################################   Differences that could affect infection severity
####################################################

keep <- ob_cond$donor!="o" & ob_cond$infected=="FALSE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~severity_pcr)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET 9ind noninf untreated CORRELATION severity_pcr")

keep <- ob_cond$donor!="o" & ob_cond$infected=="FALSE" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~severity_rnaseq)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="SUBSET 9ind noninf untreated CORRELATION severity_rnaseq")




####################################################
####################################################   Which genes correlate with infection level?
####################################################


############### qpcr covid measurement
keep <- ob_cond$donor!="o" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~donorid+log_newpcr) 
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="inflevel_cont SUBSET 9ind untreated CORRELATION covid_qpcr")


############### rnaseq covid measurement
keep <- ob_cond$donor!="o" & ob_cond$Treatment=="none"
dds <- DESeqDataSetFromMatrix(countData = ob_counts[,keep], colData = ob_cond[keep,], design = ~donorid+covidrnaseq)
dds <- DESeq(dds)
r <- as.data.frame(results(cooksCutoff=FALSE,independentFiltering=FALSE,dds))
r$geneid <- rownames(r)
r <- merge(r, ensid)
r <- r[order(r$pvalue),]
plotVolcano.save(r,title="inflevel_contSUBSET 9ind untreated CORRELATION covid_rnaseq")


### List of secreted genes that we stick to (from DE analysis with cutoff)
list_genes <- data.frame(symbol=c(
  'SERPINE1',
  'COL1A1',
  'FN1',
  'MMP10',
  'CES1',
  'VCAN',
  'LTBP1',
  'ENPP5',
  'LTBP2',
  'PXDN',
  'SERPINF1',
  'COL4A2',
  'SERPINE2',
  'PAMR1',
  'SPOCK1',
  'COL1A2',
  'AGR3',
  'APOD',
  'SERPINA1',
  'COL5A2',
  'MUC5AC'))


list_genes <- merge(list_genes,ensid_red)


###########################################################################
##################### Compare to CRISPR screen ############################
###########################################################################


##sup5   229E 
crispr_5 <- merge(ensid_red,read.csv("our_bulk/crisprscreen_tab5.csv"))
crispr_5 <-crispr_5[crispr_5$geneid %in% list_genes$geneid,]
crispr_5[order(crispr_5$p_value_neg),]
##serpinf1  rank 97



#https://www.nature.com/articles/s41588-021-00805-2
#Supplementary Table 10. SARS-CoV-2 low stringency screen - gene-level analysis SARS-CoV-2 versus control  
crispr_low <- merge(ensid_red,read.csv("our_bulk/crisprscreen_low.csv"))
crispr_low <- crispr_low[crispr_low$geneid %in% list_genes$geneid,]
crispr_low[order(crispr_low$p_value_neg),]
crispr_low[order(crispr_low$p_value_pos),]

#Supplementary Table 7. SARS-CoV-2 high stringency screen - gene-level analysis 
crispr_high <- merge(ensid_red,read.csv("our_bulk/crisprscreen_high.csv"))
crispr_high <-crispr_high[crispr_high$geneid %in% list_genes$geneid,]
crispr_high[order(crispr_high$p_value_neg),]

##sup7:
#serpinf1 rank 246ish
#pxdn rank 837
#






###########################################################################
##################### Compare to GWAS #####################################
###########################################################################

## Paper: https://www.nature.com/articles/s41586-020-03065-y

### Data here:
#GenoMICC EUR vs UK biobank controls (84.1Mb)
#https://datashare.is.ed.ac.uk/bitstream/handle/10283/3796/genomicc.EUR.PLINK2.txt.gz?sequence=1&isAllowed=y


##### Filter GWAS positions
gwas <- read.csv("~/temp/genomicc.EUR.PLINK2.txt",sep="\t")[,c("CHR","POS","Effect","Pval")]

subgwas <- gwas[gwas$Pval<0.001,c("CHR","POS","POS")]
subgwas$name <- 1:nrow(subgwas) 
write.table(subgwas,"~/temp/cord.bed", row.names = FALSE, col.names = FALSE,sep="\t")

homo <- read.csv("~/temp/homo.gtf",sep="\t")
homo <- homo[homo[,3]=="gene",]
homo <- homo[,c(1,4,5,9)]
write.table(homo, "~/temp/homo3.gtf", row.names = FALSE, col.names = FALSE, quote = FALSE,sep="\t")
colnames(homo)<-c("chr","from","to","desc")

d <- read.table("~/temp/closest_genename.txt",sep="\t")
d <- str_split_fixed(d$V8,"gene_id ",2)[,2]
gwas_ensid <- str_split_fixed(d,
  ";",2)[,1]

list_genes[list_genes$geneid %in% gwas_ensid,]

#https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/
#only SPOCK1 is near 1e-5
#adventitial fibroblast 


#### Some commands used as well:
#wget http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gff3.gz
#grep gene_id Homo_sapiens.GRCh38.97.gtf > homo.gtf
#awk -F '\t'  'BEGIN {OFS=","} { if (toupper($3) == "GENE")  print }' homo.gtf > homo2.gtf
#bedtools sort -i cord.bed > cord.bed.sorted
#bedtools sort -i homo3.gtf > homo3.gtf.sorted
#bedtools closest -a cord.bed.sorted -b homo3.gtf.sorted > closest.txt


## gwas 2
#https://www.nejm.org/doi/full/10.1056/NEJMoa2020283

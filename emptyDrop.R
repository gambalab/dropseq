library(DropletUtils)
library(Matrix)

###########################################################################################################
# Ribo Mito Genes (Code Snippet)
#
#
#library(biomaRt)
#ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
#all = biomaRt::getBM( attributes=c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), mart = ensembl)
#ix <- grepl(pattern = "^RPS", x = all$hgnc_symbol)
#ix <- ix | grepl(pattern = "^RPL", x = all$hgnc_symbol)
#ix <- ix | grepl(pattern = "^MT-", x = all$hgnc_symbol)
#ix <- ix | all$gene_biotype%in%c("rRNA","Mt_rRNA","Mt_tRNA")
#all = all[ix,]
#saveRDS(all,"~/dropseq_git/data/ribo.mito.genes.rds")
#
###########################################################################################################

rmdublecells = function(M) {cellSum <- Matrix::colSums(M);M[,cellSum < quantile(cellSum,.75) + 3*( quantile(cellSum,.75) - quantile(cellSum,.25))]}


# Input Parametres
# ----------------
args <- commandArgs(trailingOnly = TRUE)
wd = as.character(args[1]) #starSolo out
repoDir=as.character(args[2]) # git repo directory
sample=as.character(args[3]) # sample name
cpus=as.numeric(args[4]) # num. of cpus

ribomt = readRDS(paste0(repoDir,"/data/ribo.mito.genes.rds"))
dir.create(paste0(wd,"/RData"),showWarnings = F)

M = Matrix::readMM(paste0(wd,"/Gene/raw/matrix.mtx"))
bc = read.delim(paste0(wd,"/Gene/raw/barcodes.tsv"),stringsAsFactors = F,header = F)
colnames(bc) = "barcode"
colnames(M) = bc$barcode
genes = read.delim(paste0(wd,"/Gene/raw/features.tsv"),stringsAsFactors = F,header = F)
rownames(M) = substr(x = genes$V1,start = 1,stop = 15)

# Remove myco pseudo-genes
M = M[!grepl(pattern = "myco",x = rownames(M),fixed = T),]

# Run emptyDrops excluding mito/ribo
set.seed(100)
keep = !rownames(M)%in%ribomt$ensembl_gene_id
e.keep <- DropletUtils::emptyDrops(M[keep, ], retain = Inf,BPPARAM = BiocParallel::MulticoreParam(workers = cpus))
summary(e.keep$FDR < 0.01)

is.cell <- e.keep$FDR <= 0.01
#plot(e.keep$Total, -e.keep$LogProb, col=ifelse(is.cell, "red", "black"),xlab="Total UMI count", ylab="-Log Probability")

is.cell[is.na(is.cell)] = FALSE
M = M[,is.cell]

# Remove putative doublets
#M = rmdublecells(M)

# Remove low covarage cells
#M = M[,Matrix::colSums(M)>=500]

# Cut and Save results
M = M[Matrix::rowSums(M)>0,]
colnames(M) = paste(sample,colnames(M),sep = "_")
saveRDS(M,paste0(wd,"/RData/",sample,".Filtered.mtx.rds"))

        

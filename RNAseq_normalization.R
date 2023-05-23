setwd("/mnt/fs03/PROJECTS_MARIA/VIGex/ARTICLE/CODE/")

library(gdata)
library(openxlsx)

data <- read.table("raw_counts_RNA.txt", header=TRUE, sep="\t")


raw_counts <- data[,4:ncol(data)]

raw_bygene<-aggregate(raw_counts, list(data$gene_name) , sum)

rownames(raw_bygene)<-raw_bygene$Group.1
raw_bygene <- raw_bygene[,-c(1)]

library(DESeq2)

pheno_full<-data.frame(names=colnames(raw_bygene), dataset=rep("RNAv1", ncol(raw_bygene)))
pheno_full$names <- gsub("\\.","_",pheno_full$names, perl = TRUE)
rownames(pheno_full)<-pheno_full$names


dds <- DESeqDataSetFromMatrix(countData = raw_bygene, colData=pheno_full, checkDimnames=FALSE, design = ~ 1)

dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized=TRUE)

log2_norm <- log2(norm+1)





metadata<-openxlsx::read.xlsx("metadata_rna.xlsx")

c.id <- gsub("/","_",as.character(metadata$Capture.ID))
rownames(metadata)<-c.id

pheno_full$capture_id <- gsub("_counts","",pheno_full$names)
pheno_full$nano_id <- metadata[rownames(metadata),"VIGEX.ID"]












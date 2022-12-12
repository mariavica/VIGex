
library(openxlsx)

nanostring <- read.xlsx("Nanostring_data.xlsx", sheet="Raw_counts")

vigex_num <- t(nanostring[,-c(1)])

rownames(vigex_num)<-colnames(nanostring)[-c(1)]
colnames(vigex_num)<-nanostring$samples

positive <- c("POS_A","POS_B","POS_C","POS_D","POS_E","POS_F")
negative <- c("NEG_A","NEG_B","NEG_C","NEG_D","NEG_E","NEG_F","NEG_G","NEG_H")
housekeeping <- c("ACTB","ALAS1","CHMP2A.exon.1","CHMP2A.exon.3","CLTC","EMC7.exon.3","EMC7.exon.5","GAPDH.E1.E2.E3","GPI.exon.4",
                  "GPI.exon.6","MRPL19","OAZ1.E2.E3","POLR2A.E20.21","RPL19","RPLP0","SF3A1","TBP","TUBB","c1orf43.exon.2")
endo <-  c("BCL3","CCL3","CCL7","CD8A","CD24","CD44","CD68","CD163","CD274","CLCF1",
            "CSF1","CSF2RB","CTLA4","CXCL9","CXCL10","CXCL11","EBI3","FOXP3","GZMA",
            "GZMB","HLA.G","IER3","IFNg","IKZF2","IL6","IL7R","IL11","IL27","LIF",
            "LRRC32","MRC1","MSR1","NFKBIZ","OSM","PF4","PLAU","PLAUR","PRF1","SERPINE1",
            "SOCS3","STAM","TGFb1","TGFb2","TGFb3","TNFRSF12A","TNFRSF1B","TRIB1","VAV3",
            "VOPP1","MAFF","RBMS3","IL2","IL10","IL10RA","IL10RB","PDCD1","CCL5","CD27",
            "CD276","CMKLR1","CXCR6","HLA.DQA1","HLA.DRB1","HLA.E","IDO1","LAG3","NKG7",
            "PDCD1LG2","PSMB10","RAG1","RAG2","STAT1","TIGIT","TMEM173","TNFRSF9","CD4",
            "CD19","CD278","HAVCR2","HLA.A","HLA.B","HLA.C","HLA.DQB1")


pos.mean<-apply(vigex_num[positive,],2,geometric.mean)

mean.pos.mean<-mean(pos.mean)

norm.fact <- mean.pos.mean/pos.mean

vigex_num_pos <- array(0, dim=dim(vigex_num))
colnames(vigex_num_pos)<-colnames(vigex_num)
rownames(vigex_num_pos)<-rownames(vigex_num)


for ( i in 1:ncol(vigex_num_pos)) {
  vigex_num_pos[,i]<-vigex_num[,i]*norm.fact[i]
}


pos.mean.hk<-apply(vigex_num_pos[housekeeping,],2,geometric.mean)
mean.pos.mean.hk<-mean(pos.mean.hk)
norm.fact.hk <- mean.pos.mean.hk/pos.mean.hk
vigex_num_pos_hk <- vigex_num_pos


for ( i in 1:length(norm.fact.hk)) {
  vigex_num_pos_hk[c(housekeeping,endo),i]<-vigex_num_pos[c(housekeeping,endo),i]*norm.fact.hk[i]
  
  
}

library(factoextra)
endoss <- t(scale(t(data.frame(log2(vigex_num_pos_hk[endo,]))),center = TRUE, scale = FALSE))



library("impute")

endoss_imp_full <- impute.knn(endoss, k=50)$data






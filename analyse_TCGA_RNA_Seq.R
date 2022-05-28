setwd("~/Desktop/PhD_Project_related/TCGA ML PROJECT")


paad=read.table(file="TCGA_PAAD_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)
coad=read.table(file="TCGA_COAD_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)
read=read.table(file="TCGA_READ_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)
aml=read.table(file="TCGA_AML_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)

condition_df=as.data.frame(c(rep("PAAD",183),rep("COAD",521),rep("AML",255),rep("READ",177)))

all_counts_combined=cbind(paad,coad,aml,read)

colnames(condition_df)="condition"
condition_df$condition=as.factor(condition_df$condition)


library(DESeq2)
library(edgeR)
dds = DESeqDataSetFromMatrix(countData = all_counts_combined,
                             colData = condition_df,
                             design = ~ condition)

nrow(dds)
keep = rowSums(cpm(dds) > 5) >= 57 ## at least 5% samples with count of 1 or higher
dds = dds[keep,]
nrow(dds) ## 47609


vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)


rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)


*** PCA plot***
  
  plotPCA(vsd, intgroup = c("condition")) 

#### Running the differential expression pipeline
## we can run the differential expression pipeline on the raw counts with a single call to the function DESeq
dds <- DESeq(dds)

## Building the results table

Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula

res <- results(dds)
res


res_AML_READ <- results(dds, contrast=c("condition","AML","READ"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_RES, use.names = TRUE)
res_AML_READ_df=as.data.frame(res_AML_READ)
res_AML_READ_df$Ensembl=rownames(res_AML_READ_df)
summary(res_AML_READ)


res_COAD_READ <- results(dds, contrast=c("condition","COAD","READ"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_COAD_READ, use.names = TRUE)
res_COAD_READ_df=as.data.frame(res_COAD_READ)
res_COAD_READ_df$Ensembl=rownames(res_COAD_READ_df)
summary(res_COAD_READ)


res_PAAD_READ <- results(dds, contrast=c("condition","PAAD","READ"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_PAAD_READ, use.names = TRUE)
res_PAAD_READ_df=as.data.frame(res_PAAD_READ)
res_PAAD_READ_df$Ensembl=rownames(res_PAAD_READ_df)
summary(res_PAAD_READ)


res_AML_PAAD <- results(dds, contrast=c("condition","AML","PAAD"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_RES, use.names = TRUE)
res_AML_PAAD_df=as.data.frame(res_AML_PAAD)
res_AML_PAAD_df$Ensembl=rownames(res_AML_PAAD_df)
summary(res_AML_PAAD)


res_COAD_PAAD <- results(dds, contrast=c("condition","COAD","PAAD"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_COAD_PAAD, use.names = TRUE)
res_COAD_PAAD_df=as.data.frame(res_COAD_PAAD)
res_COAD_PAAD_df$Ensembl=rownames(res_COAD_PAAD_df)
summary(res_COAD_PAAD)


res_AML_COAD <- results(dds, contrast=c("condition","AML","COAD"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_RES, use.names = TRUE)
res_AML_COAD_df=as.data.frame(res_AML_COAD)
res_AML_COAD_df$Ensembl=rownames(res_AML_COAD_df)
summary(res_AML_COAD)


write.table(res_AML_READ_df,file="Diffexp_AML_READ.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(res_COAD_READ_df,file="Diffexp_COAD_READ.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(res_PAAD_READ_df,file="Diffexp_PAAD_READ.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(res_AML_PAAD_df,file="Diffexp_AML_PAAD.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(res_COAD_PAAD_df,file="Diffexp_COAD_PAAD.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(res_AML_COAD_df,file="Diffexp_AML_COAD.txt",col.names = T,row.names = F,sep="\t",quote = F)

norm.counts=counts(dds, normalized=T)

write.table(norm.counts,file="DESeq2_Normalised_Counts.txt",col.names = T,row.names = T,sep="\t",quote = F)


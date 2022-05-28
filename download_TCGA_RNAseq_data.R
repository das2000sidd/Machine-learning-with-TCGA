setwd("~/Desktop/PhD_Project_related/TCGA ML PROJECT")



library(TCGAbiolinks)
labels=read.csv(file="labels.csv",header = T,stringsAsFactors = F)
table(labels$Class) ## BRCA, COAD< KIRC LUAD PRAD

query <- GDCquery(project = "TARGET-AML",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

query2 <- GDCquery(project = "TCGA-READ",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "STAR - Counts")

query_COAD <- GDCquery(project = "TCGA-COAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "STAR - Counts")

query_KIRC <- GDCquery(project = "TCGA-KIRC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "STAR - Counts") ## 613

query_PAAD <- GDCquery(project = "TCGA-PAAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "STAR - Counts") ## 183



GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)

GDCdownload(
  query = query2, 
  method = "api", 
  files.per.chunk = 10
)

GDCdownload(
  query = query_COAD, 
  method = "client", 
  files.per.chunk = 10
)

GDCdownload(
  query = query_KIRC, 
  method = "client", 
  files.per.chunk = 10
) ## 613

GDCdownload(
  query = query_PAAD, 
  method = "client", 
  files.per.chunk = 10
) ## 183



data <- GDCprepare(query = query)
data2 <- GDCprepare(query = query2)
data_COAD <- GDCprepare(query = query_COAD)
#,directory = "/Users/siddhaduio.no/Desktop/PhD_Project_related/TCGA ML PROJECT/GDCdata/TCGA-COAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification")
data_PAAD<- GDCprepare(query = query_PAAD)
 

library(DT)
library(SummarizedExperiment)
datatable(
  assay(data)[1:20,], 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = TRUE
)

rowRanges(data)

exp_data=assay(data)
exp_data_2=assay(data2)
exp_data_COAD=assay(data_COAD)
exp_data_PAAD=assay(data_PAAD)


write.table(exp_data,file="TCGA_AML_raw_RNASeq_counts_data.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(exp_data_2,file="TCGA_READ_raw_RNASeq_counts_data.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(exp_data_COAD,file="TCGA_COAD_raw_RNASeq_counts_data.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(exp_data_PAAD,file="TCGA_PAAD_raw_RNASeq_counts_data.txt",col.names = T,row.names = T,sep="\t",quote = F)



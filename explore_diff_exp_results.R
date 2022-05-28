setwd("~/Desktop/PhD_Project_related/TCGA ML PROJECT")


r1=read.table(file="Diffexp_AML_COAD.txt",header = T,sep="\t",stringsAsFactors = F)
r2=read.table(file="Diffexp_COAD_PAAD.txt",header = T,sep="\t",stringsAsFactors = F)
r3=read.table(file="Diffexp_AML_PAAD.txt",header = T,sep="\t",stringsAsFactors = F)
r4=read.table(file="Diffexp_PAAD_READ.txt",header = T,sep="\t",stringsAsFactors = F)
r5=read.table(file="Diffexp_COAD_READ.txt",header = T,sep="\t",stringsAsFactors = F)
r6=read.table(file="Diffexp_AML_READ.txt",header = T,sep="\t",stringsAsFactors = F)

r1_diff=subset(r1,r1$padj < 0.01 & abs(r1$log2FoldChange) > 1)
r2_diff=subset(r2,r2$padj < 0.01 & abs(r2$log2FoldChange) > 1)
r3_diff=subset(r3,r3$padj < 0.01 & abs(r3$log2FoldChange) > 1)
r4_diff=subset(r4,r4$padj < 0.01 & abs(r4$log2FoldChange) > 1)
r5_diff=subset(r5,r5$padj < 0.01 & abs(r5$log2FoldChange) > 1)
r6_diff=subset(r6,r6$padj < 0.01 & abs(r6$log2FoldChange) > 1)

genes_of_interest=c(r1_diff$Ensembl,r2_diff$Ensembl,r3_diff$Ensembl,r4_diff$Ensembl,r5_diff$Ensembl,r6_diff$Ensembl)
genes_of_interest=unique(genes_of_interest)

all.norm.counts=read.table(file="DESeq2_Normalised_Counts.txt",header = T,sep="\t",stringsAsFactors = F)

norm.counts.goi=all.norm.counts[genes_of_interest,]

paad=read.table(file="TCGA_PAAD_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)
coad=read.table(file="TCGA_COAD_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)
read=read.table(file="TCGA_READ_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)
aml=read.table(file="TCGA_AML_raw_RNASeq_counts_data.txt",header = T,sep="\t",stringsAsFactors = F)


paad_samples=colnames(paad)
coad_samples=colnames(coad)
read_samples=colnames(read)
aml_samples=colnames(aml)

paad_samples=as.data.frame(paad_samples)
coad_samples=as.data.frame(coad_samples)
read_samples=as.data.frame(read_samples)
aml_samples=as.data.frame(aml_samples)

paad_samples$Class="PAAD"
coad_samples$Class="COAD"
read_samples$Class="READ"
aml_samples$Class="AML"

colnames(paad_samples)[1]="ID"
colnames(coad_samples)[1]="ID"
colnames(read_samples)[1]="ID"
colnames(aml_samples)[1]="ID"


all_samples_class=rbind(paad_samples,coad_samples,read_samples,aml_samples)
norm.counts.goi=t(norm.counts.goi)
rownames(all_samples_class)=all_samples_class$ID

class_order=all_samples_class$ID
norm.counts.goi.ordered=norm.counts.goi[class_order,]
gene_exp_order=rownames(norm.counts.goi.ordered)


check_order=0
for(index in 1:1136){
  if(gene_exp_order[index]==class_order[index]){
    check_order=check_order+1
  }
}

set.seed(96)
library(caret)
norm.counts.goi.ordered=as.matrix(norm.counts.goi.ordered)
norm.counts.goi.ordered=t(norm.counts.goi.ordered)
nzv <- nearZeroVar(norm.counts.goi.ordered, saveMetrics= TRUE)
## nothing is true in nzv

## idenityf highly correlated predictors
#descrCor <-  cor(norm.counts.goi)
#highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .999)


inTrain <- sample(1:1136, 568)
inTest = -inTrain

norm.counts.goi.ordered=t(norm.counts.goi.ordered)
training=norm.counts.goi.ordered[inTrain,]
test=norm.counts.goi.ordered[inTest,]
training_label=all_samples_class[inTrain,]
test_label=all_samples_class[inTest,]

table(training_label$Class)
table(test_label$Class)


fitControl = trainControl(method="cv",number = 5)
training_data_label=cbind(training_label,training[rownames(training_label),])
test_data_label=cbind(test_label,test[rownames(test_label),])

training_data_label$ID=NULL
test_data_label$ID=NULL


gbmFit1 <- train(Class ~ ., data = training_data_label, 
                 method = "gbm", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = TRUE)
gbmFit1


trellis.par.set(caretTheme())
plot(gbmFit1)  

plot(gbmFit1, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))


whichTwoPct <- tolerance(gbmFit1$results, metric = "Accuracy", 
                         tol = 2, maximize = TRUE) 

gbmFit1$results[whichTwoPct,]

gbmFit1.pred=predict(gbmFit1, newdata = test_data_label)
table(Pred=gbmFit1.pred,Actual=test_data_label$Class)
## As expected COAD and READ show worst performance

library(gbm)
varImp.gbm=varImp(gbmFit1)
varImp.gbm.tab=varImp.gbm$importance

write.table(varImp.gbm.tab,file="GBM_model_variable_importance.txt",col.names = T,row.names = T,sep="\t",quote = F)


rdaFit <- train(Class ~ ., data = training_data_label, 
                method = "rf", 
                trControl = fitControl, 
                tuneLength = 4,
                metric = "Accuracy")
rdaFit


library(randomForest)
library(mlbench)
library(e1071)


#Metric compare model is Accuracy
metric <- "Accuracy"
set.seed(123)
#Number randomely variable selected is mtry
mtry <- sqrt(ncol(training_data_label)-1)
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(Class~., 
                    data=training_data_label, 
                    method='rf', 
                    metric='Accuracy', 
                    tuneGrid=tunegrid, 
                    trControl=fitControl)


rf_default


trellis.par.set(caretTheme())
plot(rf_default)  

plot(rf_default, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))


whichTwoPct <- tolerance(rf_default$results, metric = "Accuracy", 
                         tol = 2, maximize = TRUE) 

rf_default$results

rf_default.pred=predict(rf_default, newdata = test_data_label)
table(Pred=rf_default.pred,Actual=test_data_label$Class)
## As expected COAD and READ show worst performance

varImp.rf=varImp(rf_default)
varImp.rf.tab=varImp.rf$importance

write.table(varImp.rf.tab,file="RF_model_variable_importance.txt",col.names = T,row.names = T,sep="\t",quote = F)


varImpPlot(rf_default,sort=TRUE,n.var = 30)

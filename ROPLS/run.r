library(ropls)
getwd()
## [1] "/media/sda/user/chengxu/project/SIKE-NASH/monkey_RNA-Seq/PLS-DA_20190604"

dataMatrix0 <- read.table("NAS_matrix.txt",header=T,sep="\t",check.names=F,row.names=1)
dataMatrix0 <- t(dataMatrix0)
sampleMetadata0 <- read.table("NAS_metadata.txt",header=F,row.names=1)
colnames(sampleMetadata0) <- c("group","NAS","fibrosis")

keep <- apply(dataMatrix0,2,mean)>0.5
summary(keep)
dataMatrix0 <- dataMatrix0[,keep]

###########################################################

### Expression
s1 <- grep("NAS_2",rownames(dataMatrix0))
s2 <- grep("NAS_3_4",rownames(dataMatrix0))
dataMatrix <- dataMatrix0[c(s1,s2),]

### NAS_2 vs. NAS_3_4
s1 <- grep("NAS_2",rownames(sampleMetadata0))
s2 <- grep("NAS_3_4",rownames(sampleMetadata0))
sampleMetadata <- sampleMetadata0[c(s1,s2),]
group <- sampleMetadata[,"group"]
NAS <- sampleMetadata[,"NAS"]

### PCC between expr and NAS
pvaVn <- apply(dataMatrix,2,function(feaVn) cor.test(NAS,feaVn)[["p.value"]])
corVn <- apply(dataMatrix,2,function(feaVn) cor.test(NAS,feaVn)[["estimate"]])
dataMatrix <- dataMatrix[,pvaVn<0.05 & !is.na(pvaVn)]
dir.create("ropls_NAS_2_NAS_3_4")
out <- cbind(corVn,pvaVn)
out <- cbind(rownames(out),out)
colnames(out)[1] <- "gene"
write.table(out,file="ropls_NAS_2_NAS_3_4/PCC_NAS_output.xls",sep="\t",quote=F,row.names=F)

### PCA
pca <- opls(dataMatrix,plotL=FALSE)
pdf("ropls_NAS_2_NAS_3_4/PCA_summary_plot.pdf",width=8,height=8)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for(typeC in c("overview","outlier","x-score","x-loading"))
plot(pca,typeVc=typeC,parDevNewL=FALSE)
dev.off()
pdf("ropls_NAS_2_NAS_3_4/PCA_score_plot.pdf",width=8,height=8)
plot(pca,typeVc="x-score",parAsColFcVn=group,parEllipsesL=TRUE,parDevNewL=FALSE)
dev.off()

### PLS-DA
plsda <- opls(dataMatrix,as.character(group),plotL=FALSE)
pdf("ropls_NAS_2_NAS_3_4/PLS-DA_summary_plot.pdf",width=8,height=8)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for(typeC in c("permutation","overview","outlier","x-score"))
plot(plsda,typeVc=typeC,parDevNewL=FALSE)
dev.off()
pdf("ropls_NAS_2_NAS_3_4/PLS-DA_score_plot.pdf",width=8,height=8)
plot(plsda,typeVc="x-score",parAsColFcVn=group,parEllipsesL=TRUE,parDevNewL=FALSE) ## @scoreMN
dev.off()
write(names(plsda@vipVn[plsda@vipVn>1]),file="ropls_NAS_2_NAS_3_4/PLS-DA_vip1.list")

### OPLS-DA
oplsda <- opls(dataMatrix,as.character(group),predI=1,orthoI=NA,plotL=FALSE)
pdf("ropls_NAS_2_NAS_3_4/OPLS-DA_summary_plot.pdf",width=8,height=8)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for(typeC in c("permutation","overview","outlier","x-score"))
plot(oplsda,typeVc=typeC,parDevNewL=FALSE)
dev.off()
pdf("ropls_NAS_2_NAS_3_4/OPLS-DA_score_plot.pdf",width=8,height=8)
plot(oplsda,typeVc="x-score",parAsColFcVn=group,parEllipsesL=TRUE,parDevNewL=FALSE) ## @scoreMN
dev.off()
write(names(oplsda@vipVn[oplsda@vipVn>1]),file="ropls_NAS_2_NAS_3_4/OPLS-DA_vip1.list")

###########################################################

### Expression
s1 <- grep("NAS_3_4",rownames(dataMatrix0))
s2 <- grep("NAS_5",rownames(dataMatrix0))
dataMatrix <- dataMatrix0[c(s1,s2),]

### NAS_3_4 vs. NAS_5
s1 <- grep("NAS_3_4",rownames(sampleMetadata0))
s2 <- grep("NAS_5",rownames(sampleMetadata0))
sampleMetadata <- sampleMetadata0[c(s1,s2),]
group <- sampleMetadata[,"group"]
NAS <- sampleMetadata[,"NAS"]

### PCC between expr and NAS
pvaVn <- apply(dataMatrix,2,function(feaVn) cor.test(NAS,feaVn)[["p.value"]])
corVn <- apply(dataMatrix,2,function(feaVn) cor.test(NAS,feaVn)[["estimate"]])
dataMatrix <- dataMatrix[,pvaVn<0.05 & !is.na(pvaVn)]
dir.create("ropls_NAS_3_4_NAS_5")
out <- cbind(corVn,pvaVn)
out <- cbind(rownames(out),out)
colnames(out)[1] <- "gene"
write.table(out,file="ropls_NAS_3_4_NAS_5/PCC_NAS_output.xls",sep="\t",quote=F,row.names=F)

### PCA
pca <- opls(dataMatrix,plotL=FALSE)
pdf("ropls_NAS_3_4_NAS_5/PCA_summary_plot.pdf",width=8,height=8)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for(typeC in c("overview","outlier","x-score","x-loading"))
plot(pca,typeVc=typeC,parDevNewL=FALSE)
dev.off()
pdf("ropls_NAS_3_4_NAS_5/PCA_score_plot.pdf",width=8,height=8)
plot(pca,typeVc="x-score",parAsColFcVn=group,parEllipsesL=TRUE,parDevNewL=FALSE)
dev.off()

### PLS-DA
plsda <- opls(dataMatrix,as.character(group),plotL=FALSE)
pdf("ropls_NAS_3_4_NAS_5/PLS-DA_summary_plot.pdf",width=8,height=8)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for(typeC in c("permutation","overview","outlier","x-score"))
plot(plsda,typeVc=typeC,parDevNewL=FALSE)
dev.off()
pdf("ropls_NAS_3_4_NAS_5/PLS-DA_score_plot.pdf",width=8,height=8)
plot(plsda,typeVc="x-score",parAsColFcVn=group,parEllipsesL=TRUE,parDevNewL=FALSE) ## @scoreMN
dev.off()
write(names(plsda@vipVn[plsda@vipVn>1]),file="ropls_NAS_3_4_NAS_5/PLS-DA_vip1.list")

### OPLS-DA
oplsda <- opls(dataMatrix,as.character(group),predI=1,orthoI=NA,plotL=FALSE)
pdf("ropls_NAS_3_4_NAS_5/OPLS-DA_summary_plot.pdf",width=8,height=8)
layout(matrix(1:4,nrow=2,byrow=TRUE))
for(typeC in c("permutation","overview","outlier","x-score"))
plot(oplsda,typeVc=typeC,parDevNewL=FALSE)
dev.off()
pdf("ropls_NAS_3_4_NAS_5/OPLS-DA_score_plot.pdf",width=8,height=8)
plot(oplsda,typeVc="x-score",parAsColFcVn=group,parEllipsesL=TRUE,parDevNewL=FALSE) ## @scoreMN
dev.off()
write(names(oplsda@vipVn[oplsda@vipVn>1]),file="ropls_NAS_3_4_NAS_5/OPLS-DA_vip1.list")

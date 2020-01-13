args  <- commandArgs(TRUE)

genesetgmt <- args[1]
geneexp <- args[2]
control <- args[3]
treat   <- args[4]
options(stringsAsFactors=F)

#载入R包
library(GSEABase)
library(GSVAdata)
library(GSVA)
library(limma)


# genesetgmt <- "batch1_Rep-6h_DEGs_up.gmt"
# geneexp <- "@protein_coding.symbol.max_fpkm.del@2+.xls"
# control <- "NC"
# treat <- "HFD-12w"


#读入基因集和表达量文件
geneSets <- getGmt(genesetgmt)
mydata <- read.table(file = geneexp,header=T,check.names=F)

#对表达量文件根据第一列进行去重
index <- duplicated(mydata[,1])
fildup=mydata[!index,]
exp=fildup[,-1]
name=fildup[,1]
row.names(exp)=name

#由于样品名称可能是包含小横线,所以作为参数输入的时候用下划线代替,在这里会替换成小横线,如果命名不包含小横线，这一行命令可以删除
control=sub(pattern = "_", replacement = "-", control)
treat=sub(pattern = "_", replacement = "-", treat)


#从多个样品中提取需要进行差异分析的比较组合生成一个新的矩阵
condition1 <- paste("^",control,"-[0-9]+$",sep="")
condition2 <- paste("^",treat,"-[0-9]+$",sep="")
sample <- colnames(exp)
colnum1=grep(condition1,sample)
colnum2=grep(condition2,sample)
sssname=sample[c(colnum1,colnum2)]
diffexp=exp[,sssname]

#计算比较组合的样品数目
cnum=length(colnum1)
print(paste("control num is:",cnum,sep=""))
tnum=length(colnum2)
print(paste("treat num is:",tnum,sep=""))


#将数据框转换成矩阵
mydata= as.matrix(diffexp)
#对表达量值进行均一化
mydata <- log2(mydata+1)
#使用GSVA得到gsva 得分矩阵,mx.diff=TRUE是得到一个近似的正态分布，如果使用limma来进行下游的差异分析，这里可以设置为默认的TRUE
es.diff <- gsva(mydata,geneSets,mx.diff=TRUE,min.sz=10,max.sz=500,parallel.sz=4)

#修改GSVA矩阵的名称并输出矩阵结果
GSVA_score <- cbind(rownames(es.diff),es.diff)
colnames(GSVA_score)[1] <- "path"
write.table(GSVA_score,file=paste(control,"_vs_",treat,".gsva_score.xls",sep=""),quote=F,sep="\t",row.names=F)
print("GSVA score matrix was created successfully......")


#构建一个差异分组的矩阵
# group <- c(rep(control,cnum),rep(treat,tnum))
# phno <- data.frame(sssname,group)
# group0 <- factor(phno$group,levels=levels(phno$group))
# design <- model.matrix(~0+group0)
# colnames(design) <- c("control","treat")


#构建实验设计矩阵，这里根据实际的情况设置（表型）分组，与表达矩阵的列对应起来
group_list <- c(rep("control",cnum),rep("treat",tnum))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))


#获取用于gsva得分的矩阵，进行差异分析
diff=es.diff[,sssname]
rownames(design)=colnames(diff)
contrast.matrix<-makeContrasts(treat-control,levels = design)
print("difference grouping matrix was created successfully......")


#进行拟合和差异分析
#线性模型拟合
fit <- lmFit(diff,design)
#根据对比模型进行差值计算
fit <- contrasts.fit(fit, contrast.matrix)
#贝叶斯检验
fit <- eBayes(fit)



#生成所有基因集的检验结果报告
allGenes <- topTable(fit,coef=1,number=Inf)

DEs <- cbind(rownames(allGenes),allGenes)
DEs_output <- merge(GSVA_score,DEs,by.x="path",by.y="rownames(allGenes)")
DEs_output <- DEs_output[order(DEs_output$P.Value),]
write.table(DEs_output,file=paste("DEs_",control,"_vs_",treat,".limma_output.xls",sep=""),quote=F,sep="\t",row.names=F)

#
#用P.Value进行筛选，得到全部差异表达基因
DEs_output_filterbyp <- DEs_output[DEs_output[, "P.Value"]<0.01,]
write.table(DEs_output_filterbyp,file=paste("DEs_",control,"_vs_",treat,".limma_output.filter.xls",sep=""),quote=F,sep="\t",row.names=F)
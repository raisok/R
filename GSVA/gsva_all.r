args  <- commandArgs(TRUE)

genesetgmt <- args[1]
geneexp <- args[2]
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
mydata <- read.table(file = geneexp,header=T,check.names=F,sep="\t")

#对表达量文件根据第一列进行去重
index <- duplicated(mydata[,1])
fildup=mydata[!index,]
exp=fildup[,-1]
exp = subset(exp,select=-c(description))
name=fildup[,1]
row.names(exp)=name

#对所有的表达量矩阵进行gsva计算
#将数据框转换成矩阵
allmydata= as.matrix(exp)
#对表达量值进行均一化
allmydata <- log2(allmydata+1)
#使用GSVA得到gsva 得分矩阵
es.diff <- gsva(allmydata,geneSets,mx.diff=TRUE,min.sz=10,max.sz=500,parallel.sz=4)

#修改GSVA矩阵的名称并输出矩阵结果
all_GSVA_score <- cbind(rownames(es.diff),es.diff)
colnames(all_GSVA_score)[1] <- "path"
write.table(all_GSVA_score,file="all.gsva_score.xls",,quote=F,sep="\t",row.names=F)
print("all GSVA score matrix was created successfully......")


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
library(pheatmap)


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



#由于样品名称可能是包含小横线,所以用下划线代替,再替换成小横线,如果命名规范可以不要这一行命令
#treat=sub(pattern = "_", replacement = "-", treat)
#control=sub(pattern = "_", replacement = "-", control)


#提取比较组合生成一个新的矩阵
condition1 <- paste("^",control,"-[0-9]+$",sep="")
condition2 <- paste("^",treat,"-[0-9]+$",sep="")
print(condition2)
sample <- colnames(exp)
print(sample)
colnum1=grep(condition1,sample)
colnum2=grep(condition2,sample)
sssname=sample[c(colnum1,colnum2)]
diffexp=exp[,sssname]

#计算组合的样品数目
cnum=length(colnum1)
tnum=length(colnum2)

print(paste("control num is:",cnum,sep=""))
print(paste("treat num is:",tnum,sep=""))


#将数据框转换成矩阵
mydata= as.matrix(diffexp)
#对表达量值进行均一化
mydata <- log2(mydata+1)
#使用GSVA得到gsva 得分矩阵
es.diff <- gsva(mydata,geneSets,mx.diff=TRUE,min.sz=10,max.sz=500,parallel.sz=4)

#修改GSVA矩阵的名称并输出矩阵结果
GSVA_score <- cbind(rownames(es.diff),es.diff)
colnames(GSVA_score)[1] <- "path"
write.table(GSVA_score,file=paste(control,"_vs_",treat,".gsva_score.xls",sep=""),quote=F,sep="\t",row.names=F)
print("GSVA score matrix was created successfully......")


#构建差异分组矩阵
group_list <- c(rep("control",cnum),rep("treat",tnum))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))


#获取用于gsva得分的矩阵，进行差异分析
diff=es.diff[,sssname]
rownames(design)=colnames(diff)
contrast.matrix<-makeContrasts(treat-control,levels = design)
print("difference grouping matrix was created successfully......")


#进行拟合和差异分析
fit <- lmFit(diff,design)

fit <- contrasts.fit(fit, contrast.matrix)

fit <- eBayes(fit)


#提取差异分析的结果并输出
allGenes <- topTable(fit,coef=1,number=Inf)
DEs <- cbind(rownames(allGenes),allGenes)

DEs_output <- merge(GSVA_score,DEs,by.x="path",by.y="rownames(allGenes)")
DEs_output <- DEs_output[order(DEs_output$P.Value),]
write.table(DEs_output,file=paste("DEs_",control,"_vs_",treat,".limma_output.xls",sep=""),quote=F,sep="\t",row.names=F)


#用P.Value进行筛选，得到全部差异表达基因
DEs_output_filterbyp <- DEs_output[DEs_output[, "P.Value"]<0.01,]
write.table(DEs_output_filterbyp,file=paste("DEs_",control,"_vs_",treat,".limma_output.filter.xls",sep=""),quote=F,sep="\t",row.names=F)


data_heatmap = read.table(paste("DEs_",control,"_vs_",treat,".limma_output.xls",sep=""),sep="\t",header=T,check.names=F)

pre_heatmap = data_heatmap[,sssname]
rownames(pre_heatmap) = data_heatmap$path

sam <- gsub("-","_",gsub("-[^-]+$","",colnames(pre_heatmap)))

anno_col <- data.frame(sample=factor(sam))
rownames(anno_col) <- colnames(pre_heatmap)
color <- c("green","yellow","pink")
anno_colors <- NA
sam_colors <- NA
for(i in 1:length(unique(sam))) {
		  set0 <- paste(unique(sam)[i],"='",color[i],"'",sep="")
  sam_colors <- c(sam_colors,set0)
}
set_sam_colors <- paste("sample=c(",paste(sam_colors[-1],collapse=","),")",sep="")
anno_colors <- c(anno_colors,set_sam_colors)
text0 <- paste("anno_colors <- list(",paste(anno_colors[-1],collapse=","),")",sep="")
eval(parse(text=text0))

pre_heatmap[pre_heatmap>=1] <- 1
pre_heatmap[pre_heatmap<=-1] <- -1

if(nrow(pre_heatmap) >1){
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
         cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
         annotation_col=anno_col,
         annotation_colors=anno_colors,
         border_color=NA,
	     color=colorRampPalette(c("blue","white","red"))(100),
         filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.pdf",sep=""))
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
			 cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
			 annotation_col=anno_col,
			 annotation_colors=anno_colors,
			 border_color=NA,
			 color=colorRampPalette(c("blue","white","red"))(100),
			 filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.png",sep=""))
}else{
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
             cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
	         annotation_col=anno_col,
	         annotation_colors=anno_colors,
             border_color=NA,
             color=colorRampPalette(c("blue","white","red"))(100),
	         filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.pdf",sep=""))
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
			 cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
		     annotation_col=anno_col,
 			 annotation_colors=anno_colors,
  			 border_color=NA,
 			 color=colorRampPalette(c("blue","white","red"))(100),
 			 filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.png",sep=""))
}
print("draw heatmap done!")






data_heatmap = read.table(paste("DEs_",control,"_vs_",treat,".limma_output.filter.xls",sep=""),sep="\t",header=T,check.names=F)

pre_heatmap = data_heatmap[,sssname]
rownames(pre_heatmap) = data_heatmap$path

sam <- gsub("-","_",gsub("-[^-]+$","",colnames(pre_heatmap)))

anno_col <- data.frame(sample=factor(sam))
rownames(anno_col) <- colnames(pre_heatmap)
color <- c("green","yellow","pink")
anno_colors <- NA
sam_colors <- NA
for(i in 1:length(unique(sam))) {
		  set0 <- paste(unique(sam)[i],"='",color[i],"'",sep="")
  sam_colors <- c(sam_colors,set0)
}
set_sam_colors <- paste("sample=c(",paste(sam_colors[-1],collapse=","),")",sep="")
anno_colors <- c(anno_colors,set_sam_colors)
text0 <- paste("anno_colors <- list(",paste(anno_colors[-1],collapse=","),")",sep="")
eval(parse(text=text0))

pre_heatmap[pre_heatmap>=1] <- 1
pre_heatmap[pre_heatmap<=-1] <- -1

if(nrow(pre_heatmap) >1){
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
         cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
         annotation_col=anno_col,
         annotation_colors=anno_colors,
         border_color=NA,
	     color=colorRampPalette(c("blue","white","red"))(100),
         filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.pdf",sep=""))
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
			 cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
			 annotation_col=anno_col,
			 annotation_colors=anno_colors,
			 border_color=NA,
			 color=colorRampPalette(c("blue","white","red"))(100),
			 filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.png",sep=""))
}else{
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
             cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
	         annotation_col=anno_col,
	         annotation_colors=anno_colors,
             border_color=NA,
             color=colorRampPalette(c("blue","white","red"))(100),
	         filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.pdf",sep=""))
    pheatmap(pre_heatmap,cluster_row=F,cluster_cols=FALSE,
			 cellwidth=12,cellheight=9,fontsize=6,fontsize_row=6,
		     annotation_col=anno_col,
 			 annotation_colors=anno_colors,
  			 border_color=NA,
 			 color=colorRampPalette(c("blue","white","red"))(100),
 			 filename=paste("DEs_",control,"_vs_",treat,".limma_output.filter.png",sep=""))
}
print("draw heatmap done!")

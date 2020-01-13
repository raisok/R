library(pheatmap)
library(RColorBrewer)

args = commandArgs()

fpkm <- read.csv(args[5],header=T,sep="\t",check.names=F)
reactome_list <- read.csv(args[6],header=T,sep="\t",check.names=F)
if(nrow(reactome_list)>=50){
reactome_list <- reactome_list[1:50,]
}else if (nrow(reactome_list)==0){
cat(c("Warning",args[6],"is a blank file","\n"))
q()
}

k <- nrow(reactome_list)
for(i in 1:k){
symbol <- reactome_list[i,]$geneID
symbol <- strsplit(as.character(symbol),"/")
symbol <- unlist(symbol)
expr1 <- fpkm[match(symbol,fpkm$gene),]

colnum1 <- grep(paste("^",args[7],"-[0-9]+$",sep=""),colnames(expr1))
control_num = length(colnum1)
colnum2 <- grep(paste("^",args[8],"-[0-9]+$",sep=""),colnames(expr1))
treat_num = length(colnum2)
expr2 <- expr1[,c(colnum1,colnum2)]
rownames(expr2) <- expr1$gene

name <- as.character(reactome_list[i,]$Description)
name1 <- gsub("/|'|\"","",name)
name1 <- gsub(" ","_",name1)
name1 <- gsub("_-_","-",name1)

#### get control and treat mean, then sort by fold
expr2$controlmean <- apply(expr2[,1:control_num], 1, mean,na.rm=T)
expr2$treatmean <- apply(expr2[,(control_num+1):(control_num+treat_num)], 1, mean,na.rm=T)
expr2$fold <- expr2$treatmean/expr2$controlmean
expr2 = expr2[order(expr2$fold,decreasing = TRUE),]
expr2 = subset(expr2,select = -c(controlmean,treatmean,fold))
#### add by yueyao at 20191121

write.csv(expr2,file=paste(name1,".csv",sep=""),quote=F,row.names=T)

expr2 = t(scale(t(expr2)))
expr2[expr2>=1] <- 1
expr2[expr2<=-1] <- -1
if(nrow(expr2)>1){
pheatmap(expr2,cluster_row=FALSE,cluster_cols=FALSE,
         cellwidth=15,cellheight=6,
         fontsize=6,fontsize_row=6,border_color = NA,
         gaps_col = c(control_num),drop_levels = TRUE,
         fontsize_col=6,
         #border_color="black",
         color=colorRampPalette(c("blue","white","red"))(100),
         main=name,filename=paste(name1,".pdf",sep=""))

pheatmap(expr2,cluster_row=FALSE,cluster_cols=FALSE,
         cellwidth=15,cellheight=6,
         fontsize=6,fontsize_row=6,fontsize_col=6,border_color = NA,
         gaps_col = c(control_num),drop_levels = TRUE,
         #border_color="black",
         color=colorRampPalette(c("blue","white","red"))(100),
         main=name,filename=paste(name1,".png",sep=""))
} else {
pheatmap(expr2,cluster_row=FALSE,cluster_cols=FALSE,
         cellwidth=15,cellheight=6,
         fontsize=6,fontsize_row=6,fontsize_col=6,border_color = NA,
         drop_levels = TRUE,
         #border_color="black",
         color=colorRampPalette(c("blue","white","red"))(100),
         main=name,filename=paste(name1,".pdf",sep=""))

pheatmap(expr2,cluster_row=FALSE,cluster_cols=FALSE,
         cellwidth=15,cellheight=6,
         fontsize=6,fontsize_row=6,fontsize_col=6,border_color = NA,
         #border_color="black",
         color=colorRampPalette(c("blue","white","red"))(100),
         main=name,filename=paste(name1,".png",sep=""))
}
}
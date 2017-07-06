#利用循环同时对多个文件做折线图的R脚本
#gene1 = log10(gene)
#gene1[gene1 == "-Inf"] =0
#gene1[gene1 >= 20] =20
setwd("c:/Users/yueyao/Desktop/latex_cluster/")
for(j in 1:10){
input=paste("cluster",j,".fpkm.txt",sep="")
output=paste("cluster",j,".pdf",sep="")
gene1 = read.table(file=input,header = T,row.names = 1)
x=c(1:27)
y= c(gene1[1,])
pdf(output)
plot(x,y,type = "l",xlab = "Species", ylab="fpkm")
len = nrow(gene1)
for(i in 2:len){
a=c(gene1[i,])
lines(x,a)
}
dev.off()
}
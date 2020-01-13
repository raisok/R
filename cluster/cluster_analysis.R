library(gmodels)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(ggrepel)

args <- commandArgs()

#### Only for Test
setwd("E:/IMA/Rscript/cluster/")
args[5] = "Ctrl_vs_XM1845_OE.fpkm.txt"
args[6] = "Ctrl"
args[7] = "XM1845_OE"
#windows set font
#windowsFonts(Times = windowsFont("Times New Roman"))
#windowsFonts(Arial = windowsFont("Arial"))
###

expr <- read.table(args[5],sep="\t",header=T,check.names=F,quote="")
if(length(args)==7){

colnum1 <- grep(paste("^",args[6],"-[0-9]+$",sep=""),colnames(expr))
colnum2 <- grep(paste("^",args[7],"-[0-9]+$",sep=""),colnames(expr))
expr <- expr[,c(1,colnum1,colnum2)]
keep <- apply(expr[,-1],1,sum) > 0
write.table(expr[keep,],file=paste(args[6],"_vs_",args[7],".fpkm.txt",sep=""),sep="\t",quote=F,row.names=F)
expr <- expr[,-1]
} else {
expr <- expr[,c(-1,-2)]
}

group <- rep(0,ncol(expr))
sorted <- rep(0,(length(args)-5))
for(m in 6:length(args)) {
	M = m-5
	sorted[M] = args[m]
	colnum = grep(paste("^",args[m],sep=""),colnames(expr))
	for(n in 1:length(colnum)) {
	group[colnum[n]] = args[m]
	}
}

data = t(as.matrix(log2(expr+1)))

data.pca = fast.prcomp(data,retx=T,scale=F,center=T)

a = summary(data.pca)
tmp = a[4]$importance
pro1 = as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2 = as.numeric(sprintf("%.3f",tmp[2,2]))*100

pc = as.data.frame(a$x)
pc$group = group
pc$names = rownames(pc)

xlab = paste("PC1(",pro1,"%)",sep="")
ylab = paste("PC2(",pro2,"%)",sep="")
pca_nolab = ggplot(pc,aes(PC1,PC2)) +geom_point(size=2,aes(color=group)) +
  labs(x=xlab,y=ylab) +geom_hline(yintercept=0,linetype="dashed",color="black") +
  geom_vline(xintercept=0,linetype="dashed",color="black") +theme_bw()+
  theme(axis.title.x=element_text(colour = "black",size=6),
        axis.title.y=element_text(colour = "black",size=6),
        panel.grid=element_blank(),
        axis.line=element_line(size=0.25,colour="black"),
        axis.text = element_text(colour = "black",size=6)
        ) +
  scale_color_discrete(limits=sorted,guide=FALSE) +labs(color="Sample") +
  stat_ellipse(geom="polygon",alpha=1/4,aes(fill=group)) +
  scale_fill_discrete(limits=sorted,guide=FALSE) +
  geom_rug(aes(color=group))

#### change by yueyao at 2019.11.19
pca = pca_nolab+geom_text_repel(aes(label=names),size=2) 
#### change by yueyao at 2019.11.19

if(length(args)==7){
out1 = paste(sorted,collapse="_vs_")
} else {
out1 = "allsamples"
}
out2 = paste("PCA.",out1,sep="")
ggsave(paste(out2,".pdf",sep=""),pca,width=2,height=2)
ggsave(paste(out2,".png",sep=""),pca,width=2,height=2)

####change by yueyao at 2019.11.19
ggsave(paste(out2,".nolab.pdf",sep=""),pca_nolab,width=2,height=2)
ggsave(paste(out2,".nolab.png",sep=""),pca_nolab,width=2,height=2)
####change by yueyao at 2019.11.19

###################################################################################################


hc = hclust(dist(data),method="average")

color = rep(0,length(labels(hc)))
for(m in 6:length(args)) {
	M = m-5
	colnum = grep(paste("^",args[m],"-[0-9]+",sep=""),labels(hc))
	for(n in 1:length(colnum)) {
	color[colnum[n]] = rainbow(length(args)-5)[M]
	}
}

dend = as.dendrogram(hc) %>% set("labels_col",color) %>% set("labels_cex",1)  %>% set("leaves_pch",19) %>% set("leaves_cex",1) %>% set("leaves_col",color)

out3 = paste("hclust",".",out1,sep="")
rw = length(labels(hc))/2.5
pdf(file=paste(out3,".pdf",sep=""),width=rw,height=3)
par(mar=c(10,3,2,0),cex=1)
plot(dend,ylab="Height",cex=0.5)
dev.off()

rw = rw*80
png(filename=paste(out3,".png",sep=""),width=rw,height=300)
par(mar=c(10,3,2,0),cex=1)
plot(dend,ylab="Height",cex=0.5)
dev.off()

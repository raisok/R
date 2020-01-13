library(pheatmap)




args<-commandArgs(T)
ff=args[1]
pheatout=args[2]
corout=args[3]
fpkm=read.table(file =ff,header = T,row.names = 1,sep = "\t")
draw=fpkm[,!names(fpkm) %in% c("description","total_count","baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
data_plot<-t(scale(t(draw)))
data_plot[data_plot>1.5]<-1.5
data_plot[data_plot<(-1.5)]<-(-1.5)

pheatmap(data_plot,filename=pheatout)

cordraw=cor(draw)

pheatmap(cordraw,filename=corout)


library(reshape2)
library(ggplot2)
library(RColorBrewer)

xx=as.matrix(cordraw)
xx=melt(xx)
names(xx)=c("Var1","Var2","pearson_value");
pdf(corout,width=8,height=8)
ggplot(xx, aes(Var1, Var2, fill=pearson_value))+
  geom_tile(color='black')+
  geom_text(label=round(xx$pearson_value, 3))+
  scale_fill_gradient(low='#DEEBF7',high='#08519C')+
  theme(axis.text = element_text(angle=30, hjust=1,size=11,vjust=0,color='black'),
        panel.background = element_rect(fill='transparent'),
        panel.grid=element_line(color='grey'),legend.title = element_text(size = 13))+ 
  labs(x="",y="")
dev.off()



setwd("E:/IMA/Pro/2019_075/fold1.5/Report/3.DEG/")
fpkm=read.table(file = "Ctrl_vs_XM1905_OE.differentially_expressed_genes.xls",header = T,row.names = 1,sep = "\t")
draw=fpkm[,!names(fpkm) %in% c("description","total_count","baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
data_plot<-t(scale(t(draw)))
data_plot[data_plot>1.5]<-1.5
data_plot[data_plot<(-1.5)]<-(-1.5)

pheatmap(data_plot)
cordraw=cor(draw)
pheatmap(cordraw)

draw=fpkm[,!names(fpkm) %in% c("norm")]
head(draw)


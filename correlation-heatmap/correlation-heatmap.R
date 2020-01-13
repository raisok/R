library(reshape2)
library(ggplot2)
library(RColorBrewer)
x <- read.table("/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/example/result.SE/Analysis_Report/BGI_result/Quantify/GeneExpression/CorrelationHeatmap/AllSamples.correlation.xls", sep = "	", head = T)
xx = as.matrix(x[,-1])
rownames(xx) = names(x)[-1]
xx = melt(xx)
names(xx)=c("Var1","Var2","pearson_value");
pdf("/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/example/result.SE/Analysis_Report/BGI_result/Quantify/GeneExpression/CorrelationHeatmap/AllSamples.CorrelationHeatmap.pdf",width=8,height=8)
ggplot(xx, aes(Var1, Var2, fill=pearson_value))+
 #geom_tile(width=0.8, height=0.8)+
  geom_tile(color='black')+
  geom_text(label=round(xx$pearson_value, 3))+
  scale_fill_gradient(low='#DEEBF7',high='#08519C')+
  theme(axis.text = element_text(angle=30, hjust=1,size=11,vjust=0,color='black'),
  panel.background = element_rect(fill='transparent'),
  panel.grid=element_line(color='grey'),legend.title = element_text(size = 13))+ 
  labs(x="",y="")
dev.off()

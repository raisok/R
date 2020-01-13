library(ggplot2)

setwd("E:\\IMA\\R\\barplot")
gene_fold ="cell_death_related_gene_list.xls"
out_pdf="cell_death_related_gene_list.pdf"

args  <- commandArgs(TRUE)

#gene_fold = args[1]
#out_pdf = args[2]


data=read.table(file=gene_fold,sep="\t",header=T,check.names = F)
data = data[order(data[,2], decreasing =T),]

data$Gene=factor(data$Gene,levels=data$Gene)

p = ggplot(data,aes(x=Gene,y=data$`log2(Foldchange)`))+
  geom_bar(stat = "identity",fill ="red",colour="black")+
  xlab("Gene")+ylab("Log2Folchange")+
  theme_classic()+
  theme(
        axis.text.x=element_text(size=6,color="black",angle = 30,hjust = 1)
  )+geom_hline(aes(yintercept=1), colour="black", linetype="dashed")

ggsave(out_pdf,p,height=5,width=6)



out_png = sub(".pdf",".png",out_pdf)
ggsave(out_png,p,height=5,width=6)
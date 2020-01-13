args = commandArgs()

library(ReactomePA)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)
library(stringr)

#setwd("E:/IMA/Rscript/Function/")
###
#args <- c(1,1,1,1,"Sham_vs_12-HETE.differentially_expressed_genes.xls","Sham_vs_12-HETE")
###

num<- function(data,n){
  if(nrow(data)>n) {
    data <- data[1:n,]
  } else {
    data <- data
  }
  return(data)
}

shorten_names <- function(x, n_word=4, n_char=40){
  shortname = x["Description"]
  if (length(strsplit(shortname, " ")[[1]]) > n_word || (nchar(shortname) > n_char))
  {
    if (nchar(shortname) > n_char) shortname <- substr(shortname, 1, n_char)
    y <- paste(paste(strsplit(shortname, " ")[[1]][1:min(length(strsplit(shortname," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(y)
  } 
  else
  {
    return (shortname)
  }
}



DEGs <- read.table(args[5],sep="\t",header=T,check.names=F,quote="")
entrezID2symbol <- AnnotationDbi::select(org.Mm.eg.db,keys=as.character(DEGs$gene),columns='ENTREZID',keytype='SYMBOL')

da <- enrichPathway(entrezID2symbol$ENTREZID,organism="mouse",pvalueCutoff=0.05,qvalueCutoff=1,readable=T)

write.table(da,file=paste(args[6],".reactome_enrichment.xls",sep=""),sep="\t",quote=F,row.names=F)

da<- da[order(da$pvalue),]
da <- num(da,15)

da$pretty_varname = apply(da,1,shorten_names)

bar <- ggplot(da,aes(y=Count,x=reorder(pretty_varname,-log10(pvalue)))) +
  geom_bar(aes(fill=-log10(pvalue)),stat="identity",width=0.8,position = position_dodge(0.9)) +
  scale_fill_gradientn(colours=c(rev(brewer.pal(11, "RdYlBu")[2:10]))) +
  coord_flip() +labs(color=expression(-log10(pvalue)),y="Gene Number",x="Pathway") +
  theme_bw() +
  theme(axis.text.y=element_text(size=6,color = "black"),
        axis.text.x=element_text(size=6,color = "black"),
        axis.title.y=element_text(size=6,color = "black"),
        axis.title.x=element_text(size=6,color = "black"),
        
        panel.grid=element_blank(),
        
        legend.key.size=unit(0.3,'cm'),
        legend.text = element_text(colour="black", size = 6),
        legend.title = element_text(colour="black", size=6)
        )

w_h <-function(data,data1){
  max1<- max(nchar(c(as.character(data1))))
  w<- max1*0.09
  if(nrow(data)>10) {
    h <- 8*nrow(data)/60
  } 
  else {
    h <- 1.6
  }
  return(c(w,h))
}


wh = w_h(da,da$pretty_varname)
w=wh[1]
h=wh[2]


ggsave(file=paste(args[6],".reactome_enrichment.pdf",sep=""),bar,width=w,height=h)
ggsave(file=paste(args[6],".reactome_enrichment.png",sep=""),bar,width=w,height=h)



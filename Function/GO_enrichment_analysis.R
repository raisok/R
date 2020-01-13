args = commandArgs()

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)

###
#setwd("E:/IMA/Rscript/Function/")
#args[5]="Sham_vs_12-HETE.differentially_expressed_genes.xls"
#args[6]="Sham_vs_12-HETE"
#args <- c(1,1,1,1,"Sham_vs_12-HETE.differentially_expressed_genes.xls","Sham_vs_12-HETE")
###


### change GO Term name length
shorten_names2 <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

shorten_names <- function(x, n_word=4, n_char=40){
  shortname = x["Description"]
  Go_id = x["ID"]
  if (length(strsplit(shortname, " ")[[1]]) > n_word || (nchar(shortname) > n_char))
  {
    if (nchar(shortname) > n_char) shortname <- substr(shortname, 1, n_char)
    y <- paste(paste(strsplit(shortname, " ")[[1]][1:min(length(strsplit(shortname," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
	y <- paste(Go_id,y,sep=" / ")
    return(y)
  } 
  else
  {
    return (shortname)
  }
}
### add by yueyao at 20191121

### 
set_width <- function(para1){
  #para<- para1[order(para1$p.adjust),][1:15,]
  Des <- para1$pretty_varname
  max1<- max(nchar(c(as.character(Des))))
  width<- max1*0.09
  if(nrow(para1)>10) {
   height <- 8*nrow(para1)/60
  } else {
	height <- 1.6

  }
  return(c(width,height))
}
### add by hxl at 20191121

num<- function(data,n){
  if(nrow(data)>n) {
    data <- data[1:n,]
  } else {
    data <- data
  }
  return(data)
}


###

draw_bar_plot <- function(goinput){
  
  bar <- ggplot(goinput,aes(y=Count,x=reorder(pretty_varname,-log10(pvalue)))) +
    geom_bar(stat="identity", width=0.7,aes(fill=-log10(pvalue))) + coord_flip() + 
    scale_fill_gradientn(colours=c(rev(brewer.pal(9, "RdYlBu")))) +
    labs(color=expression(-log10(pvalue)),x="GO term",y="Gene Number",title="Most Top 15 enrich GO Term")+
    theme_bw() +
    theme(
      title =element_text(size=6,face="plain",colour="black"),
      axis.title = element_text(size=6,face="plain",colour="black"),
      panel.grid=element_blank(),
      
      axis.line = element_line(size=0.25,colour = "black"),
      axis.text = element_text(color = "black",size = 6,face="plain"),
	  legend.key.size=unit(0.3,'cm'),
      legend.title=element_text(size=6,colour = "black",face="plain"),
      legend.text = element_text(colour = 'black',face = 'plain',size=6)
      )
    return (bar)
}

###


DEGs <- read.table(args[5],sep="\t",header=T,check.names=F,quote="")$gene

ont <- c("MF","BP","CC")
type <- c("molecular_function","biological_process","cellular_component")

ego_MF <- enrichGO(gene=as.vector(DEGs),OrgDb=org.Mm.eg.db,ont=ont[1],pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=1,keyType='SYMBOL')
write.table(ego_MF,file=paste(args[6],".",type[1],"_enrichment.xls",sep=""),sep="\t",quote=F,row.names=F)


ego_MF<- ego_MF[order(ego_MF$pvalue),]
ego_MF <- num(ego_MF,15)

ego_MF$pretty_varname = apply(ego_MF,1,shorten_names)
w <- set_width(ego_MF)

ego_MF<- draw_bar_plot(ego_MF)

ggsave(file=paste(args[6],".",type[1],"_enrichment.pdf",sep=""),ego_MF,width=w[1],height=w[2])
ggsave(file=paste(args[6],".",type[1],"_enrichment.png",sep=""),ego_MF,width=w[1],height=w[2])



ego_BP <- enrichGO(gene=as.vector(DEGs),OrgDb=org.Mm.eg.db,ont=ont[2],pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=1,keyType='SYMBOL')
write.table(ego_BP,file=paste(args[6],".",type[2],"_enrichment.xls",sep=""),sep="\t",quote=F,row.names=F)

ego_BP<- ego_BP[order(ego_BP$pvalue),]
ego_BP <- num(ego_BP,15)


ego_BP$pretty_varname = apply(ego_BP,1,shorten_names)
w <- set_width(ego_BP)
ego_BP<-draw_bar_plot(ego_BP)

ggsave(file=paste(args[6],".",type[2],"_enrichment.pdf",sep=""),ego_BP,width=w[1],height=w[2])
ggsave(file=paste(args[6],".",type[2],"_enrichment.png",sep=""),ego_BP,width=w[1],height=w[2])



ego_CC <- enrichGO(gene=as.vector(DEGs),OrgDb=org.Mm.eg.db,ont=ont[3],pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=1,keyType='SYMBOL')
write.table(ego_CC,file=paste(args[6],".",type[3],"_enrichment.xls",sep=""),sep="\t",quote=F,row.names=F)

ego_CC<- ego_CC[order(ego_CC$pvalue),]
ego_CC <- num(ego_CC,15)


ego_CC$pretty_varname = apply(ego_CC,1,shorten_names)
w <- set_width(ego_CC)
ego_CC<-draw_bar_plot(ego_CC)

ggsave(file=paste(args[6],".",type[3],"_enrichment.pdf",sep=""),ego_CC,width=w[1],height=w[2])
ggsave(file=paste(args[6],".",type[3],"_enrichment.png",sep=""),ego_CC,width=w[1],height=w[2])
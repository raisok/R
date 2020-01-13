del_human = T
gmtout = F
getko = F
outkopath = F
p = 0.05
fdr = 1
zs = 0
library(hash)
library(ggplot2)
library(RColorBrewer)

degfile <- "TAC_WT_vs_TAC_KO.differentially_expressed_genes.list"

species <- "mmu"

groupname <- "TAC_WT_vs_TAC_KO"


time <- gsub("-","",Sys.Date())
species_ko <- tolower(species)
species_list <- data.frame(matrix(c("mouse","mmu","human","hsa"),ncol=2,byrow=T),stringsAsFactors=F)
names(species_list) <- c("Species","Species_kegg")
if(species_ko%in%species_list[,1]){
  species_ko <- species_list[grep(species_ko,species_list[,1]),2]
}

if(getko==T){
  ko_gene <- kegg_ko[grep("^D",kegg_ko[,1]),1:2]
  names(ko_gene) <- c("gene","ko")
  ko_gene[,1] <- sub(";.*","",ko_gene[,1])
  ko_gene[,1] <- sub("^D.* ","",ko_gene[,1])
  ko_gene[,2] <- sub(" +.*","",ko_gene[,2])
  ko_gene <- ko_gene[!duplicated(ko_gene[,1]),]
  write.table(ko_gene,paste(species_ko,"_ko_gene_",time,".xls",sep=""),sep="\t",quote=F,row.names=F)
}


ko <- data.frame(kegg_ko[grep("^A.+|^B.+|^C.+|^D.+",kegg_ko[,1]),1:2],stringsAsFactors=F)
ko$level1 <- NA
ko$level2 <- NA
ko$level3 <- NA


level1_list <- grep("^A",ko[,1])
n <- length(level1_list)
for(i in 1:(n-1)){
  ko$level1[level1_list[i]:(level1_list[i+1]-1)] <- ko[level1_list[i],1]
}
ko$level1[level1_list[n]:nrow(ko)] <- ko[level1_list[n],1]


## level2
level2_list <- grep("^B",ko[,1])
n <- length(level2_list)
for(i in 1:(n-1)){
  ko$level2[level2_list[i]:(level2_list[i+1]-1)] <- ko[level2_list[i],1]
}
ko$level2[level2_list[n]:nrow(ko)] <- ko[level2_list[n],1]

## level3
level3_list <- grep("^C",ko[,1])
n <- length(level3_list)
for(i in 1:(n-1)){
  ko$level3[level3_list[i]:(level3_list[i+1]-1)] <- ko[level3_list[i],1]
}
ko$level3[level3_list[n]:nrow(ko)] <- ko[level3_list[n],1]


ko_path <- ko[grep("^D",ko[,1]),]
ko_path <- ko_path[grep(paste("PATH:",species_ko,sep=""),ko_path[,5]),]
ko_path[,1] <- sub(";.*","",ko_path[,1])
ko_path[,1] <- sub("^D +[0-9].*? ","",ko_path[,1])
ko_path[,2] <- sub(" +.*","",ko_path[,2])
ko_path[,3] <- sub("^A.*? ","",ko_path[,3])
ko_path[,4] <- sub("^B +[0-9].*? ","",ko_path[,4])
ko_path[,5] <- sub("^C +[0-9].*? ","",ko_path[,5])
ko_path[,5] <- sub("/"," ",ko_path[,5])
ko_path[,5] <- sub(" +"," ",ko_path[,5])
ko_path$pathway <- sub(".*:","",ko_path[,5])
ko_path$pathway <- sub("]","",ko_path$pathway)
ko_path[,5] <- sub(" \\[.*\\]","",ko_path[,5])
names(ko_path)[1:2] <- c("Gene","ko")
if(length(grep("^$",ko_path$Gene))>0){
  ko_path<-ko_path[-grep("^$",ko_path$Gene),]}
## ??????????????????????????????ko_path
if(del_human==T){
  ko_path <- ko_path[!(ko_path[,3]=="Human Diseases"),]
  if(outkopath==T){
    write.table(ko_path,"ko_path_symbol_nohuman.xls",sep="\t",quote=F,row.names=F)
  }
} else {
  if(outkopath==T){write.table(ko_path,"ko_path_symbol.xls",sep="\t",quote=F,row.names=F)
  }
}

pre_hash <- ko_path[,c(1,5)]
hash_gene <- unique(pre_hash[,1])
#hash_gene<- hash_gene[-c(grep('^$',hash_gene))]
hash_gene_path <- hash()
for(i in 1:length(hash_gene)){
  .set(hash_gene_path,keys=hash_gene[i],values=list(pre_hash[pre_hash[,1]==hash_gene[i],2]))
}

hash_path <- unique(pre_hash[,2])
hash_path_gene <- hash()
for(i in 1:length(hash_path)){
  .set(hash_path_gene,keys=hash_path[i],values=list(unique(pre_hash[pre_hash[,2]==hash_path[i],1])))
}


if(gmtout==T){
  gmt <- data.frame(keys(hash_path_gene),stringsAsFactors=F)
  gmt$gene <- NA
  for(i in 1:nrow(gmt)){
    gmt[i,]$gene <- paste(unlist(values(hash_path_gene,gmt[i,1])),collapse=";",sep="")
  }
  gmt[,1] <- gsub("-","_",gmt[,1])
  gmt$out <- paste(gmt[,1],gmt$gene,sep=";")
  gmt_out <- data.frame(gmt$out,stringsAsFactors=F)
  gmt_out <- gsub(";","\t",gmt_out[,1])
  write.table(gmt_out,paste(species_ko,"_",time,".gmt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}


gene_annotation <- unique(ko_path$Gene)
deg <- read.table(degfile,header=T,sep="\t",quote="",stringsAsFactors=F)
deg_annotation <- intersect(deg[,1],gene_annotation)
path_counts <- unique(unlist(values(hash_gene_path,deg_annotation)))
path_counts <- data.frame(path_counts,stringsAsFactors=F)
kegg_out <- ko_path[,c(3:6)]
kegg_out <- kegg_out[!duplicated(kegg_out[,3]),]
kegg_out <- merge(kegg_out,path_counts,by.x="level3",by.y="path_counts")
kegg_out <- kegg_out[,c(4,2,3,1)]
kegg_out$pathway_degs_num <- 0
degs_num_inpathway <- length(deg_annotation)
kegg_out$degs_num_inpathway <- degs_num_inpathway
kegg_out$pathway_gs_num_genome <- 0
gs_num_inpathway_genome <- length(gene_annotation)
kegg_out$gs_num_inpathway_genome <- gs_num_inpathway_genome
kegg_out$pvalue <- 1
kegg_out$FDR <- 1
kegg_out$zscore <- 0
kegg_out$pathway_degs <- NA
kegg_out$pathway_nondegs <- NA


for(i in 1:nrow(kegg_out)){
  path_gene <- unlist(values(hash_path_gene,kegg_out[i,]$level3))
  deg_inpath <- intersect(deg_annotation,path_gene)
  kegg_out[i,]$pathway_degs_num <- length(deg_inpath)
  kegg_out[i,]$pathway_gs_num_genome <- length(path_gene)
  kegg_out[i,]$pvalue <- phyper(kegg_out[i,]$pathway_degs_num-1,kegg_out[i,]$pathway_gs_num_genome,gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome,degs_num_inpathway,lower.tail=FALSE)
  kegg_out[i,]$zscore <- (kegg_out[i,]$pathway_degs_num-degs_num_inpathway*kegg_out[i,]$pathway_gs_num_genome/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome)) / (sqrt(degs_num_inpathway*(kegg_out[i,]$pathway_gs_num_genome/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome))*(1-kegg_out[i,]$pathway_gs_num_genome/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome))*(1-(degs_num_inpathway-1)/(kegg_out[i,]$pathway_gs_num_genome+gs_num_inpathway_genome-kegg_out[i,]$pathway_gs_num_genome-1))))
  kegg_out[i,]$pathway_degs <- paste(deg_inpath,collapse=";")
  kegg_out[i,]$pathway_nondegs <- paste(path_gene,collapse=";")
}
kegg_out$FDR <- p.adjust(kegg_out$pvalue,method="fdr")
kegg_out <- kegg_out[order(kegg_out$pvalue),]
kegg_out_filter <- kegg_out[kegg_out$pvalue<p,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$FDR<fdr,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$zscore>zs,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$pathway_degs_num>=3,]
write.table(kegg_out_filter,paste(groupname,".kegg_enrichment.xls",sep=""),sep="\t",row.names=F,quote=F)
write.table(kegg_out,paste("unfiltered.",groupname,".kegg_enrichment.xls",sep=""),sep="\t",row.names=F,quote=F)


pathway <- kegg_out_filter
if(nrow(pathway)>15) {
  data <- pathway[1:15,]
} else {
  data <- pathway
}



if(nrow(data)>10) {
  h <- 8*nrow(data)/60
} else if(nrow(data)<=5){
  h <- 1.2
}else{
  h <- 1.6
}

shorten_names <- function(x, n_word=4, n_char=40){
  shortname = x["level3"]
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

data$pretty_varname = apply(data,1,shorten_names)
max1<- max(nchar(c(as.character(data$pretty_varname))))
w<- max1*0.09




bar <- ggplot(data,aes(y=pathway_degs_num,x=reorder(pretty_varname,-log10(pvalue)))) +
  geom_bar(aes(fill=-log10(pvalue)),stat="identity",width=0.8,position = position_dodge(0.9)) +
  scale_fill_gradientn(colours=c(rev(brewer.pal(11, "RdYlBu")[2:10]))) +
  #scale_fill_gradientn(colours=c(rev(brewer.pal(11, "Spectral")[1:5]))) +
  coord_flip() +labs(color=expression(-log10(pvalue)),y="Gene Number",x="Pathway") +
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
ggsave(paste(groupname,".kegg_enrichment_bar.png",sep=""),bar,width=w,height=h)
ggsave(paste(groupname,".kegg_enrichment_bar.pdf",sep=""),bar,width=w,height=h)


## ?????????
colors <- length(data$level2[!duplicated(data$level2)])
category <-  ggplot(data,aes(x=-log10(pvalue),y=reorder(pretty_varname,-log10(pvalue)))) +
  geom_segment(aes(yend=reorder(pretty_varname,-log10(pvalue))),xend=0,colour="grey50") +
  geom_point(aes(size=pathway_degs_num,colour=level2)) +
  guides(colour=FALSE) +scale_size_area(max_size=5) +
  facet_grid(level2~.,scales="free_y",space="free_y") +
  labs(color=expression(-log10(pvalue)),size="Gene number",x="-log10(P-value)",y="Pathway") +
  theme_bw()+
  theme(
        strip.text.y=element_text(angle=0,size=6),
        strip.background = element_rect(fill=brewer.pal(11, "RdYlBu")[9]),
        
        legend.position="top",
        title =element_text(size=6,face="plain",colour="black"),
        axis.title = element_text(size=6,face="plain",colour="black"),
        panel.grid=element_blank(),
        
        axis.text = element_text(color = "black",size = 6,face="plain"),
        legend.title=element_text(size=6,colour = "black",face="plain"),
        legend.text = element_text(colour = 'black',face = 'plain',size=6)
  )

w=w*1.5
#h=h*2

if(nrow(data)>10) {
  h=h*2
} else if(nrow(data)<=5){
  h <- h*1.5
}else{
  h <- h*1.8
}

ggsave(paste(groupname,".kegg_enrichment_category.png",sep=""),category,width=w,height=h)
ggsave(paste(groupname,".kegg_enrichment_category.pdf",sep=""),category,width=w,height=h)



library("gplots")
library("pheatmap")
data.inter<-read.table("/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/example/result.SE/process/Cluster_Pheatmap/HBRR1-VS-UHRR1.PossionDis-HBRR2-VS-UHRR2.PossionDis.inter.tmp",head=TRUE)
len<-length(data.inter)
row <- data.inter$V1
colnames(data.inter) <- c("geneID","HBRR1-VS-UHRR1.PossionDis  ","HBRR2-VS-UHRR2.PossionDis  ")
rownames(data.inter) <- as.vector(row) 
data.inter<-data.inter[,2:len]
mycolors <- colorRampPalette(c("#377EB8","white","#E41A1C"))(1000)
#mycolors <- colorRampPalette(c("blue","white","red"))(1000)
pdf("/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/example/result.SE/process/Cluster_Pheatmap/HBRR1-VS-UHRR1.PossionDis-HBRR2-VS-UHRR2.PossionDis.inter.pdf",width=10,height=10)
obj <- pheatmap(
        data.inter,
        show_rownames=FALSE,
        show_colnames=TRUE,
        col=mycolors,
        cluster_rows=TRUE,
        cluster_cols=TRUE,
        legend=TRUE,
        fontsize=15,
        main="Hierarchical Clustering of DEGs(Inter)",
        display_numbers=FALSE,
        cellwidth = 200,
        border_color=NA
)
order <- obj$tree_row$order
write.table(order,file="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/example/result.SE/process/Cluster_Pheatmap/HBRR1-VS-UHRR1.PossionDis-HBRR2-VS-UHRR2.PossionDis.inter.cluster.order",quote=FALSE)
dev.off()

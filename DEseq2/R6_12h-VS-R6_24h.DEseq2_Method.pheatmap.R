library("pheatmap")
data <- read.table("/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.DEseq2_Method.prepare4pheatmap-plot",head=T,row.names=1)
aa <- head(data,1)
a <- t(aa)
sample <- row.names(a)[-1]
height <- 4
if (length(data) >10){
width <- 2 + length(data)*0.2
}else{
width <- 4
}
group <- data.frame(Var1 = c("R6_12h","R6_12h","R6_12h","R6_24h","R6_24h","R6_24h")) 
rownames(group) = paste(sample)
colnames(group) = paste("Group")

Up_Down <- data.frame(data[,1]) 
rownames(Up_Down) = paste(row.names(data))
colnames(Up_Down) = paste("Up_Down")

################################################
###you can custom colors,or use default in R.###
################################################
ann_colors = list (Group = c("R6_12h"="#BEBADA","R6_24h"="#FDC086"),Up_Down=c("Up"="#E41A1C","Down"="#377EB8"))

pheatmap(
	data[,-1],
        annotation_col = group,
        annotation_row = Up_Down,
        annotation_names_col =F,
	annotation_names_row =F,
        annotation_colors = ann_colors,##group colors
        fontsize_col=9,
        fontsize_row=2.7,
        border_color = "transparent",
        color = colorRampPalette(c("#FFFFD9","yellow2","orange","orangered2","gray10"))(200), #####colors and gradient for pheatmap 
	main="Pheatmap for R6_12h&R6_24h",
	show_rownames = F,
	file = '/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.DEseq2_Method.pheatmap-plot.pdf',
	width=width,
	height=height
)
dev.off()

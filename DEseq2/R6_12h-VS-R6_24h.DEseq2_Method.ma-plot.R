pdf("/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.DEseq2_Method.MA-plot.pdf",width=12,height=10,pointsize=15)
data <- read.table("/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.DEseq2_Method.prepare4MA-plot", sep = "\t", skip = 1)
no_degs_x <- data$V2[ grepl("\\*",data$V3,perl=TRUE) ]
no_degs_y <- data$V1[ grepl("\\*",data$V3,perl=TRUE) ]
degs_x <- data$V2[ grepl("[^\\*]",data$V3,perl=TRUE) ]
degs_y <- data$V1[ grepl("[^\\*]",data$V3,perl=TRUE) ]

no_degs_x <- 0.5 * log2(no_degs_x)
degs_x <- 0.5 * log2(degs_x)

xmin <- min(c(no_degs_x,degs_x)[grepl("\\d+",c(no_degs_x,degs_x),perl=TRUE)])
xmax <- max(c(no_degs_x,degs_x)[grepl("\\d+",c(no_degs_x,degs_x),perl=TRUE)])
ymin <- min(c(no_degs_y,degs_y)[grepl("\\d+",c(no_degs_y,degs_y),perl=TRUE)])
ymax <- max(c(no_degs_y,degs_y)[grepl("\\d+",c(no_degs_y,degs_y),perl=TRUE)])
aa <- factor(data$V3)
number <- as.data.frame(table(aa))
par(mar=c(6,6,6,14)+0.1,xpd=TRUE)
plot(0,0,pch="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),ylab="M",xlab="A",cex.lab = 1.6,main="MA plot of R6_12h-VS-R6_24h.DEseq2_Method",cex.main=1.5)
for (i in 1 : length(no_degs_x)) {
	points(no_degs_x[i], no_degs_y[i], col = "grey65", pch = 20, cex = 0.5)
}
for (j in 1 : length(degs_x)) {
	if (degs_y[j] > 0) {
		points(degs_x[j], degs_y[j], col = "#E41A1C", pch = 20, cex = 0.5)
	}else {
		points(degs_x[j], degs_y[j], col = "#377EB8", pch = 20, cex = 0.5)
	}
}
lines(c(xmin,xmax),c(0,0),col="black",lty=2)
legend(xmax+(xmax-xmin)/12,(ymax+ymin)/2+(ymax-ymin)/5-(ymax-ymin)/15, legend=expression(atop(italic("log")[2]*FoldChange >=1.~ italic(",Padj") <=0.05)), pch="",bty = "n",pt.cex=1,cex=0.7)
legend(xmax+(xmax-xmin)/12,(ymax+ymin)/2-(ymax-ymin)/15, legend=expression(atop(italic("log")[2]*FoldChange <=-1.~ italic(",Padj") <=0.05)), pch="",bty = "n", pt.cex=1.8, cex=0.7)
legend(xmax+(xmax-xmin)/12,(ymin+ymax)/2-(ymax-ymin)/5-(ymax-ymin)/15, legend=expression(atop(italic("abs(log")[2]*"FoldChange)" <1.~ italic(" or Padj") >0.05)), pch="",bty = "n", pt.cex=1.8, cex=0.7)
legend(xmax+(xmax-xmin)/15,(ymax+ymin)/2+(ymax-ymin)/5, legend=paste("Up:",number[3,2]),col="#E41A1C",pch=20,bty = "n", pt.cex=2, cex=1.3)
legend(xmax+(xmax-xmin)/15,(ymax+ymin)/2, legend=paste("Down:",number[2,2]),col="#377EB8",pch=20,bty = "n", pt.cex=2, cex=1.3)
legend(xmax+(xmax-xmin)/15,(ymin+ymax)/2-(ymax-ymin)/5, legend=paste("no-DEGs:",number[1,2]),col="grey65",pch=20,bty = "n", pt.cex=2, cex=1.3)
dev.off()

library(plotrix)
pdf(file="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/example/result.SE/Analysis_Report/BGI_result/Quantify/GeneExpression/VennDiagram/HBRR1-HBRR2-HBRR3/HBRR1-HBRR2-HBRR3.venn.pdf",width=8.8,height=8)
color <- c("#E41A1C","#377EB8","#FDB462")
color_transparent <- adjustcolor(color, alpha.f = 0.2) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################
p_x   <- c(0,13,-13,  10.4,-10.4,0, 0)
p_y   <- c(13,-9,-9,  4,4,-13.5, -1)
p_lab <- c(733,617,574,420,375,1378,15039)

title_x <- c(0,17,-17)
title_y <- c(19.8,-18,-17)
title_lab <- c("HBRR1\n(16567)","HBRR2\n(17454)","HBRR3\n(17366)")
########################################
par(mar=c(7,10,7,8)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",,xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(0,4,-4), c(3.02,-3.912,-3.912), c(14,14,14), c(14,14,14),border=color_transparent1,
angle=c(60,120,0), lty=1,col = color_transparent,lwd = 2)
text (p_x,p_y,p_lab,cex=1,col="grey20")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
dev.off()

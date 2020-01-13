library("DESeq2")
countdata <- read.table("/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.deseq2.list",skip=1)
len <- length(countdata)
rownames(countdata) <- countdata[,1]
countdata <- countdata[,2:len]
type <- c("Control","Control","Control","Treat","Treat","Treat")
coldata <- data.frame(type)
rownames(coldata) <- c("R6_12h_1","R6_12h_2","R6_12h_3","R6_24h_1","R6_24h_2","R6_24h_3")
colnames(countdata) <- c("R6_12h_1","R6_12h_2","R6_12h_3","R6_24h_1","R6_24h_2","R6_24h_3")

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~ type)
dds <- DESeq(dds,quiet=TRUE)
sizefactor <- sizeFactors(dds)
result <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE, pAdjustMethod="BH")
write.table(result, file="/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.deseq2.output", quote=FALSE, sep="\t")
write.table(sizefactor, file="/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneDiffExp_Allin/DEseq2/R6_12h-VS-R6_24h.deseq2.sizefactors", quote=FALSE, sep="\t")

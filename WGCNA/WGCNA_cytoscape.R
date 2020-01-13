library(WGCNA)
library(flashClust)
library(iterators)
options(stringsAsFactors = FALSE)
options(expressions = 50000)##
enableWGCNAThreads()

femData = read.table('/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneCoExpression_WGCNA/geneExpMatrix.xls',header=TRUE)
datExpr = as.data.frame(t(femData[, -c(1:1)]))

# gene names
names(datExpr) = femData$GeneID

#sample names
rownames(datExpr) = names(femData)[-c(1:1)]

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

gene.names=names(datExpr)

# Choosing a soft-threshold to fit a scale-free topology to the network
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = 'unsigned')
Rsquare<-(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
softPower <- sft$fitIndices[which(Rsquare==max(Rsquare)),1]

TOM=TOMsimilarityFromExpr(datExpr,networkType = 'unsigned', TOMType = 'unsigned', power = softPower)

# Export the network into edge and node list files Cytoscape can read
probes = names(datExpr)
altNodeNames = NA
#=======
#filter data for generate result 
#adjMat = as.matrix(TOM)
#adjMat[is.na(adjMat)] = 0
#nRow = nrow(adjMat)
#checkAdjMat(adjMat, min = -1, max = 1)
#adjDst = as.dist(adjMat)
#edges = abs(adjDst) > threshold
#while (length(adjDst[abs(adjDst)>threshold])<num){
#threshold = threshold - 0.01
#}
#==filter by nodes
threshold = 0.9
#0.9
cyt = exportNetworkToCytoscape(TOM, weighted = TRUE, threshold = threshold, nodeNames = probes, altNodeNames = altNodeNames)
if (length(unique(c(cyt$edgeData$fromNode,cyt$edgeData$toNode))) >50){
0.9
}else{
num = 1000
if (threshold>=0){
while (length(unique(c(cyt$edgeData$fromNode,cyt$edgeData$toNode))) <num ){
cyt = exportNetworkToCytoscape(TOM, weighted = TRUE, threshold = threshold, nodeNames = probes, altNodeNames = altNodeNames)
threshold = threshold - 0.01
}
}else{
threshold = 0
}
}
threshold
#=======
#cyt = exportNetworkToCytoscape(TOM, edgeFile = paste('/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneCoExpression_WGCNA/CytoscapeInput-edges', '.txt', sep=''), weighted = TRUE, threshold = 0.9, nodeNames = probes, altNodeNames = altNodeNames);
cyt = exportNetworkToCytoscape(TOM, edgeFile = paste('/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneCoExpression_WGCNA/CytoscapeInput-edges', '.txt', sep=''), weighted = TRUE, threshold = threshold, nodeNames = probes, altNodeNames = altNodeNames);


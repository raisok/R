library(WGCNA)
library(flashClust)
library(iterators)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

femData = read.table('/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneCoExpression_WGCNA/geneExpMatrix.InCytoscape.txt',header=TRUE)
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

colnames(TOM) =rownames(TOM) =gene.names
dissTOM=1-TOM

# Module detection
geneTree = flashClust(as.dist(dissTOM),method='average')

# Set the minimum module size
minModuleSize = 20

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method='hybrid', deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unasunsigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unasunsigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#discard the unasunsigned genes, and focus on the rest
restGenes= (dynamicColors != 'grey')
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)
hier1=flashClust(as.dist(diss1), method='average')

#set the diagonal of the dissimilarity to NA
diag(diss1) = NA

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
pdf('/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneCoExpression_WGCNA/Modules_heatmap.pdf')
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))
dev.off()

# Extract modules
module_colors= setdiff(unique(dynamicColors), 'grey')

for (color in module_colors){
    module=gene.names[which(dynamicColors==color)]
    write.table(module, paste('/ldfssz1/ST_BIGDATA/USER/yueyao/12.Pro/test/result/process/GeneCoExpression_WGCNA/Modules_gene_',color, '.txt',sep=''), sep= , row.names=FALSE, col.names=FALSE,quote=FALSE)
}



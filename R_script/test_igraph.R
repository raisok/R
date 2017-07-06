#!/usr/bin/R
library("igraph")
setwd("c:/Users/yueyao/Desktop/")
data = read.table(file = "test.txt", header = T)
g <- graph.data.frame(data);

color_index = V(g)$name %in% data[,1]
color = character(0)
color[color_index]="#FFFF00FF"
color[!color_index]="#FF0FFFFF"
V(g)$color <- color;

size_index = V(g)$name %in% data[,1]
V(g)$size=4;
V(g)[size_index]$size = 12

q=0.999;
V(g)$bte=betweenness(g,directed=FALSE);
or=order(V(g)$bte);
ind=round(length(or)*q);
btemax=V(g)$bte[or[ind]];
V(g)[bte>btemax]$size=16;

label_index = V(g)$name %in% data[,1]
V(g)$label=NA;
V(g)[label_index]$label= V(g)[label_index]$name;

#pdf("/home/Project/KC2017-A05/miRNA/11.kegg/UACC812.3440_UACC812/sig_Down_miRNAs_target_gene_network.pdf")
plot(g,layout=layout.fruchterman.reingold, vertex.size=V(g)$size,
     vertex.color= color, vertex.label=V(g)$label,
     vertex.label.cex=V(g)$cex,
     edge.color=grey(0.5),edge.arrow.mode="-")
#dev.off()

#png(filename ="/home/Project/KC2017-A05/miRNA/11.kegg/UACC812.3440_UACC812/sig_Down_miRNAs_target_gene_network.png", width = 800, height = 800, pointsize = 20)
plot(g,layout=layout.fruchterman.reingold, vertex.size=V(g)$size,
     vertex.color= color, vertex.label=V(g)$label,
     vertex.label.cex=V(g)$cex,
     edge.color=grey(0.5),edge.arrow.mode="-")
#dev.off()

#!/usr/bin/R
library("igraph")

data = read.table("/home/Project/KC2017-A05/miRNA/11.kegg/UACC812.3440_UACC812/sig_Down_miRNAs_target_gene.txt", header = T)
g <- graph.data.frame(data);

color_index = V(g)$name %in% data[,1]
color = character(0)
color[color_index]="#FFFF00FF"
color[!color_index]="#FF0FFFFF"
V(g)$color <- color;

size_index = V(g)$name %in% data[,1]
V(g)$size=4;
V(g)[size_index]$size = 12

q=0.999;
V(g)$bte=betweenness(g,directed=FALSE);
or=order(V(g)$bte);
ind=round(length(or)*q);
btemax=V(g)$bte[or[ind]];
V(g)[bte>btemax]$size=16;

label_index = V(g)$name %in% data[,1]
V(g)$label=NA;
V(g)[label_index]$label= V(g)[label_index]$name;

pdf("/home/Project/KC2017-A05/miRNA/11.kegg/UACC812.3440_UACC812/sig_Down_miRNAs_target_gene_network.pdf")
plot(g,layout=layout.fruchterman.reingold, vertex.size=V(g)$size,
        vertex.color= color, vertex.label=V(g)$label,
        vertex.label.cex=V(g)$cex,
        edge.color=grey(0.5),edge.arrow.mode="-")
dev.off()

png(filename ="/home/Project/KC2017-A05/miRNA/11.kegg/UACC812.3440_UACC812/sig_Down_miRNAs_target_gene_network.png", width = 800, height = 800, pointsize = 20)
plot(g,layout=layout.fruchterman.reingold, vertex.size=V(g)$size,
        vertex.color= color, vertex.label=V(g)$label,
        vertex.label.cex=V(g)$cex,
        edge.color=grey(0.5),edge.arrow.mode="-")
dev.off()
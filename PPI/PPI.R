#args  <- commandArgs(TRUE)
#input <- args[1]
#input1 <- args[2]
#out   <- args[3]

input <- "R6_12h-VS-R6_24h.DEseq2_Method.network.relation.xls.100"
input1 <- "R6_12h-VS-R6_24h.DEseq2_Method.network.relation.xls.100.fd"
out   <- "R6_12h-VS-R6_24h.DEseq2_Method.network.pdf"

options(stringsAsFactors=F)
library(igraph)
t   = read.table(input,header=F)
a = read.table(input1,header=F)
pdf(out,width=12, height=12)
t1  = unique(t[,c(1,2)])
d = data.frame(p1=t1[,1],p2=t1[,2],weight=0.1)
g = graph.data.frame(d, directed = FALSE)
v = c(unique(d[,1]),unique(setdiff(d[,2],d[,1])))
#a$V1[order(V(g)$name)]
E(g)$color <- "black"
Vla = ifelse(vcount(g)<100,0.5,ifelse(vcount(g)>500,0.1,0.3))
V(g)$lable.cex = rep(Vla,vcount(g))
shape = rep("circle",vcount(g))
if(vcount(g)<60){
  size = rep(15,vcount(g))
}

if(vcount(g)<300 & vcount(g) >=60){
  size = rep(6,vcount(g))
}
if(vcount(g) < 800 & vcount(g) >= 300){
  size = rep(2,vcount(g))
}
if(vcount(g) >= 800){
  size = rep(2,vcount(g))
}
a<-a[match(V(g)$name,a[,1]),]
a$V4[a$V3<0]="#377EB8"
a$V4[a$V3>0]="#E41A1C"
#d[match(c[,2],d[,2]),] 
#V(g)$name
#size1=size/max(a$V2)/0.8*a$V2
a$V5<-log(a$V2)*1.8+6
#round(log(a$V2)*2+6)
#size1=log(a$V2)*2+6
V(g)$size = a$V5
V(g)$shape = shape
col=a$V4

#col = rep(rgb(127/255,196/255,228/255,0.5),vcount(g))
V(g)$color = col
vertex.frame.color = col
#plot(g,vertex.size=V(g)$size,vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,vertex.label.color="gray40",vertex.frame.color=vertex.frame.color,edge.color="gray85",layout=layout_in_circle,vertex.label.cex=.7,vertex.label.font=2)
plot(g,vertex.size=V(g)$size,vertex.label.cex=V(g)$lable.cex,
     edge.width = E(g)$weight,edge.arrow.size=0,vertex.label.color="black",
     vertex.frame.color=vertex.frame.color,edge.color="black",
     layout=layout_with_dh ,vertex.label.cex=.5,vertex.label.font=1,
     vertex.label.family="Times")
dev.off()
